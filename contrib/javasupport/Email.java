////////////////////////////////////////////////////////
//
// Simple Email Class
//
// SMPT mail sender to email a message to port 25 (the standard SMPT socket)
// Adheres to  TCP/SMTP specs (RFC 821)
//
// Written in Java 1.1 by Chris Parkinson, 1997
//
////////////////////////////////////////////////////////


import java.awt.event.*;
import java.io.*;
import java.net.*;
import java.util.Vector;

public class Email {
	
	//ID's for action messages
	public static final int MESSAGE = 0;
	public static final int ERROR = 1;

	//Mail server and port information
	private int mailPort;			//SMPT PORT - should usually be 25
	private String mailServer;  	//Host name or IP address

	//Individual mail information
	private String mailFrom;   		//Who is sending this message
	private String[] mailTo;   		//Who are the recipients
	private String mailSubject; 	//The subject line of the message
	private String mailMessage;		//The actual mail message

	//Transient variables
	private Socket socket;                  //Socket to talk to host
	private BufferedReader dataInput;  		//Receive from sockets
	private PrintWriter dataOutput; 		//Send to socket
	private String heloHost;                //The server identity
	private ActionListener actionListener;  //Who is monitoring us

 
////////////////////////////////////////////////////////
//
// Constructor
//
////////////////////////////////////////////////////////
/**
 * Construct a simple emailer
 */
 public Email() {
 	setMailPort(25);   //Set up some defaults
	setMailServer("pnl.gov");  
  }

////////////////////////////////////////////////////////
//
// Set the mail port number
//
////////////////////////////////////////////////////////
/**
 * Set the mail port number
 */
 public void setMailPort(int portNumber) {
 	mailPort = portNumber;
}

////////////////////////////////////////////////////////
//
// Set the mail server (a string of name or number)
//
////////////////////////////////////////////////////////
/**
 * Set the mail server
 */
 public void setMailServer(String server) {
 	mailServer = new String(server);
}

////////////////////////////////////////////////////////
//
// Set the from email address
//
////////////////////////////////////////////////////////
/**
 * Set the mail sender
 */
 public void setMailFrom(String from) {
 	mailFrom = new String(from);
}

////////////////////////////////////////////////////////
//
// Set the to email address(s)
//
////////////////////////////////////////////////////////
/**
 * Set the 'to' email address(s)
 */
 public void setMailTo(String to) {
 	mailTo = new String[1];
	mailTo[0] = new String(to);
}
   
/**
 * Set the 'to' email address(s)
 */
 public void setMailTo(String to[]) {
	for (int i =0; i<to.length; i++) addMailTo(to[i]);
}

/**
 * Add a 'to' address
 */
 public void addMailTo(String to) {
   	String[] newList = new String[mailTo.length +1];
	for (int i =0; i<mailTo.length; i++) newList[i+1] = mailTo[i];
	newList[0]= new String(to);
	mailTo = newList;
}

////////////////////////////////////////////////////////
//
// Set the email subject line
//
////////////////////////////////////////////////////////
/**
 * Set the email subject line
 */
 public void setMailSubject(String subject) {
 	mailSubject = new String(subject);
}


////////////////////////////////////////////////////////
//
// Set the mail message
//
////////////////////////////////////////////////////////
/**
 * Set the mail message
 */
 public void setMailMessage(String message) {
 	mailMessage = new String(message);
}

/**
 * Add text to the mail message
 */
 public void addMailMessage(String message) {
 	mailMessage = mailMessage + "\r\n" + new String(message);
}


////////////////////////////////////////////////////////
//
// Set action listener
//
////////////////////////////////////////////////////////
/**
 * Set a listener to monitor what is going on
 */
 public void setActionListener(ActionListener listener) {
 	actionListener = listener;
}



////////////////////////////////////////////////////////
//
// Send the Message, putting the results so far in the 
// text panel
//
////////////////////////////////////////////////////////
/**
 * Send the message
 */
public boolean sendMail() {
	//Check the validity of our data
	if (isDataValid()==false) return false;
  
	sendActionMessage("Connecting to server (" + mailServer+ ")...");
	if (openConnection()==false) return false; 	   //Initiate sockets & streams

	sendActionMessage("Waiting for response from server (" + mailServer + ")...");
	if (greetServer()==false)return false;		   //Start conversation with server
  
  	if (sendHeader()==false)return false;		   //Send the mail headers
	if (sendData()==false)return false;	 		   //Send the data and quit commands
 
	closeConnection();							   //Finish up
	sendActionMessage("Connection closed");
	sendActionMessage("Mail sent !");
	return true;
}


////////////////////////////////////////////////////////
//
// Send the Headers to the message
//
////////////////////////////////////////////////////////
/**
 * Send the headers
 */
private boolean sendHeader() {
	//Say Helo back to server
	if (sendToServer("HELO", heloHost, 250,true)==false)return false;
   
	//Mail From
    if (sendToServer("MAIL FROM:",  mailFrom, 250,true)==false)return false;

	//Mail To (Can have multiple recipients here)
	for (int i =0; i<mailTo.length; i++)
 	   	if (sendToServer("RCPT TO:",  mailTo[i], 250,true)==false)return false;
	return true;
}

					 
////////////////////////////////////////////////////////
//
// Send the actual message data
//
////////////////////////////////////////////////////////
/**
 * Send the message data
 */
private boolean sendData() {
	//Tell server we are sending data over
	if (sendToServer("DATA",  null, 354,true)==false)return false;

	//Protocol says that first should be the subject line + blank line
	String[] lines = validateMailData(mailSubject);
	sendActionMessage("Sending data to server...");
	if (sendToServer("Subject:",  lines[0], 0,false)==false)return false;
	if (sendToServer("", null, 0,false)==false)return false;      //Blank line 
 
	//Now the main body of the message (Note any single '.' should be changed to '..')
   	lines = validateMailData(mailMessage);
	for (int i =0; i< lines.length; i++)  {
		sendActionMessage("Sending data to server...");
 		if (sendToServer(lines[i], null, 0,false)==false)return false;
  	}	 

	//End DATA section with '\r\n.\r\n'
	if (sendToServer("\r\n.", null, 250,false)==false)return false;
	sendActionMessage("Sent data");

	//Finally send a QUIT message to close the conversation
   	return sendToServer("QUIT", null, 221,true);
}
	   

////////////////////////////////////////////////////////
//
// Greet the server
// The server should send a '220' service ready message
//
////////////////////////////////////////////////////////
/**
 * Greet the Server
 */
 private boolean greetServer() {
 	String input = getServerLine(true);
	if (input.startsWith("220")==false)
		return sendActionError("Server not able to accept connection - Try again later");

	//Now we have to get the server's identity, which is the next word only
	try {
		int startIndex = input.indexOf(" ");
		int endIndex = input.indexOf(" ", startIndex+1);
		heloHost = input.substring(startIndex, endIndex).trim();
	}catch(Exception e) {
		return sendActionError("Unknown reply from mail server");
	}
  
   	//Now skip any other 220 messages that are sent
 	while (input.startsWith("220"))input = getServerLine(false);
	return true;
}


 



////////////////////////////////////////////////////////
//
// Send a message to the mail server
//
////////////////////////////////////////////////////////
/**
 * Send a single message to the server
 */
 private boolean sendToServer(String command, String data, int code, boolean report) {
 	String message = new String(command);
	if (data!=null)message = message + " " + data;
 	if (report)sendActionMessage(message);

    try {
    	dataOutput.println(message);
    	dataOutput.flush();
    }
	catch(Exception e) {
		return sendActionError("Could not write to mail server");
	}
	if (code!=0) return receiveFromServer("" + code, command);
	return true;
 }

////////////////////////////////////////////////////////
//
// Receive a message from the server and check it
// against our anticipated code
//
////////////////////////////////////////////////////////
/**
 * Receive and check the response from the server
 */
 private boolean receiveFromServer(String expectedCode, String error) {
	String input = getServerLine(true);
   	if (input.startsWith(expectedCode)==false)
   		return sendActionError(input);
    while (input.startsWith(expectedCode))input = getServerLine(false);
	return true;
}





////////////////////////////////////////////////////////
//
// Read a Line of Data from the server
//
////////////////////////////////////////////////////////
/**
 * Read a line of data from the server
 */
 private String getServerLine(boolean waitForReply) {
 	try {
 		if (waitForReply==false && dataInput.ready()==false)return "";
 		return  dataInput.readLine();
  	}catch(Exception e) {}
  	return "";
}



////////////////////////////////////////////////////////
//
// Open up a connection to the server
//
////////////////////////////////////////////////////////
/**
 * Open Connection to the Server
 */
 private boolean openConnection() {
    try{
    	socket = new Socket(mailServer, mailPort);
 	 	dataInput = new BufferedReader(new InputStreamReader(socket.getInputStream()));
		dataOutput = new PrintWriter(socket.getOutputStream());
		return true;
    } 
    catch (Exception e){
		return sendActionError("Can't connect to mail server " + e.toString());
	}
}

////////////////////////////////////////////////////////
//
// Close connection to server
//
////////////////////////////////////////////////////////
/**
 * Close Connection to the Server
 */
 private void closeConnection() {
 	try {
 		dataOutput.close();
		dataInput.close();
		socket.close();
 	}catch(Exception e) {}
}
        




////////////////////////////////////////////////////////
//
// Is our data valid
//
////////////////////////////////////////////////////////
/**
 * Is our data valid
 */
 private boolean isDataValid() {
	String result = null;
	if (mailMessage==null)result= "No Mail Message";
	if (mailSubject==null)result= "No Mail Subject Line";
	if (mailTo.length<=0)result= "No Mail 'To' Address";
	if (mailFrom ==null)result= "No Mail 'From' Address";
 	if (mailServer==null)result= "No Server Name";

	if (result!=null)return sendActionError(result);
	return true;
}


////////////////////////////////////////////////////////
//
// Send an Error Message to our listener and close the
// connection
//
////////////////////////////////////////////////////////
/**
 * Send a terminal error message
 */
 private boolean sendActionError(String message) {
  	//Make sure our connections are closed before returning
	closeConnection();
	ActionEvent event = new ActionEvent(this, ERROR, message);
	if (actionListener!=null)actionListener.actionPerformed(event);
	return false;
}

////////////////////////////////////////////////////////
//
// Send a message to our action listener
//
////////////////////////////////////////////////////////
/**
 * Send a message to out action listener
 */
 private boolean sendActionMessage(String message) {
 	ActionEvent event = new ActionEvent(this, MESSAGE, message);
	if (actionListener!=null)actionListener.actionPerformed(event);
	return false;
}
    
////////////////////////////////////////////////////////
//
//  Validate email data
//
////////////////////////////////////////////////////////
/** 
 * Split a piece of text into an array of lines
 * suitable for transmission as DATA i.e., making sure
 * that there are no lines that contain just a full stop '.'
 */
 private String[] validateMailData(String message) {
	Vector lineList = new Vector();
	try {
 	 	StringReader reader = new StringReader(message);
		BufferedReader lineReader = new BufferedReader(reader);
		String aLine=null;
		do {
			aLine = lineReader.readLine();
			if (aLine!=null) {
				if (aLine.equals("."))aLine = ". ";
				lineList.addElement(aLine);
			}
		}while(aLine!=null);
	} catch(Exception e) {}
	String[] results = new String[lineList.size()];
	lineList.copyInto(results);
	return results;
}

}