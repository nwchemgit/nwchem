/**
 * This Java Source Code and its associated Class files are 
 * <P> 
 * (c) 1998  Battelle Pacific Northwest National Laboratories
 * <P>
 * For further information, contact :
 * <P>
 * Chris Parkinson
 * Environmental Molecular Sciences Laboratory
 * Battelle Pacific Northwest National Laboratories
 * Richland 
 * WA 99352
 * USA
 */

 
import java.awt.*;
import java.io.*;
import java.net.*;
import java.util.*;
import java.text.*;

/**
 *
 * Communicate Mail class is the abstract class that provide basic 
 * SMPT (for outgoing) and POP3 (for incoming) communication methods.
 *
 * This class is subclasses by the SendMail and ReceiveMail classes, which make
 * use of these general methods.
 * <p>
 *
 * @version 2.0 August 1998
 * @author Chris Parkinson
 */
public abstract class CommunicateMail extends CommunicateGUI 
{
	
	private int mailPort;					//SMPT PORT - should usually be 25
	private String mailServer;  			//Host name or IP address

	private Socket socket;                  //Socket to talk to host
	private BufferedReader dataInput; 		//Receive from sockets
	private PrintWriter dataOutput; 		//Send to socket

   
////////////////////////////////////////////////////////
//
// Constructor
//
////////////////////////////////////////////////////////
/**
 * Construct an empty communicator
 */
 public CommunicateMail(Frame parentFrame) {
	super(parentFrame);
}

////////////////////////////////////////////////////////
//
// Set and Get the mail port number
//
////////////////////////////////////////////////////////
/**
 * Set the mail port number. The default and usual port to use for
 * SMTP connections is 25, and for POP3 is 110
 * <P>
 * @param portNumber the port to use for SMTP or POP3 communcations.
 */
 public void setMailPort(int portNumber) {
 	mailPort = portNumber;
}

/**
 * Get the mail port number. The default and usual port to use for
 * SMTP connections is 25, and for POP3 is 110
 * <P>
 * @return the port to use for SMTP or POP3 communcations.
 */
 public int getMailPort() {
 	return mailPort;
}


////////////////////////////////////////////////////////
//
// Set the mail server (a string of name or number)
//
////////////////////////////////////////////////////////
/**
 * Set the mail server for either POP3 or SMTP
 * <P>
 * @param the name (or IP address) of the outgoing SMTP mail server or 
 * incoming POP3 server.
 */
 public void setMailServer(String server) {
 	mailServer = new String(server);
}

/**
 * Get the mail server for either POP3 or SMTP
 * <P>
 * @return the name (or IP address) of the outgoing SMTP mail server or 
 * incoming POP3 server.
 */
 public String getMailServer() {
 	return mailServer;
}


					 
					 

	   

 



////////////////////////////////////////////////////////
//
// Send a message to the mail server
//
////////////////////////////////////////////////////////
/**
 * Send a single message to the server. The message automatically has
 * a CRLF pair added to the end of the string.
 * <P>
 * @param command the command to send
 * @param data optional data to send with the command (message = command + ' ' + data)
 */
 public void sendToServer( String command,  String data,  boolean report) {
 	String message = new String(command);
	if (data!=null)message = message + " " + data;
 	if (report)setProgress(message);
    try {
 		dataOutput.print(message + "\r\n");
    	dataOutput.flush();	  
    }
	catch(Exception e) { setError("Could not write to mail server");  }
 }


////////////////////////////////////////////////////////
//
// Receive a message from the server and check it
// against our anticipated code
//
////////////////////////////////////////////////////////
/**
 * Receive a message from the server, and check the response from the server
 * using the given 'expectedResponse'. If the response is as expected, this will
 * return the whole input line. If not it returns NULL and throws an error
 * <P>
 * @param expectedResponse the expected reponse string
 * @return the message received, or NULL if error
 */
 public String receiveFromServer(String expectedResponse) {
 	String input = getServerLine(true);
	if (input==null)input="";
   	if (input.startsWith(expectedResponse))return input;
   	setError(input);
   	return null;
}


////////////////////////////////////////////////////////
//
// Read a Line of Data from the server
//
////////////////////////////////////////////////////////
/**
 * Read a line of data from the server, optionally waiting for a line
 * or returning if one is not ready.
 * <P>
 * @param waitForReply ... if true this method will block until a line of data
 * is ready to read from the server (or a time out occurs). If this is false and there
 * is no data ready to be read, the method will return a NULL.
 */
 public String getServerLine(boolean waitForReply) {
 	try {
 		if (waitForReply==false && dataInput.ready()==false)return null;
		String result = dataInput.readLine();
		return result;
  	}catch(Exception e) {
	  	if (waitForReply)setError("Could not read from mail server");
	}
  	return null;
}





////////////////////////////////////////////////////////
//
// Open up a connection to the server
//
////////////////////////////////////////////////////////
/**
 * Open Connection to the Server. Tries to open a connection to the
 * port and host as already set. 
 * <P>
 */
 public void openConnection() {
    try{
 		socket = new Socket(mailServer, mailPort);
		socket.setSoTimeout(1000*10);
 	 	dataInput = new BufferedReader(new InputStreamReader(socket.getInputStream()));
		dataOutput = new PrintWriter(socket.getOutputStream());
    } 
    catch (Exception e){
		setError("Can't connect to mail server " + e.toString());
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
 public void closeConnection() {
 	try {
 		dataOutput.close();
		dataInput.close();
		socket.close();
 	}catch(Exception e) {}
}
        


					




////////////////////////////////////////////////////////
//
//  Convert large string into a Multi Line string array
//
////////////////////////////////////////////////////////
/** 
 * Convert a large message into a multi-line message suitable
 * for use in the email system.
 * Whenever a CR/LF is found in the message, it is split into a new
 * line. If a line exists with just a single '.' in it, it is converted
 * to a '. ' so that it is not recognised as an END OF MESSAGE signal.
 * <P>
 * @param message the multiline message to split
 * @return an array of String
 */
 public String[] toMultiLine(String message) {
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

////////////////////////////////////////////////////////
//
// Tokenize a String into component words
//
////////////////////////////////////////////////////////
/** 
 * Convert a line of text into words, split by white spaces
 * <P>
 * @param message the multiword  message to split
 * @return an array of String words
 */
 public String[] toMultiWord(String message) {
	if (message==null)message="";

	Vector wordList = new Vector();
 	StringTokenizer st = new StringTokenizer(message);
    while (st.hasMoreTokens()) {
		wordList.addElement(st.nextToken());
	}
	String[] results = new String[wordList.size()];
	wordList.copyInto(results);
	return results;
}










}
    
