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


import java.util.*;
import java.text.*;
import java.awt.*;
/**
 *
 * SendMail class ... extends from the CommunicateMail class to provide
 * method for sending an Email message via a SMTP mail server.
 * <p>
 * SMTP mail sender uses port 25 (the standard SMTP socket)
 * and adheres to  TCP/SMTP specs (RFC 821)
 *
 * @version 2.0 August 1998
 * @author Chris Parkinson
 */

public class SendMail extends CommunicateMail 
 {
	
	private Email email;				//The email message to send
 	private String heloHost;            //The server identity
  
////////////////////////////////////////////////////////
//
// Constructor
//
////////////////////////////////////////////////////////
/**
 * Construct a new SendMail object, by default setting the SMTP port to 25
 * and the mail server to 'pnl.gov'
 */
 public SendMail(Frame parentFrame) {
	super(parentFrame);
 	setMailPort(25);   //Set up some defaults
	setMailServer("pnl.gov");  
}

////////////////////////////////////////////////////////
//
// Set and Get the email object to send
//
////////////////////////////////////////////////////////
/**
 * Set the email object to send
 * <P>
 * @param email the email object to send
 */
 public void setEmail(Email email) {
 	this.email = email;
}

/**
 * Get the email object to send
 * <P>
 * @return the email object to send
 */
 public Email getEmail() {
 	return email;
}



  
////////////////////////////////////////////////////////
//
// Send the Message
//
////////////////////////////////////////////////////////
/**
 * Send the message. 
 */
public void mainMethod() {	

	//Check the validity of our data
	isDataValid();
	if (isCancelled())return;
	  
	//Open up the socket connection
	setProgress("Connecting to server (" + getMailServer()+ ")...");
	openConnection();
	if (isCancelled())return;

	//Start talking 'server talk'
	setProgress("Waiting for response from server (" + getMailServer() + ")...");

	greetServer();  		//Greet the server
 	if (isCancelled())return;

	sendHeader();	   		//Send the mail headers
	if (isCancelled())return;
	
	sendData();	 	   		//Send the data and quit commands
 	if (isCancelled())return;
	
	sendQuit();		   		//Quit communication
	if (isCancelled())return;

	//Close the connection
	closeConnection();							  
	setProgress("Connection closed");
	setProgress("Mail sent !");
}

 

						





////////////////////////////////////////////////////////
//
// Is our data valid
//
////////////////////////////////////////////////////////
/**
 * Is our data valid
 */
 private void isDataValid() {

  	if (getEmail()==null) {
  		setError("No Email Message to Send");
		return;
	}
 	String result = null;
	if (getEmail().getMailMessage()==null)result= "No Mail Message";
	if (getEmail().getMailMessage().equals(""))result= "No Mail Message";
	if (getEmail().getMailSubject()==null)result= "No Mail Subject Line";
	if (getEmail().getMailSubject().equals(""))result= "No Mail Subject Line";
	if (getEmail().getMailTo()==null)result ="No Mail 'To' Address";
	else if (getEmail().getMailTo().length<=0)result= "No Mail 'To' Address";
	if (getEmail().getMailFrom() == null)result= "No Mail 'From' Address";
	if (getEmail().getMailFrom().equals(""))result= "No Mail 'From' Address";
 	if (getMailServer()==null)result= "No Server Name";

	if (result!=null)setError(result);
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
 private void greetServer() {
   	try {
	 	String words[] = toMultiWord(getServerLine(true));
		if (words[0].equals("220")==false)  {
	  		setError("Server not able to accept connection - Try again later");
			return;
		}

		//Now we get the server's identity, which is the next word
	  	heloHost = words[1];
	 }catch(Exception e) {
	 	setError("Unknown reply from mail server");
	 	return;
	 }
  
   	//Now skip any other 220 messages that are sent	
   	//(this is just extra info from server that we cna ignore)
	String input="";
 	while (input!=null)input = getServerLine(false);
}


////////////////////////////////////////////////////////
//
// Send the Headers to the message
//
////////////////////////////////////////////////////////
/**
 * Send the headers
 */
private void sendHeader() {
	//Say Helo back to server
	sendToServer("HELO", heloHost, true);
	if (isCancelled())return;
	receiveFromServer("250"); 
	if (isCancelled())return;

   
	//Mail From
    sendToServer("MAIL FROM:",  getEmail().getMailFrom(), true);
	if (isCancelled())return;
	receiveFromServer("250"); 
	if (isCancelled())return;
 
	//Mail To (Can have multiple recipients here)
	String to[] = getEmail().getMailTo();
	for (int i =0; i<to.length; i++)   {
 	   	sendToServer("RCPT TO:",  to[i], true);
		if (isCancelled())return;

		//Reply to this can be '250' or '251', so we shall check to '25x'
 		receiveFromServer("25"); 
		if (isCancelled())return;
	}
	//Mail CC (Can have multiple recipients here)
	String cc[] = getEmail().getMailCC();
	for (int i =0; i<cc.length; i++)   {
 	   	sendToServer("RCPT TO:",  cc[i], true);
		if (isCancelled())return;

		//Reply to this can be '250' or '251', so we shall check to '25x'
 		receiveFromServer("25"); 
		if (isCancelled())return;
	}

}

					 
////////////////////////////////////////////////////////
//
// Send the actual message data
//
////////////////////////////////////////////////////////
/**
 * Send the message data
 */
private void sendData() {
	//Tell server we are sending data over
	sendToServer("DATA",  null, true);
	if (isCancelled())return;
	receiveFromServer("354"); 
	if (isCancelled())return;


	//Protocol says that first should be the subject line
	String[] lines = toMultiLine(getEmail().getMailSubject());
	setProgress("Sending data to server.");
	sendToServer("Subject:",  lines[0], false);
	if (isCancelled())return;

 
 	//Date line
	Date now = new Date();
	SimpleDateFormat formatter = new SimpleDateFormat("EEE, d MMM yyy hh:mm:ss z");
	String date = formatter.format(now);
	setProgress("Sending data to server..");
 	sendToServer("Date:",  date, false);
	if (isCancelled())return;	  

	//To list if they exist
	lines = getEmail().getMailTo();
	setProgress("Sending data to server.");
	String tos="";
	for (int i=0; i<lines.length; i++) {
		if (tos.equals("")==false) tos= tos + ", ";
		tos = tos + lines[i];
	}
	if (tos.equals("")==false) {
		sendToServer( "To:", tos,  false);
		if (isCancelled())return;
	}



  	//CC list if they exist
	lines = getEmail().getMailCC();
	setProgress("Sending data to server.");
	String ccs="";
	for (int i=0; i<lines.length; i++) {
		if (ccs.equals("")==false) ccs= ccs + ", ";
		ccs = ccs + lines[i];
	}
	if (ccs.equals("")==false) {
		sendToServer( "cc:", ccs,  false);
		if (isCancelled())return;
	}


    //Alternative headers if they exist
	lines = getEmail().getMailHeader();
	setProgress("Sending data to server.");
	for (int i=0; i<lines.length; i++) {
		sendToServer( lines[i], "", false);
		if (isCancelled())return;
	}

    




	//Now a blank line to separate the main message
 	setProgress("Sending data to server...");
	sendToServer("", null, false);      //Blank line 
	if (isCancelled())return;

	
	//Now the main body of the message (Note any single '.' should be changed to '..')
   	lines = toMultiLine(getEmail().getMailMessage());
	StringBuffer message = new StringBuffer("Sending data to server");
	for (int i =0; i< lines.length; i++)  {
		message.append(".");
		setProgress(message.toString());
 		sendToServer(lines[i], null, false);
		if (isCancelled())return;
  	}	 

	//End DATA section with '\r\n.\r\n'
	sendToServer("\r\n.", null, false);
	if (isCancelled())return;

    //receiveFromServer("250");	 //Seems to be unneccessary and locks up sometimes
  	//if (isCancelled())return;

	setProgress("Sent data");
}
	   
////////////////////////////////////////////////////////
//
// Send a Quit message
//
////////////////////////////////////////////////////////
/**
 * Send a Quit Message
 */
private void sendQuit() {
   	//Finally send a QUIT message to close the conversation
   	sendToServer("QUIT", null, true);
   	if (isCancelled())return;

	//Commenting this out  'cos its hangs and seems to be unnecessary
   //	receiveFromServer("221");
   //	if (isCancelled())return;
}
	   



}