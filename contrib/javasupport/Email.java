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

/**
 *
 * Email class
 * <p>
 * Simple email class that stores all parts of an Email message.
 * Messages can either be created and Sent using SendMail, or can
 * be received and interrogated using ReceiveMail.
 * <p>
 *
 * @version 2.0 August 1998
 * @author Chris Parkinson
 */


public class Email 
{

	private String 			mailFrom;   	//Who is sending this message
	private Vector 			mailTo;   		//Who are the recipients
	private Vector			mailCC;			//Who are we CCing to
	private Vector			mailHeader;	    //Other headers, such as X-...
			
	private String 			mailSubject; 	//The subject line of the message
	private StringBuffer 	mailMessage;	//The actual mail message

  
////////////////////////////////////////////////////////
//
// Constructor
//
////////////////////////////////////////////////////////
/**
 * Construct an empty email message
 */
 public Email() {
	mailTo = new Vector();
	mailCC = new Vector();
	mailHeader = new Vector();
	mailMessage = new StringBuffer();
 }



////////////////////////////////////////////////////////
//
// Set and Get the from email address
//
////////////////////////////////////////////////////////
/**
 * Set the mail sender.
 * <P>
 * @param from the email address of the user sending the mail
 */
 public void setMailFrom(String from) {
 	mailFrom = removeIdentifier("From:", from);
}

/**
 * Get the mail sender.
 * <P>
 * @return the email address of the user sending the mail
 */
 public String getMailFrom() {
    return mailFrom;
}
 

////////////////////////////////////////////////////////
//
// Set and Get the to email recipient vector
//
////////////////////////////////////////////////////////
/**
 * Get the 'to' email address(s) as an array of String addresses
 * <P>
 * @return a String array containing the recipients.
 */
 public String[] getMailTo() {
	String result[] = new String[ mailTo.size()];
	mailTo.copyInto(result);
	return result;
}


/**
 * Set the 'to' email address(s) as an array of String addresses
 * <P>
 * @param a String array containing the recipients.
 */
 public void setMailTo(String array[]) {
	if (array==null)return;
	for (int i=0; i<array.length; i++)	setMailTo(array[i]);
}

/**
 * Add a 'to' email address as a String address
 * <P>
 * @param a String containing a recipient.
 */
 public void setMailTo(String recipient) {
	String address = removeIdentifier("To:", recipient);
	if (address==null)return;
	if (address.equals(""))return;
	for (int i=0; i<mailTo.size(); i++) {
		String existing = (String)mailTo.elementAt(i);
		if (existing.equals(address))return;
	}
	mailTo.addElement(address);
}

   

////////////////////////////////////////////////////////
//
// Set and Get the CC email recipient vector
//
////////////////////////////////////////////////////////
/**
 * Get the 'CC' email address(s) as an array of String addresses
 * <P>
 * @return a String array containing the CC recipients.
 */
 public String[] getMailCC() {
	String result[] = new String[ mailCC.size()];
	mailCC.copyInto(result);
	return result;
}


/**
 * Set the 'CC' email address(s) as an array of String addresses
 * <P>
 * @param a String array containing the CC recipients.
 */
 public void setMailCC(String array[]) {
	if (array==null)return;
	for (int i=0; i<array.length; i++)  setMailCC(array[i]);
}
   

/**
 * Add a 'CC' email address as a String address
 * <P>
 * @param a String containing a CC recipient.
 */
 public void setMailCC(String recipient) {
 	String address = removeIdentifier("cc:", recipient);
	if (address==null)return;
	if (address.equals(""))return;
	for (int i=0; i<mailCC.size(); i++) {
		String existing = (String)mailCC.elementAt(i);
		if (existing.equals(address))return;
	}
	for (int i=mailTo.size()-1; i>=0; i--) {
		String existing = (String)mailTo.elementAt(i);
		if (existing.equals(address))return;
	}

	mailCC.addElement(address);
}



////////////////////////////////////////////////////////
//
// Set and Get the additional Headers for this email
//
////////////////////////////////////////////////////////
/**
 * Get the additional headers as an array of String addresses
 * <P>
 * @return a String array containing the headers.
 */
 public String[] getMailHeader() {
	String result[] = new String[ mailHeader.size()];
	mailHeader.copyInto(result);
	return result;
}


/**
 * Set the additional mail headers as an array of String addresses
 * <P>
 * @param a String array containing the headers.
 */
 public void setMailHeader(String array[]) {
	if (array==null)return;
	for (int i=0; i<array.length; i++)	setMailHeader(array[i]);
}
   

/**
 * Add an extra mail header as a String
 * <P>
 * @param a String containing a header.
 */
 public void setMailHeader(String header) {
	if (header==null)return;
	mailHeader.addElement(header);
}





////////////////////////////////////////////////////////
//
// Set and Get the email subject line
//
////////////////////////////////////////////////////////
/**
 * Set the email subject line.
 * <P>
 * @param subject the subject line of the mail message
 */
 public void setMailSubject(String subject) {
	if (subject==null)subject="";
 	mailSubject = removeIdentifier("Subject:", subject);
}


/**
 * Get the email subject line.
 * <P>
 * @return the subject line of the mail message
 */
 public String getMailSubject() {
    return mailSubject;
}



////////////////////////////////////////////////////////
//
// Set and Get the mail message
//
////////////////////////////////////////////////////////
/**
 * Set the mail message
 * <P>
 * @param message the text of the mail message
 */
 public void setMailMessage(String message) {
	if (message!=null)	mailMessage = new StringBuffer(message);
}

/**
 * Add text to the mail message
 * <P>
 * @param message a line of text to add to the existing mail message
 */
 public void addMailMessage(String message) {
	mailMessage.append(message);
	mailMessage.append("\n");
}	  


/**
 * Get the mail message
 * <P>
 * @return  the text of the mail message
 */
 public String getMailMessage() {
 	return mailMessage.toString();
}



////////////////////////////////////////////////////////
//
// Remove Identifier from Line
//
////////////////////////////////////////////////////////
/**
 * Utility method to remove a line identifier. This is used
 * mainly for decoding email messages from the POP3 server. Here,
 * lines will arrive such as "From: chris@pnl.gov". So we use this method
 * to get rid of the "From" and just return the data
 */
 public String removeIdentifier(String identifier, String rawString) {
	if (rawString==null)return null;
	if (identifier==null)return rawString.trim();

	String newRawString = new String(rawString);

 	//Make everything lowercase
	String newID = identifier.toLowerCase();
	String newRaw = newRawString.toLowerCase();
	if (newRaw.startsWith(newID)==false)return newRawString;
	String result = newRawString.substring(identifier.length());
   	result = result.trim();
	if (result.equals(""))return null;
	return result;
}


////////////////////////////////////////////////////////
//
// To String
//
////////////////////////////////////////////////////////
/**
 * toString
 * <P>
 * @return a string representation of the email
 */
 public String toString() {
 	return getMailSubject();
}


}
  
