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
import java.awt.event.*;
import java.util.*;



/**
 *
 * CommunicateGUI provides the underlying GUI framework for both
 * sending emails, and retreiving them.
 * <P>
 * Note that the Send or Receiving class does not need to call 'displayWindow'
 * if it doesn't want a GUI. In this case, all messages will be written to 
 * stdout
 *
 *
 * @version 1.0 August 1998
 * @author Chris Parkinson
 */


public class CommunicateGUI  implements Runnable, ActionListener
{		
	private ActionListener	actionListener;		//Who is listening
	private boolean			isGUI;

	private String 			errorMessage;		//Error message we received
	private boolean			isCancelled;		//Should we cancel operations
	private boolean			isClosed;			//Have we closed the window yet

 	private Thread 			motor;
	    
////////////////////////////////////////////////////////
//
// Constructor
//
////////////////////////////////////////////////////////
/**
 * Constructor will initialize a new CommunciateGUI Window
 * <p>
 * @param parentFrame the parent frame that owns this object
 */	
public CommunicateGUI(Frame parentFrame)  {

 	setCancelled(false);
	setError(null);
	isClosed = false;
}	

////////////////////////////////////////////////////////
//
// Are we running as a GUI
//
////////////////////////////////////////////////////////
/**
 * Are we running a GUI
 */
 public boolean isGUI() {
 	return false;
}
						   	  
								 	   	
////////////////////////////////////////////////////////
//
// Get the ProgressBar object
//
////////////////////////////////////////////////////////
/**
 * Get the progress bar object
 * <P>
 * @return a JProgressBar object
 */		
 public Component getProgressBar() {
	return null;
}
	 
								 	   	
////////////////////////////////////////////////////////
//
// Get the Text Message object
//
////////////////////////////////////////////////////////
/**
 * Get the textArea for text messages
 * <P>
 * @return a GUILabel object
 */		
 public Component getTextArea() {
 	return null;
}

 						   

////////////////////////////////////////////////////////
//
// Get and Set the 'Cancel' status
//
////////////////////////////////////////////////////////
/**
 * Get the cancel status. If this returns true, then the email
 * process has been cancelled (either by user or by error)
 * and the current process should stop and exit. 
 * <P>
 * @return a boolean, true if process has been cancelled or should stop
 */		
 public boolean isCancelled() {
 	return isCancelled;
}

/**
 * Set the cancel status. If this is set to be true then the email
 * processes should stop at their next convenience and exit.
 * <P>
 * @param isCancelled a boolean, true if process should be stopped
 */		
 public void setCancelled(boolean isCancelled) {
 	this.isCancelled =  isCancelled;
	if (isCancelled==true) {
		motor = null;
		closeWindow();
	}
}



////////////////////////////////////////////////////////
//
// Set and Get an error message to display on exit
//
////////////////////////////////////////////////////////
/**
 * Set the error message to display on exit and also set our cancelled status
 * <P>
 * @param error the message to display on exit	
 * @returns a boolean result which is always false
 */		
 public boolean setError(String error) {
	if (error!=null && errorMessage!=null)return false;  //Once an error is set, keep it
 	errorMessage = error;
	setCancelled(errorMessage!=null);
	return false;
}

/**
 * Get the error message to display on exit. This will return 'null'
 * if no error message is present
 * <P>
 * @return the message to display on exit
 */		
 public String getError() {
 	return errorMessage;
}

								 	   	
////////////////////////////////////////////////////////
//
// Increment the Progress Bar
//
////////////////////////////////////////////////////////
/**
 * Increment the ProgressBar by one unit and display the given message
 * If the GUI is not active, the message will be send to stdout
 * <P>
 * @param message the text message to display in the progress bar
 */		
 public void setProgress(String message) {
  //java.lang.System.out.println(message);
}




////////////////////////////////////////////////////////
//
// Thread runner
//
////////////////////////////////////////////////////////
public void start() {
	run();
}


public void run() {
  	//System.out.println("Starting Run...");
	mainMethod();
	closeWindow();
	fireActionEvent("Email Communication Successful");
  //	System.out.println("*** Thread finished");
 }


public void mainMethod() {
}


						   	  
								 	   	
	  
////////////////////////////////////////////////////////
//
// Listen to exit buttons being pressed
//
////////////////////////////////////////////////////////
 /**
  * ActionPerformed
  */
public void actionPerformed(ActionEvent event) {
	String command = event.getActionCommand();
	if (command.equals("Cancel"))setError("Operation Cancelled by User");
	if (command.equals("Close"))setError("Operation Cancelled by User");
}
////////////////////////////////////////////////////////
//
// Display Window
//
////////////////////////////////////////////////////////
/**
 * Display the Window if we have one
 */
 public void displayWindow(Component owner) {
} 

	
////////////////////////////////////////////////////////
//
// Close the Window	and Exit
//
////////////////////////////////////////////////////////
/**
 * Close the Window and Exit. If an error occurred up to this point
 * we will display an error message
 */
 public boolean closeWindow() {
	if (isClosed==false) {
		isClosed=true;
	    //Do we have an error to display
		if (getError()!=null)displayErrorMessage(getError());
	}
	return true;
}

	
////////////////////////////////////////////////////////
//
// Display a message
//
////////////////////////////////////////////////////////
/**
 * Display an error message. If the GUI is not active, this will
 * write the message to stdout
 * <P>
 * @param errorMessage the error Message to display
 */
 public void displayErrorMessage(String errorMessage) {
	System.err.println("\nError...");
	System.err.println("    An error occurred while talking to the Mail Server :");
	System.err.println("    " + errorMessage);
}
	
 
  

/////////////////////////////////////////////////////////////////
//
// Event Processing
//
/////////////////////////////////////////////////////////////////
/** 
* Add Action Listener
*/
public void addActionListener(ActionListener l) {
	actionListener = AWTEventMulticaster.add(actionListener, l);
}
 
/** 
 * Fire an action event
 */
public void fireActionEvent(String message) {
	if (actionListener==null)return;
	ActionEvent event = new ActionEvent(this, 0, message);
	actionListener.actionPerformed(event);
}

 					
}
			











							  