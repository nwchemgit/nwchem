// $Id$

import java.net.Socket;
import java.io.*;
import java.net.*;
import java.util.*;
import java.lang.*;
import java.awt.*;
import java.awt.Color;
import java.awt.event.*;
import java.awt.image.*;

public class NWChem extends Frame
		 {
/******all window event handling**********/

public void windowClosed(WindowEvent event) {}
public void windowDeiconified(WindowEvent event) {}
public void windowIconified(WindowEvent event) {}
public void windowActivated(WindowEvent event) {}
public void windowOpened(WindowEvent event) {}
public void windowClosing(WindowEvent event) {
	setVisible(false);
	super.dispose();
	System.exit(0);
}


public NWChem() {

	super("NWChem support");
	
/***************The menubar currently has only an "Exit" button***********/
	super.setBackground(new Color(202,225,255));
	super.setForeground(new Color(139,0,0));
	MenuBar mb = new MenuBar();
	Menu fileMenu = new Menu("File");
	fileMenu.add(new MenuItem("Exit "));
	mb.add(fileMenu);
	setMenuBar(mb);

	setLayout(null);
	setSize(700,700);
	frame = new Frame("");
//	frame.setSize(600,360);	
	
/****************Labels 1-14 are  the text on the frame******************/
	label1=new Label("nwchem-support@emsl.pnl.gov  Request"); //Main title
	label1.setFont(new Font("Dialog",Font.BOLD,14));
	add(label1);
	label1.setBounds(12,55,650,30);

	label2=new Label("Complete the following request form:");
	label2.setFont(new Font("Dialog",Font.BOLD,12));
	add(label2);
	label2.setBounds(12,85,650,20);

	label3=new Label("To       :   nwchem-support@emsl.pnl.gov");
	label3.setFont(new Font("Dialog",Font.BOLD,12));
	add(label3);
	label3.setBounds(12,115,318,13);

	label4=new Label("From   :");
	label4.setFont(new Font("Dialog",Font.BOLD,12));
	add(label4);
	label4.setBounds(12,140,60,20);

	label5=new Label("(e-mail address required)");
	label5.setFont(new Font("Dialog",Font.PLAIN,12));
	add(label5);
	label5.setBounds(306,140,200,20);

	label6=new Label("Subject:");
	label6.setFont(new Font("Dialog",Font.BOLD,12));
	add(label6);
	label6.setBounds(12,170,60,20);

	label7=new Label("(required)");
	add(label7);
	label7.setBounds(390,170,60,20);

	label8=new Label("cc         :");
	label8.setFont(new Font("Dialog",Font.BOLD,12));
	add(label8);
	label8.setBounds(12,205,60,20);

	label10=new Label("Receipt:");
	label10.setFont(new Font("Dialog",Font.BOLD,12));
	add(label10);
	label10.setBounds(12,240,60,20);

	label11=new Label("Help us to process your request faster by selecting annotations:");
	label11.setFont(new Font("Dialog",Font.BOLD,12));
	add(label11);
	label11.setBounds(12,280,650,20);

	label12=new Label("Enter your support request: (Include reason and need-by-date if priority is vital or emergency)");
	label12.setFont(new Font("Dialog",Font.BOLD,12));
	add(label12);
	label12.setBounds(12,390,655,20);

	label13 = new Label("Input attachment:");
	label13.setFont(new Font("Dialog",Font.BOLD,12));
	add(label13);
	label13.setBounds(7,550,150,30);

	label14 = new Label("Output attachment:");
	label14.setFont(new Font("Dialog",Font.BOLD,12));
	add(label14);
	label14.setBounds(7,585,150,30);

	label16 = new Label("Other attachment:");
	label16.setFont(new Font("Dialog",Font.BOLD,12));
	add(label16);
	label16.setBounds(7,620,150,30);
	
	label15 = new Label("Danielle Farrar, 1998");
	label15.setFont(new Font("SansSerif",Font.BOLD,12));
	add(label15);
	label15.setBounds(550,650,120,12);

	String picture = "pc22.gif";
	ImageViewer newimage = new ImageViewer(picture);
	add(newimage);
	newimage.setBounds(580,586,56,64);

	

/*****************Each text field is a part of the message****************/
	edit1 = new TextField(37);
	add(edit1);
	edit1.setBounds(78,135,210,30);

try {
	String user = System.getProperty("user.home");
	String nwpath = user + "/nwchemsupport";
	File nwfile = new File(nwpath);

	RandomAccessFile emailfile = new RandomAccessFile(nwfile, "r");
	String sentfrom = emailfile.readLine();
	edit1.setText(sentfrom);
	} catch (IOException e) {System.out.println("Remember to set your email address to the 'nwchemsupport' file in your home directory!");}
	
	edit2=new TextField(37);
////////Subject?
	add(edit2);
	edit2.setBounds(78,165,306,30);

	edit3=new TextField(18);
////////cc?
	add(edit3);
	edit3.setBounds(78,200,150,30);

	edit4=new TextArea(7,64);
////////additional message field	
	add(edit4);
	edit4.setBounds(12,420,650,100);

	edit5=new TextField(68);
////////input file attachment
	add(edit5);
	edit5.setBounds(155,550,210,30);

	edit6= new TextField(68);
////////output file attachment
	add(edit6);
	edit6.setBounds(155,585,210,30);

	edit7 = new TextField(18);
////////ccarea
	edit7.setBounds(250,200,150,30);
	add(edit7);

	edit8 = new TextField(18);
////////ccarea
	edit8.setBounds(420,200,150,30);
	add(edit8);

	edit9 = new TextField(18);
////////another attachment directory/filename
	edit9.setBounds(155,620,210,30);
	add(edit9);

	group1= new CheckboxGroup();
////////receipt request checkboxes
	check1=new Checkbox("No",group1, true);
	add(check1);
	check1.setBounds(78,240,36,20);
	check2=new Checkbox("Yes",group1, false);
	add(check2);
	check2.setBounds(120,240,48,20);


/****************choice boxes**************************/

	choice1= new Choice();
	add(choice1);
	choice1.setBounds(12,308,132,26);
	choice1.addItem("user_priority");
	choice1.addItem("low");
	choice1.addItem("normal");
	choice1.addItem("secondary-importance    ");
	choice1.addItem("primary-importance");
	choice1.addItem("vital");
	choice1.addItem("emergency");

	choice2= new Choice();
	add(choice2);// second choice box, describes the problem
	choice2.setBounds(200,308,132,26);
	choice2.addItem("Problem?");
	choice2.addItem("distribution");
	choice2.addItem("execution");
	choice2.addItem("fatal-bug");
	choice2.addItem("file");
	choice2.addItem("help-request");
	choice2.addItem("info-request");
	choice2.addItem("input");
	choice2.addItem("install");
	choice2.addItem("non-fatal-bug");
	choice2.addItem("nwchem-configuration");
	choice2.addItem("output");
	choice2.addItem("performance");
	choice2.addItem("suggestion");
	choice2.addItem("web_page_problem");
	choice2.addItem("unknown");
	choice2.addItem("developer-note          ");

	choice3= new Choice();
	add(choice3);
	choice3.setBounds(12,344,132,26);
	choice3.addItem("Module/Theory?");
	choice3.addItem("all");
	choice3.addItem("atomic_scf");
	choice3.addItem("basis_object");
	choice3.addItem("blas_lapack");
	choice3.addItem("ccsd");
	choice3.addItem("chemio");
	choice3.addItem("columbus");
	choice3.addItem("cphf");
	choice3.addItem("dft");
	choice3.addItem("dplot");
	choice3.addItem("driver_geometry_opt");
	choice3.addItem("ecp");
	choice3.addItem("esp_charge_fit");
	choice3.addItem("fft");
	choice3.addItem("gapss");
	choice3.addItem("geometry_object");
	choice3.addItem("global_array");
	choice3.addItem("input_module");
	choice3.addItem("integrals");
	choice3.addItem("lapi");
	choice3.addItem("makefiles");
	choice3.addItem("mcscf");
	choice3.addItem("md_nwargos");
	choice3.addItem("memory_allocator");
	choice3.addItem("moints");
	choice3.addItem("mp2");
	choice3.addItem("mrpt");
	choice3.addItem("nbo");
	choice3.addItem("peigs");
	choice3.addItem("prepar");
	choice3.addItem("properties");
	choice3.addItem("pstat");
	choice3.addItem("qmmm");
	choice3.addItem("rimp2");
	choice3.addItem("riscf");
	choice3.addItem("rtdb");
	choice3.addItem("scf");
	choice3.addItem("selected_ci");
	choice3.addItem("stepper_geometry_opt");
	choice3.addItem("symmetry");
	choice3.addItem("tasks");
	choice3.addItem("tcgmsg_or_tcgmsg_mpi");
	choice3.addItem("util_routine");
	choice3.addItem("vib");
	choice3.addItem("independent");

	choice4= new Choice(); 
////////fourth choice box, describes operation
	add(choice4);
	choice4.setBounds(380,308,132,26);
	choice4.addItem("Operation?");
	choice4.addItem("dynamics");
	choice4.addItem("energy");
	choice4.addItem("frequencies");
	choice4.addItem("gradient");
	choice4.addItem("independent                ");
	choice4.addItem("optimize");
	choice4.addItem("saddle");
        choice4.addItem("shell");
	choice4.addItem("thermodynamics");
	choice4.addItem("unknown");

	choice5= new Choice();
////////fifth choice box, describes platform
	add(choice5);
	choice5.setBounds(200,344,132,26);
	choice5.addItem("Platform?");
	choice5.addItem("c90_cray");
	choice5.addItem("dec");
	choice5.addItem("hp");
	choice5.addItem("ibm");
	choice5.addItem("ibmsp");
	choice5.addItem("j90_cray");
	choice5.addItem("linux");
	choice5.addItem("paragon");
	choice5.addItem("sgi32");
	choice5.addItem("sgi64");
	choice5.addItem("sun");
	choice5.addItem("t3d_cray");
	choice5.addItem("t3e_cray");
	choice5.addItem("other");
	choice5.addItem("independent                 ");

	button1 = new Button("Submit"); 
////////the submit button, sends e-mail by clickedButton1() method
	add(button1);
	button1.setBounds(12,650,78,26);

	button2=new Button("Reset");
////////reset button, clears screen by clickedButton2() method
	add(button2);
	button2.setBounds(99,650,78,26);

	button3= new Button("Attach File");
////////attach file calls up dialog box to pick directory/file for attachment, by clickedButton3() method,
////////inserts into textfield edit5
	button3.setBounds(380,550,85,30);
	add(button3);

	button4= new Button("Attach File");
////////calls up dialog box, exactly the same as button3, except inserts directory/file into textfield edit6
	button4.setBounds(380,585,85,30);
	add(button4);

	button5= new Button("Attach File");
////////calls up dialog box, exactly the same as button3, except inserts directory/file into textfield edit6
	button5.setBounds(380,620,85,30);
	add(button5);


	setVisible(true); 	



/*******************Event handling for all awt tools**************/

/////////listens for when textfield value changed, but no action is performed///////////////
	edit1.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e)
{ }});

	edit2.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e)
{ }});
	edit3.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e)
{ }});
	edit4.addTextListener(new TextListener() {
	public void textValueChanged(TextEvent t){ 
	 } });
	
	edit5.addTextListener(new TextListener() {
	public void textValueChanged(TextEvent t){ 
	 } });

	edit6.addTextListener(new TextListener() {
	public void textValueChanged(TextEvent t){ 
	 } });


	button1.addActionListener(new ActionListener()
	{
 	public void actionPerformed(ActionEvent e) {
	if (e.getSource() == button1) {
	clickedButton1(); }}});

	button2.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e){
	if (e.getSource() == button2) {
	clickedButton2(); }}});

	button3.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e){
	if (e.getSource() == button3) {
	clickedButton3();
}}});

	button4.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e){
	if (e.getSource() == button4) {
	clickedButton4(); }}});

	button5.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e){
	if (e.getSource() == button5) {
	clickedButton5(); }}});

	choice1.addItemListener(new ItemListener() {
	public void itemStateChanged(ItemEvent e) {
}});

	choice2.addItemListener(new ItemListener() {
	public void itemStateChanged(ItemEvent e) { 

}});

	choice3.addItemListener(new ItemListener() {
	public void itemStateChanged(ItemEvent e) { 

}});

	choice4.addItemListener(new ItemListener() {
	public void itemStateChanged(ItemEvent e) { 

}});

	
	choice5.addItemListener(new ItemListener() {
	public void itemStateChanged(ItemEvent e) { 

}});	


/////listens for when "File; Exit" selected, then closes frame//////////
	fileMenu.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) 
		{
		if ( e.getSource() instanceof MenuItem){
			if (e.equals("Exit")); {
				Exit();}}}});}


    public static void main(String args[]) {
	NWChem NWC = new NWChem();
	}


 
/*********declaring all variables**************/
    	Label label1;
    	Label label2;
    	Label label3;
    	Label label4;
    	Label label5;
    	Label label6;
    	Label label7;
    	Label label8;
    	Label label9;
    	Label label10;
    	Label label11;
    	Label label12;
    	Label label13;
    	Label label14;
	Label label15;
	Label label16;
	Label variable;
    	TextField edit1;
    	TextField edit2;
    	TextField edit3;
    	TextArea edit4;
    	TextField edit5;
    	TextField edit6;
    	TextField edit7;
    	TextField edit8;
	TextField edit9;
    	CheckboxGroup group1;
    	Checkbox check1;
    	Checkbox check2;
    	Choice choice1;
    	Choice choice2;
    	Choice choice3;
    	Choice choice4;
    	Choice choice5;
    	Button button1;
    	Button button2;
    	Button button3;
    	Button button4;
	Button button5;
    	FileDialog openDialog;
    	Frame frame;
	String sent;
	Image image;
	int width, height;

    	public synchronized void setVisible() {
		setLocation(50, 50);
		}// sets the visibility of the window.  change location??

public void Exit() { 
////////closes window when "Exit" has been selected///////////
	super.dispose();
	System.exit(0);
	setVisible(false);
}



public void clickedButton3() {
////////calls up dialog box, when directory/file selected, they are inserted into edit5 textField
	FileDialog openDialog;
	openDialog = new FileDialog(frame, "Text:Open", FileDialog.LOAD);
	openDialog.setVisible(true);
	String FileName = openDialog.getDirectory() + openDialog.getFile();
	File file = new File(FileName);
	if (file.exists()) {
	edit5.setText(FileName);}
	else {edit5.setText("");}

	
}

public void clickedButton4() {
////////also calls up dialog box, edit6 textfield set to selected directory/file
	FileDialog openDialog;
	openDialog = new FileDialog(frame, "Text:Open", FileDialog.LOAD);
	openDialog.setVisible(true);
	String FileName = openDialog.getDirectory() + openDialog.getFile();
	File file = new File(FileName);
	if (file.exists()) {
	edit6.setText(FileName);}
	else { edit6.setText("");}
}


public void clickedButton5() {
////////also calls up dialog box, edit6 textfield set to selected directory/file
	FileDialog openDialog;
	openDialog = new FileDialog(frame, "Text:Open", FileDialog.LOAD);
	openDialog.setVisible(true);
	String FileName = openDialog.getDirectory() + openDialog.getFile();
	File file = new File(FileName);
	if (file.exists()) {
	edit9.setText(FileName);}
	else { edit9.setText("");}
}




public void clickedButton1() 
    { 
///////sends an e-mail to nwchem with information from textfields, choiceboxes, and attachments
///////must be in a try-catch loop to send file information
    try {
	waiting waiting2 = new waiting();
	waiting2.setVisible(true);

////////class "Email" sends to designated address

	Email support = new Email();

////////e-mail address from edit1 is set as sender
	String sending = edit1.getText().toString();
	support.setMailFrom(sending);
	
	notsent notsentwin = new notsent();



/////////string "to" can be set to any e-mail address	
	String to = "nwchem-support@emsl.pnl.gov";
	support.setMailTo(to);

////////sends the courtesy copies if textfields have been filled
	if (edit3.getText().length() > 0)
	{	support.setMailCC(edit3.getText()); 	}
	if (edit7.getText().length() > 0)
	{	support.setMailCC(edit7.getText());	}
	if (edit8.getText().length() > 0)
	{ 	support.setMailCC(edit8.getText());	}

////////subject is set to text from textfield edit2
	String subject = edit2.getText().toString();
	support.setMailSubject(subject);

	if (edit1.getText().toString() == "") 
		{
		waiting2.setVisible(false);
		notsentwin.setVisible(true);
		} 

	if (edit2.getText().length() <= 0)
		{
		waiting2.setVisible(false);
		notsentwin.setVisible(true);
		} 


	String start = "";
	support.setMailMessage(start);

////////if "Yes" has been checked off under "Receipt:", "X-Special: Receipt" flag will be added to message
	if (group1.getSelectedCheckbox() == check2 )
	{
	String receipt = "X-Special: Receipt";
	support.setMailHeader(receipt);
	}
////////all other flags
	String bug = "X-Bug: " + choice2.getSelectedItem();
	support.setMailHeader(bug);

	String Module = "X-Module: " + choice3.getSelectedItem();
	support.setMailHeader(Module);

	String Operation = "X-Operation: " + choice4.getSelectedItem();
	support.setMailHeader(Operation);

	String Platform = "X-Platform: " + choice5.getSelectedItem();
	support.setMailHeader(Platform);

	String Priority = "X-Priority: " + choice1.getSelectedItem();
	support.setMailHeader(Priority);
	
    	String message = "--- QMH Annotations ---";
//    	support.addMailMessage(message);

////////choice box information will be sent only if choice has been changed from default.	
/*    	if (choice1.getSelectedItem() != "user_priority")	{
        	String message0 = "User Priority\t:\t" + choice1.getSelectedItem().toString();
       		support.addMailMessage(message0); 
        	}
    	if (choice2.getSelectedItem() != "Problem?")	{
        	String message1 = "Problem\t\t:\t" + choice2.getSelectedItem().toString();
        	support.addMailMessage(message1);
    		}
    	if (choice3.getSelectedItem() != "Module/Theory?")	{
        	String message2 = "Module/Theory\t:\t" + choice3.getSelectedItem().toString();
        	support.addMailMessage(message2);
    		}
    	if (choice4.getSelectedItem() != "Operation?")	{
        	String message3 = "Operation\t:\t" + choice4.getSelectedItem().toString();
        	support.addMailMessage(message3);
    		}
    	if (choice5.getSelectedItem() != "Platform?")	{
        	String message4 = "Platform\t:\t" + choice5.getSelectedItem().toString();
        	support.addMailMessage(message4);
    		}
*/
////////message body, including date, sender, reciever, cc, subject, etc
	String message5 ="\n" + "--- Request Message ---";
//	support.addMailMessage(message5);
    	Date thedate = new Date();
   	String message6 = "Date\t\t:\t" + thedate.toString();
//   	support.addMailMessage(message6);
    	String message7 = "To\t\t\t:\t nwchem-support@emsl.pnl.gov";
//    	support.addMailMessage(message7);
 //   	String message8 = "From\t\t:\t" + sender;
//   	support.addMailMessage(message8);
        String message9 = "cc\t\t\t:\t" + edit3.getText().toString() + edit7.getText() + edit8.getText();
//        support.addMailMessage(message9);
    	String message10 = "Subject\t\t:\t" + subject;
//    	support.addMailMessage(message10);
    	String message11 = "\n" + edit4.getText().toString() + "\r\n\r\n";
    	support.addMailMessage(message11); 

////////retrieves filename from textfield edit5 to send "input file" as attachment appended to e-mail,
////////attachment doesn't need to be included for message to be sent
	File filename = new File(edit5.getText());
	if (filename.exists()) {
		String message12 = "<>*<>*<>*<>*<>*<>*<>*<>*<>*<>*<>*<>*INPUT FILE<>*<>*<>*<>*<>*<>*<>*<>*<>*<>*<>*<>*\r\n\r\n\r\n\r\n";
		support.addMailMessage(message12);
		RandomAccessFile File = new RandomAccessFile(filename, "rw");
		String attach0 = edit5.getText().toString();
		while (File.getFilePointer() < File.length())	{
			String message13 =  File.readLine();
			support.addMailMessage(message13);  
			}
		}

////////retrieves filename from textfield edit6 to send "output file" as an additional attachment appended to e-mail
	File filename1 = new File(edit6.getText());
	if (filename1.exists()) {
		String messagenull = "\r\n\r\n";
		support.addMailMessage(messagenull);
		String message14 = "<>*<>*<>*<>*<>*<>*<>*<>*<>*<>*<>*<>*OUTPUT FILE<>*<>*<>*<>*<>*<>*<>*<>*<>*<>*<>*<>*\r\n\r\n\r\n\r\n";
		support.addMailMessage(message14);
		RandomAccessFile File1 = new RandomAccessFile(filename1, "rw");
		String attach = edit6.getText();
		while (File1.getFilePointer() < File1.length())	{
			String message15 = File1.readLine();
			support.addMailMessage(message15);
	}}

	File filename2 = new File(edit9.getText());
	if (filename2.exists()) {
		String message16 = "\r\n\r\n";
		support.addMailMessage(message16);
		String message17 = "<>*<>*<>*<>*<>*<>*<>*<>*<>*<>*<>*<>*OTHER FILE<>*<>*<>*<>*<>*<>*<>*<>*<>*<>*<>*<>*\r\n\r\n\r\n\r\n";
		support.addMailMessage(message17);
		RandomAccessFile File2 = new RandomAccessFile(filename2, "rw");
		String attach1 = edit9.getText();
		while (File2.getFilePointer() < File2.length()) {
			String message18 = File2.readLine();
			support.addMailMessage(message18);}}

/*
////////if all has been entered correctly, dialog box reading "You're message has been submitted" appears
////////and screen is reset

	if  (support.sendMail() == true) {
		waiting2.setVisible(false);
        submitted thesubmitted = new submitted();
		thesubmitted.setVisible(true);
        clickedButton2();
    }else if (support.sendMail() == false) {
		waiting2.setVisible(false);
		notsent notsentyet = new notsent();
//		notsentyet.setVisible(true);
}	*/
	
/**
 * Modification by Chris Parkinson  26th March 1999
 * Bug fixing old Email classes. The Email class has now
 * been split into 2 : a Email class that just stores the
 * contents of an Email; and a SendMail class that sends
 * a single Email object. 
 * To accomodare these changes the following code change has been made
 */

//Send the Email
	SendMail sendMail = new SendMail(null);
	sendMail.setEmail(support);
	sendMail.start(); 	    //Start the mailing process
	waiting2.setVisible(false);

	//Did this get sent successfully
	if (sendMail.isCancelled()==false) {
	  	submitted thesubmitted = new submitted();
		thesubmitted.setVisible(true);
        clickedButton2();
	}

} catch (IOException e) {System.out.println("Error:" + e.toString());}
		
} 

public void clickedButton2() {
////////clears all textfields and choiceboxes back to default

try {
	String user = System.getProperty("user.home");
	String nwpath = user + "/nwchemsupport";
	File nwfile = new File(nwpath);

	RandomAccessFile emailfile = new RandomAccessFile(nwfile, "rw");
	String sentfrom = emailfile.readLine();
	edit1.setText(sentfrom);
	} catch (IOException e) {}

    	edit2.setText("");
    	edit3.setText("");
    	edit4.setText("");
    	edit5.setText("");
    	edit6.setText("");
	edit7.setText("");
	edit8.setText("");
	edit9.setText("");
    	group1.setSelectedCheckbox(check1);
    	choice1.select("user_priority");
    	choice2.select("Problem?");
    	choice3.select("Module/Theory?");
    	choice4.select("Operation?");
    	choice5.select("Platform?");


   
 }

}     
   
