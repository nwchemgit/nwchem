// $Id$

import java.awt.*;
import java.awt.event.*;

public class notsent extends Frame {
	Label label1;
	Label label2;
	Button button1;

	public notsent() {
	super("    ");

	setLayout(null);
	addNotify();
	setSize(500,200);
	label1 = new Label("You have not completed all the necessary information");
	label1.setFont(new Font("Dialog",Font.BOLD,14));
	add(label1);
	label1.setBounds(45,80,500,20);
	label2 = new Label("Please go back and fill out all required fields.");
	label2.setFont(new Font("Dialog",Font.BOLD,14));
	add(label2);
	label2.setBounds(70,100,500,20);
	button1= new Button("Okay");
	add(button1);
	button1.setBounds(200,130,80,30);
	button1.addActionListener(new ActionListener() { 	
		public void actionPerformed(ActionEvent e) { 
			clickedButton1(); }});}

	public synchronized void setVisible() {
		setLocation(100,100);
}

	public void clickedButton1() {
		setVisible(false);
		}
	public void windowClosing(WindowEvent event){}
	public void windowClosed(WindowEvent event) {}
	public void windowDeiconified(WindowEvent event) {}
	public void windowIconified(WindowEvent event){}
	public void windowActivated(WindowEvent event) {}
	public void windowOpened(WindowEvent event) {}
}
