// $Id$

import java.awt.*;
import java.awt.event.*;

public class submitted extends Frame {

    public submitted() {

        super("   ");

        //{{INIT_CONTROLS
        setLayout(null);
        addNotify();
        setSize(400,200);
        label1=new Label("Your request has been submitted");
        label1.setFont(new Font("Dialog",Font.BOLD,14));
        add(label1);
        label1.setBounds(65,75,270,25);
        button1=new Button("Okay");
        add(button1);
        button1.setBounds(150,100,70,30);
        //}}
	button1.addActionListener(new ActionListener()
	{
	public void actionPerformed(ActionEvent e) {
	clickedButton1();
	}});
 

    }

    public synchronized void setVisible() {
    	setLocation(50, 50);
  
    }

    Label label1;
    Button button1;

    public void clickedButton1() {
	
	     setVisible(false);
    }
	public void windowClosing(WindowEvent event){ }
	public void windowClosed(WindowEvent event) { }
	public void windowDeiconified(WindowEvent event) {}
	public void windowIconified(WindowEvent event) {}
	public void windowActivated(WindowEvent event) {}
	public void windowOpened(WindowEvent event) {}
}


