import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.net.*;

class nwchem_Socket extends JFrame implements ActionListener, ChangeListener, WindowListener, Runnable {
  
  Font defaultFont;
  Thread serverThread = null;
  ServerSocket server = null;

  Graph tmpPlot = new Graph();
  Graph prsPlot = new Graph();
  Graph volPlot = new Graph();
  Graph enePlot = new Graph();

  JLabel lstring = new JLabel("-----------------------------------------");

  public nwchem_Socket(){

    super("NWChem Socket");

    defaultFont = new Font("Dialog", Font.BOLD,12);

    super.getContentPane().setLayout(new GridBagLayout());
    super.getContentPane().setForeground(Color.black);
    super.getContentPane().setBackground(Color.lightGray);
    super.getContentPane().setFont(defaultFont);
    super.addWindowListener(this);

    addComponent(super.getContentPane(),tmpPlot,0,0,1,1,1,1,
    		 GridBagConstraints.BOTH,GridBagConstraints.NORTHWEST);
    addComponent(super.getContentPane(),prsPlot,0,1,1,1,1,1,
    		 GridBagConstraints.BOTH,GridBagConstraints.NORTHWEST);
    addComponent(super.getContentPane(),volPlot,1,0,1,1,1,1,
    		 GridBagConstraints.BOTH,GridBagConstraints.NORTHWEST);
    addComponent(super.getContentPane(),enePlot,1,1,1,1,1,1,
    		 GridBagConstraints.BOTH,GridBagConstraints.NORTHWEST);

    tmpPlot.init();
    tmpPlot.resize(100,100);
    tmpPlot.setNumSets(1);
    tmpPlot.setTitle("Temperature");
    tmpPlot.setXLabel("Time");

    prsPlot.init();
    prsPlot.resize(100,100);
    prsPlot.setNumSets(1);
    prsPlot.setTitle("Pressure");
    prsPlot.setXLabel("Time");

    volPlot.init();
    volPlot.resize(100,100);
    volPlot.setNumSets(1);
    volPlot.setTitle("Volume");
    volPlot.setXLabel("Time");

    enePlot.init();
    enePlot.resize(100,100);
    enePlot.setNumSets(1);
    enePlot.setTitle("Energy");
    enePlot.setXLabel("Time");

    validate();

    setLocation(25,225);	
    setSize(900,700);
    setVisible(true);

    // setup server side socket

    try{
      serverThread = new Thread(this);
      System.out.println("Thread start");
      serverThread.start();
      System.out.println("Thread done");
    } catch (Exception s) { s.printStackTrace(); };

  }

  public void start(){

  }

  public void stop(){
      System.out.println("Stop");
  }

  public void run(){
    try{
      System.out.println("Server");
      server = new ServerSocket(3333);
      System.out.println("Server thread running");
      lstring.setText("Server thread");
      validate();
      Socket client = server.accept();
      System.out.println("Client connected");
      lstring.setText("Client connected");
      validate();
      InputStreamReader isr = new InputStreamReader(client.getInputStream());
      BufferedReader is = new BufferedReader(isr);
      boolean first=true;
      double tim, tmp, prs, vol, ene;
      while(true){
	String inLine = is.readLine();
	if(inLine.length()>0){
	  tim=Double.valueOf(inLine.substring(6,17)).doubleValue();
	  ene=Double.valueOf(inLine.substring(18,29)).doubleValue();
	  tmp=Double.valueOf(inLine.substring(30,41)).doubleValue();
	  prs=Double.valueOf(inLine.substring(42,53)).doubleValue();
	  vol=Double.valueOf(inLine.substring(54,65)).doubleValue();
	  tmpPlot.addData(0,tim,tmp,!first,true);
	  prsPlot.addData(0,tim,prs,!first,true);
	  enePlot.addData(0,tim,ene,!first,true);
	  volPlot.addData(0,tim,vol,!first,true); first=false;
	};
      };
    } catch (Exception s) { s.printStackTrace(); };
  }
  
  static void addComponent(Container container, Component component,
			   int gridx, int gridy, int gridwidth, 
			   int gridheight, double weightx, 
			   double weighty, int fill, int anchor) {
	 LayoutManager lm = container.getLayout();
	 if(!(lm instanceof GridBagLayout)){
	   System.out.println("Illegal layout"); System.exit(1);
	 } else {
	   GridBagConstraints gbc = new GridBagConstraints();
	   gbc.gridx=gridx;
	   gbc.gridy=gridy;
	   gbc.gridwidth=gridwidth;
	   gbc.gridheight=gridheight;
	   gbc.weightx=weightx;
	   gbc.weighty=weighty;
	   gbc.fill=fill;
	   gbc.anchor=anchor;
	   container.add(component,gbc);
	 }
  }

  public void actionPerformed(ActionEvent e) {}

  public void stateChanged(ChangeEvent e) {}

  public void windowClosing(WindowEvent event) {}
  
  public void windowClosed(WindowEvent event) { }
  
  public void windowDeiconified(WindowEvent event) {}
  
  public void windowIconified(WindowEvent event) {}
  
  public void windowActivated(WindowEvent event) {}
  
  public void windowDeactivated(WindowEvent e) {}
  
  public void windowOpened(WindowEvent event) {}
  
  public void mouseClicked(MouseEvent mouse) {}

}
