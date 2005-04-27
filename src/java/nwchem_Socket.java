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
  boolean first=true;

  public nwchem_Socket(){

    super("NWChem Socket");

    defaultFont = new Font("Dialog", Font.BOLD,12);

    super.getContentPane().setLayout(new GridBagLayout());
    super.getContentPane().setForeground(Color.black);
    super.getContentPane().setBackground(Color.lightGray);
    super.getContentPane().setFont(defaultFont);
    super.addWindowListener(this);

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
      validate();
      Socket client = server.accept();
      System.out.println("Client connected");
      validate();
      InputStreamReader isr = new InputStreamReader(client.getInputStream());
      BufferedReader is = new BufferedReader(isr);
      while(true){
	String inLine = is.readLine();
	if(inLine.length()>0){
	  //	  System.out.println(inLine.substring(0,5));
	  // System.out.println(inLine.substring(1,5));
	  //	  plotData(inLine);
	  if(inLine.startsWith("tETpV")) { plotData(inLine); };
	};
      }
    } catch (Exception s) { s.printStackTrace(); }
  }
  
  public void plotData(String card){
    double tim, tmp, prs, vol, ene;

    if(first){
      super.getContentPane().removeAll();
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
      tmpPlot.setTitle("Temperature");
      tmpPlot.setXLabel("Time");
      
      prsPlot.init();
      prsPlot.resize(100,100);
      prsPlot.setTitle("Pressure");
      prsPlot.setXLabel("Time");
      
      volPlot.init();
      volPlot.resize(100,100);
      volPlot.setTitle("Volume");
      volPlot.setXLabel("Time");
      
      enePlot.init();
      enePlot.resize(100,100);
      enePlot.setTitle("Energy");
      enePlot.setXLabel("Time");
      
      validate();
    };

    tim=Double.valueOf(card.substring(6,17)).doubleValue();
    ene=Double.valueOf(card.substring(18,29)).doubleValue();
    tmp=Double.valueOf(card.substring(30,41)).doubleValue();
    prs=Double.valueOf(card.substring(42,53)).doubleValue();
    vol=Double.valueOf(card.substring(54,65)).doubleValue();
    tmpPlot.addData(0,tim,tmp,!first,true);
    prsPlot.addData(0,tim,prs,!first,true);
    enePlot.addData(0,tim,ene,!first,true);
    volPlot.addData(0,tim,vol,!first,true); first=false;
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
