import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

class nwchem_NWChem implements ActionListener, ChangeListener, WindowListener {
  
  Font defaultFont;
  int setnumber=0;
  JFileChooser chooser;
  ExtensionFilter pdbFilter;
  JFrame dialogFrame;
  String nwFile = new String();
  String command = new String();
  nwchem_Socket nwchemSocket;
  Runtime nwchemRun;
  Process nwchemProcess;

  public nwchem_NWChem(){

    //    super("NWChem running");

    //    defaultFont = new Font("Dialog", Font.BOLD,12);

    //    super.getContentPane().setLayout(new GridBagLayout());
    //    super.getContentPane().setForeground(Color.black);
    //    super.getContentPane().setBackground(Color.lightGray);
    //    super.getContentPane().setFont(defaultFont);
    //    super.addWindowListener(this);

    chooser = new JFileChooser("./");
    pdbFilter = new ExtensionFilter(".nw");
    chooser.setFileFilter(pdbFilter);
    dialogFrame = new JFrame();
    dialogFrame.setSize(300,400);
    chooser.showOpenDialog(dialogFrame);

    //    JPanel header = new JPanel();
    //    header.setLayout(new GridBagLayout());
    //    header.setForeground(Color.black);
    //    header.setBackground(Color.lightGray);
    //    addComponent(super.getContentPane(),header,0,0,2,1,1,1,
    //		 GridBagConstraints.NONE,GridBagConstraints.WEST);

    nwFile=chooser.getSelectedFile().toString();

    //    JLabel systemLabel = new JLabel(chooser.getSelectedFile().toString());
    //    addComponent(header,systemLabel,0,0,1,1,1,1,
    //    		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    //    systemLabel.setForeground(Color.black);
    
    //    JButton doneButton = new JButton("Done");
    //    addComponent(header,doneButton,5,0,1,1,1,1,
    //		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    //    doneButton.addActionListener(new ActionListener(){
    //      public void actionPerformed(ActionEvent e){ 
    //	setVisible(false); }});

    nwchemSocket = new nwchem_Socket();

    command="$NWCHEM_TOP/bin/$NWCHEM_TARGET/nwchem "+nwFile;
    System.out.println("Command: "+command);
    try{
      nwchemRun = Runtime.getRuntime();
      nwchemProcess = nwchemRun.exec(command);
      //      nwchemProcess.waitFor();
    } catch (Exception e) { e.printStackTrace(); };

    //    setLocation(25,225);	
    //    setSize(900,700);
    //    setVisible(true);

  }	

  void buildConstraints(GridBagConstraints gbc, int gx, int gy, int gw, int gh, 
			int wx, int wy){
    
    gbc.gridx = gx;
    gbc.gridy = gy;
    gbc.gridwidth = gw;
    gbc.gridheight = gh;
    gbc.weightx = wx;
    gbc.weighty = wy;
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
