import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

class nwchem_Synch extends JFrame implements ActionListener, ChangeListener, WindowListener {
  
  Font defaultFont;
  int setnumber=0;
  JFileChooser chooser;
  ExtensionFilter synchFilter;
  JFrame dialogFrame;
  int type = 0;
  int node = 0;
  double time;
  double synt;
  double stpt;
  Graph synPlot = new Graph();
  Graph timPlot = new Graph();
  Graph nodPlot = new Graph();

  public nwchem_Synch(){

    super("Synchronization");

    defaultFont = new Font("Dialog", Font.BOLD,12);

    super.getContentPane().setLayout(new GridBagLayout());
    super.getContentPane().setForeground(Color.black);
    super.getContentPane().setBackground(Color.lightGray);
    super.getContentPane().setFont(defaultFont);
    super.addWindowListener(this);

    chooser = new JFileChooser("./");
    synchFilter = new ExtensionFilter(".syn");
    chooser.setFileFilter(synchFilter);
    dialogFrame = new JFrame();
    dialogFrame.setSize(300,400);
    chooser.showOpenDialog(dialogFrame);

    JPanel header = new JPanel();
    header.setLayout(new GridBagLayout());
    header.setForeground(Color.black);
    header.setBackground(Color.lightGray);
    addComponent(super.getContentPane(),header,0,0,2,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.WEST);

    JLabel systemLabel = new JLabel(chooser.getSelectedFile().toString());
    addComponent(header,systemLabel,0,0,1,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    systemLabel.setForeground(Color.black);
    
    JButton doneButton = new JButton("Done");
    addComponent(header,doneButton,5,0,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    doneButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	setVisible(false); }});
    
    JButton addButton = new JButton("Add");
    addComponent(header,addButton,6,0,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    addButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	chooser.showOpenDialog(dialogFrame);
	try{
	  setnumber++;
	  BufferedReader br = new BufferedReader(new FileReader(chooser.getSelectedFile().toString()));
	  String card;
	  card=br.readLine();
	  int np = Integer.parseInt(card.substring(1,5).trim());
	  boolean first=true;
	  while((card=br.readLine()) != null){
	    card=br.readLine();
	    type=Integer.parseInt(card.substring(1,5).trim());
	    time=Double.valueOf(card.substring(6,17)).doubleValue();
	    synt=Double.valueOf(card.substring(18,29)).doubleValue();
	    stpt=Double.valueOf(card.substring(30,41)).doubleValue();
	    synPlot.addData(setnumber,time,synt,!first,false);
	    timPlot.addData(setnumber,time,stpt,!first,false); first=false;
	    if(type==1){card=br.readLine();};
	    for(int i=0; i<np; i++) {
		card=br.readLine();
		if(i==0){
		    node=Integer.parseInt(card.substring(1,5).trim());
		    nodPlot.addData(1,time,node,false,false);
		};
		if(i==1){
		    node=Integer.parseInt(card.substring(1,5).trim());
		    nodPlot.addData(2,time,node,false,false);
		};
	    };
	    if(type==1){card=br.readLine();};
	    if(type==1){card=br.readLine();};
	    if(type==1){card=br.readLine();};
	  };
	  synPlot.fillPlot();
	  timPlot.fillPlot();
	  nodPlot.fillPlot();
	  br.close();
	} catch(Exception ee) {ee.printStackTrace();};
      }});

    setLocation(25,225);	
    setSize(900,700);
    setVisible(true);

    try{
      addComponent(header,synPlot,0,1,5,5,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
      addComponent(header,timPlot,0,7,5,5,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
      addComponent(header,nodPlot,0,13,5,5,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
      synPlot.init();
      synPlot.resize(500,300);
      synPlot.setNumSets(5);
      synPlot.setTitle("Accumulated Synchronization Time");
      synPlot.setXLabel("Time");
      timPlot.init();
      timPlot.resize(500,300);
      timPlot.setNumSets(5);
      timPlot.setTitle("MD Step Wall Clock Time");
      timPlot.setXLabel("Time");
      nodPlot.init();
      nodPlot.resize(500,300);
      nodPlot.setNumSets(5);
      nodPlot.setTitle("Nodes");
      nodPlot.setXLabel("Time");
      nodPlot.setMarksStyle("points");
      validate();
    } catch(Exception e) {e.printStackTrace();};

    try{
      BufferedReader br = new BufferedReader(new FileReader(chooser.getSelectedFile().toString()));
      String card;
      card=br.readLine();
      int np = Integer.parseInt(card.substring(1,5).trim());
      boolean first=true;
      while((card=br.readLine()) != null){
	  while(!card.startsWith("synchronization")){card=br.readLine();};
	card=br.readLine();
	type=Integer.parseInt(card.substring(1,5).trim());
	time=Double.valueOf(card.substring(6,17)).doubleValue();
	synt=Double.valueOf(card.substring(18,29)).doubleValue();
	stpt=Double.valueOf(card.substring(30,41)).doubleValue();
	synPlot.addData(0,time,synt,!first,false);
	timPlot.addData(0,time,stpt,!first,false); first=false;
	if(type==1){card=br.readLine();};
	for(int i=0; i<np; i++) {
	    card=br.readLine();
		if(i==0){
		    node=Integer.parseInt(card.substring(1,5).trim());
		    nodPlot.addData(1,time,node,false,false);
		};
		if(i==1){
		    node=Integer.parseInt(card.substring(1,5).trim());
		    nodPlot.addData(2,time,node,false,false);
		};
	};
	if(type==1){card=br.readLine();};
	if(type==1){card=br.readLine();};
	if(type==1){card=br.readLine();};
      };
      synPlot.fillPlot();
      timPlot.fillPlot();
      nodPlot.fillPlot();
      br.close();
    } catch(Exception e) {e.printStackTrace();};
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
