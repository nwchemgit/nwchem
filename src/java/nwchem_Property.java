import java.util.*;
import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

class nwchem_Property extends JFrame implements ActionListener, ChangeListener, WindowListener, MouseListener {
  
  Font defaultFont;
  int setnumber=0;
  JFileChooser chooser;
  ExtensionFilter properFilter;
  JFrame dialogFrame;
  
  BufferedReader br;
  String card;
  int numprop=0;
  int numframes=0;
  Vector frames = new Vector(10,5);
  double data[];
  
  DefaultListModel propList = new DefaultListModel();

  JList pList = new JList(propList);
  JScrollPane propPane = new JScrollPane(pList);
  boolean selectx = true;
  JLabel xLabel = new JLabel("x:                                   ");
  JLabel yLabel = new JLabel("y:                                   ");
  int xIndex=0, yIndex=0;
  Graph prp_plot = new Graph();
  JLabel systemLabel = new JLabel();
  JLabel sizeLabel = new JLabel("Number of frames is       ");
  JButton plotButton = new JButton("plot");
  JButton plotpButton = new JButton("plot+");
  JButton addButton = new JButton("add");
  JButton addpButton = new JButton("add+");
  JButton doneButton = new JButton("done");

  int type = 0;
  double time;
  double synt;
  double stpt;
  Graph synPlot = new Graph();
  Graph timPlot = new Graph();

  double offstep = 0.0;
  double offtime = 0.0;

  public nwchem_Property(){

    super("Property Viewer");

    defaultFont = new Font("Dialog", Font.BOLD,12);

    super.getContentPane().setLayout(new GridBagLayout());
    super.getContentPane().setForeground(Color.black);
    super.getContentPane().setBackground(Color.lightGray);
    super.getContentPane().setFont(defaultFont);
    super.addWindowListener(this);

    chooser = new JFileChooser();
    properFilter = new ExtensionFilter(".prp");
    chooser.setFileFilter(properFilter);
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

    addComponent(header,xLabel,0,1,2,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    addComponent(header,yLabel,0,2,2,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    xLabel.setBackground(Color.yellow);
    
    addComponent(header,propPane,0,3,1,1,10,10,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    pList.addMouseListener(this);

    addComponent(header,sizeLabel,2,0,6,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    addComponent(header,plotButton,2,1,1,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    addComponent(header,plotpButton,3,1,1,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    addComponent(header,addButton,4,1,1,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    addComponent(header,addpButton,5,1,1,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    addComponent(header,doneButton,6,1,1,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    doneButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	setVisible(false); }});

    addComponent(header,prp_plot,2,2,5,2,10,10,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);

    prp_plot.init();
    prp_plot.resize(500,500);
    prp_plot.setTitle("Property Viewer");

    plotButton.addActionListener(this);
    plotpButton.addActionListener(this);

    addButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	chooser.showOpenDialog(dialogFrame);
	appendData(chooser.getSelectedFile().toString());
      }});

    /*  
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
	    synPlot.addData(setnumber,time,synt,!first,true);
	    timPlot.addData(setnumber,time,stpt,!first,true); first=false;
	    if(type==1){card=br.readLine();};
	    for(int i=0; i<np; i++) {card=br.readLine();};
	  };
	  synPlot.fillPlot();
	  timPlot.fillPlot();
	  br.close();
	} catch(Exception ee) {ee.printStackTrace();};
      }});
      */

    setLocation(25,225);	
    setSize(900,700);
    setVisible(true);

    appendData(chooser.getSelectedFile().toString());

    /*
    try{
      br = new BufferedReader(new FileReader(chooser.getSelectedFile().toString()));
      card=br.readLine();
      numprop = Integer.parseInt(card.substring(1,7).trim());
      for(int i=0; i<numprop; i++){ card=br.readLine(); propList.addElement(card.substring(1,card.indexOf("  ").trim())); };
      numframes=0;
      while((card=br.readLine()) != null){
	int num=0;
	data = new double[numprop];
	for(int i=0; i<numprop; i++){
	  if(num==4) { card=br.readLine(); num=0; };
	  data[i]=Double.valueOf(card.substring(num*12+1,(num+1)*12)).doubleValue();
	  num++;
	};
        data[0]=data[0]+offstep;
        data[1]=data[1]+offstep;
	frames.add(numframes,data);
	//	frames.setElementAt(data,numframes);
	numframes++;
      };
      offstep=data[0];
      offtime=data[1];
      sizeLabel.setText("Number of frames is "+frames.size());
    } catch(Exception e) {e.printStackTrace();};
    */

    pList.setVisibleRowCount(15);

    /*
    try{
      addComponent(header,synPlot,0,1,5,5,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
      addComponent(header,timPlot,0,7,5,5,1,1,
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
      validate();
    } catch(Exception e) {e.printStackTrace();};

    try{
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
	synPlot.addData(0,time,synt,!first,true);
	timPlot.addData(0,time,stpt,!first,true); first=false;
	if(type==1){card=br.readLine();};
	for(int i=0; i<np; i++) {card=br.readLine();};
      };
      synPlot.fillPlot();
      timPlot.fillPlot();
      br.close();
    } catch(Exception e) {e.printStackTrace();};
    */
  }	

  void appendData(String dataFile){

    try{
      systemLabel.setText(dataFile);
      br = new BufferedReader(new FileReader(dataFile));
      card=br.readLine();
      numprop = Integer.parseInt(card.substring(1,7).trim());
      for(int i=0; i<numprop; i++){ 
	card=br.readLine();
	if(numframes==0) propList.addElement(card.substring(0,card.indexOf("   "))+" / "+card.substring(card.lastIndexOf(" "),card.length()).trim()); 
      };
      while((card=br.readLine()) != null){
	int num=0;
	data = new double[numprop];
	for(int i=0; i<numprop; i++){
	  if(num==4) { card=br.readLine(); num=0; };
	  data[i]=Double.valueOf(card.substring(num*12+1,(num+1)*12)).doubleValue();
	  num++;
	};
        data[0]=data[0]+offstep;
        data[1]=data[1]+offtime;
	frames.add(numframes,data);
	//	frames.setElementAt(data,numframes);
	numframes++;
      };
      offstep=data[0];
      offtime=data[1];
      sizeLabel.setText("Number of frames is "+frames.size());
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

  public void actionPerformed(ActionEvent e) {
    if(e.getSource()==plotButton||e.getSource()==plotpButton){
      //     e.getSource()==addButton||e.getSource()==addpButton){
      boolean running = (e.getSource()==plotpButton||e.getSource()==addpButton);
      //    boolean adddata = (e.getSource()==addButton||e.getSource()==addpButton);
      boolean adddata=false;
      if(xIndex>0&&yIndex>0) {
	int ix1 = propList.getElementAt(xIndex-1).toString().indexOf("    ");
	int ix2 = propList.getElementAt(xIndex-1).toString().lastIndexOf("    ")+3;
	int lx =  propList.getElementAt(xIndex-1).toString().length();
	String pLab = new String();
	if(ix1>0&&ix2>ix1) {
	  if(ix2<lx-1) {
	    pLab=propList.getElementAt(xIndex-1).toString().substring(0,ix1)+" /"+propList.getElementAt(xIndex-1).toString().substring(ix2,lx);
	  } else {
	    pLab=propList.getElementAt(xIndex-1).toString().substring(0,ix1);
	  };
	} else {
	  pLab=propList.getElementAt(xIndex-1).toString();
	};
	prp_plot.setXLabel(pLab);
	ix1 = propList.getElementAt(yIndex-1).toString().indexOf("    ");
	ix2 = propList.getElementAt(yIndex-1).toString().lastIndexOf("    ")+3;
	lx =  propList.getElementAt(yIndex-1).toString().length();
	if(ix1>0&&ix2>ix1) {
	  if(ix2<lx-1) {
	    pLab=propList.getElementAt(yIndex-1).toString().substring(0,ix1)+" /"+propList.getElementAt(yIndex-1).toString().substring(ix2,lx);
	  } else {
	    pLab=propList.getElementAt(yIndex-1).toString().substring(0,ix1);
	  };
	} else {
	  pLab=propList.getElementAt(yIndex-1).toString();
	};
	prp_plot.setYLabel(pLab);
	if(!adddata) { try{ prp_plot.removeSet(0); prp_plot.removeSet(1);} catch (Exception ee){}; };
	boolean first=true;
        double yave, ysum=0.0;
	for(int i=0; i<numframes; i++){
	  data = new double[numprop];
	  data=(double[])frames.elementAt(i);
	  prp_plot.addData(0,data[xIndex-1],data[yIndex-1],!first,false);
	  if(running){
	    ysum=ysum+data[yIndex-1];
	    yave=ysum/(i+1);
	    prp_plot.addData(1,data[xIndex-1],yave,!first,false);
	  };
	  first=false;
	};
	prp_plot.fillPlot();
	//	prp_plot.drawPlot(prp_plot.getGraphics(),true);
      };    
    };
  }

  public void stateChanged(ChangeEvent e) {}

  public void windowClosing(WindowEvent event) {}
  
  public void windowClosed(WindowEvent event) { }
  
  public void windowDeiconified(WindowEvent event) {}
  
  public void windowIconified(WindowEvent event) {}
  
  public void windowActivated(WindowEvent event) {}
  
  public void windowDeactivated(WindowEvent e) {}
  
  public void windowOpened(WindowEvent event) {}
  
  public void mouseClicked(MouseEvent mouse) {}

  public void mousePressed(MouseEvent mouse){}

  public void mouseReleased(MouseEvent mouse){
    if(mouse.getModifiers()==MouseEvent.BUTTON3_MASK){
      selectx=!selectx;
      if(selectx){
	xLabel.setBackground(Color.yellow);
	yLabel.setBackground(Color.lightGray);
	xLabel.setForeground(Color.blue);
	yLabel.setForeground(Color.darkGray);
      } else {
	xLabel.setBackground(Color.lightGray);
	yLabel.setBackground(Color.yellow);
	xLabel.setForeground(Color.darkGray);
	yLabel.setForeground(Color.blue);
      };
    }
    if(mouse.getModifiers()==MouseEvent.BUTTON1_MASK){
      if(mouse.getSource()==pList){
	if(selectx){
	  xIndex=pList.getSelectedIndex()+1;
	  xLabel.setText("x: "+pList.getSelectedValue());
	} else {
	  yIndex=pList.getSelectedIndex()+1;
	  yLabel.setText("y: "+pList.getSelectedValue());
	};
      };
    };
  }

  public void mouseEntered(MouseEvent mouse){
  }

  public void mouseExited(MouseEvent mouse){
  }
  
}
