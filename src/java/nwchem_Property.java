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

    JFileChooser dataChooser;
    PrintWriter prop;
    ExtensionFilter dataFilter;
    
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
    JButton clearButton = new JButton("clear");
    JButton plotButton = new JButton("plot");
    JButton plotpButton = new JButton("plot+");
    JButton writeButton = new JButton("write");
    JButton newButton = new JButton("new");
    JButton appendButton = new JButton("append");
    JButton doneButton = new JButton("done");
    
    int plotnum = 10;
    int plotcur = 0;
    
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

    chooser = new JFileChooser("./");
    properFilter = new ExtensionFilter(".prp");
    chooser.setFileFilter(properFilter);
    dialogFrame = new JFrame();
    dialogFrame.setSize(300,400);
    chooser.showOpenDialog(dialogFrame);

    dataChooser = new JFileChooser("./");
    dataFilter = new ExtensionFilter(".dat");
    dataChooser.setFileFilter(dataFilter);

    JPanel header = new JPanel();
    header.setLayout(new GridBagLayout());
    header.setForeground(Color.black);
    header.setBackground(Color.lightGray);
    addComponent(super.getContentPane(),header,0,0,2,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.WEST);

    JLabel systemLabel = new JLabel(chooser.getSelectedFile().toString());

    addComponent(header,systemLabel,0,0,8,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    systemLabel.setForeground(Color.black);
    addComponent(header,sizeLabel,0,1,6,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);

    addComponent(header,xLabel,0,2,4,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    addComponent(header,yLabel,0,3,4,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    xLabel.setBackground(Color.yellow);

    addComponent(header,prp_plot,0,5,4,4,10,10,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    
    addComponent(header,propPane,4,5,4,4,10,10,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    pList.addMouseListener(this);

    addComponent(header,plotButton,4,3,1,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    addComponent(header,plotpButton,5,3,1,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    addComponent(header,writeButton,6,3,1,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    addComponent(header,clearButton,7,3,1,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    addComponent(header,newButton,4,2,1,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    addComponent(header,appendButton,5,2,1,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    addComponent(header,doneButton,7,2,1,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.WEST);

    doneButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	setVisible(false); }});

    prp_plot.init();
    prp_plot.resize(500,500);
    prp_plot.setTitle("Property Viewer");

    plotButton.addActionListener(this);
    plotpButton.addActionListener(this);
    writeButton.addActionListener(this);
    clearButton.addActionListener(this);

    appendButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	chooser.showOpenDialog(dialogFrame);
	appendData(chooser.getSelectedFile().toString());
      }});

    newButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	chooser.showOpenDialog(dialogFrame);
        numframes=0; propList.removeAllElements();
	offstep = 0.0;
	offtime = 0.0;
	appendData(chooser.getSelectedFile().toString());
      }});

    setLocation(25,225);	
    setSize(900,700);
    setVisible(true);

    appendData(chooser.getSelectedFile().toString());

    pList.setVisibleRowCount(15);

  }	

  void appendData(String dataFile){

    try{
      systemLabel.setText(dataFile);
      br = new BufferedReader(new FileReader(dataFile));
      card=br.readLine();
      numprop = Integer.parseInt(card.substring(1,7).trim());
      for(int i=0; i<numprop; i++){ 
	card=br.readLine();
	if(numframes==0) propList.addElement(card); 
      };
      while((card=br.readLine()) != null){ 
	  if(card.startsWith("frame")){card=br.readLine(); };
	  int num=0;
	  data = new double[numprop];
	  for(int i=0; i<numprop; i++){
	      if(num==4) { card=br.readLine(); num=0;};
	      data[i]=Double.valueOf(card.substring(num*12+1,(num+1)*12).trim()).doubleValue();
	      //
	      // following is patch to correct sign of data: substring appears to lose the minus sign
	      //
	      if(card.charAt(num*12)=='-'){data[i]=-Math.abs(data[i]);};
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
      systemLabel.setText(dataFile);
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
      if(e.getSource()==clearButton){
	  prp_plot.clear(true);
	  //	  for(int i=0; i<plotnum; i++){
	  //	      try{  prp_plot.clear(i); prp_plot.removeSet(i);} catch (Exception ee){};
	  //	  };  
	  prp_plot.fillPlot();
	  plotcur=0;
	  plotnum=0;
      };
    if(e.getSource()==plotButton||e.getSource()==plotpButton){
      boolean running = (e.getSource()==plotpButton);
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
	boolean first=true;
        double yave, ysum=0.0;
	for(int i=0; i<numframes; i++){
	  data = new double[numprop];
	  data=(double[])frames.elementAt(i);
	  prp_plot.addData(plotcur,data[xIndex-1],data[yIndex-1],!first,false);
	  if(running){
	    ysum=ysum+data[yIndex-1];
	    yave=ysum/(i+1);
	    prp_plot.addData(plotcur+1,data[xIndex-1],yave,!first,false);
	  };
	  first=false;
	};
	plotcur++; if(running){plotcur++;};
	plotnum++; if(running){plotnum++;};
	prp_plot.fillPlot();
	//	prp_plot.drawPlot(prp_plot.getGraphics(),true);
      };    
    };
    if(e.getSource()==writeButton){
      if(xIndex>0&&yIndex>0) {
	int ix1 = propList.getElementAt(xIndex-1).toString().indexOf("    ");
	int ix2 = propList.getElementAt(xIndex-1).toString().lastIndexOf("    ")+3;
	int lx =  propList.getElementAt(xIndex-1).toString().length();
	String xLab = new String();
	String yLab = new String();
	if(ix1>0&&ix2>ix1) {
	  if(ix2<lx-1) {
	    xLab=propList.getElementAt(xIndex-1).toString().substring(0,ix1)+" /"+propList.getElementAt(xIndex-1).toString().substring(ix2,lx);
	  } else {
	    xLab=propList.getElementAt(xIndex-1).toString().substring(0,ix1);
	  };
	} else {
	  xLab=propList.getElementAt(xIndex-1).toString();
	};
	ix1 = propList.getElementAt(yIndex-1).toString().indexOf("    ");
	ix2 = propList.getElementAt(yIndex-1).toString().lastIndexOf("    ")+3;
	lx =  propList.getElementAt(yIndex-1).toString().length();
	if(ix1>0&&ix2>ix1) {
	  if(ix2<lx-1) {
	    yLab=propList.getElementAt(yIndex-1).toString().substring(0,ix1)+" /"+propList.getElementAt(yIndex-1).toString().substring(ix2,lx);
	  } else {
	    yLab=propList.getElementAt(yIndex-1).toString().substring(0,ix1);
	  };
	} else {
	  yLab=propList.getElementAt(yIndex-1).toString();
	};
	try{ 
	    dataChooser.showOpenDialog(dialogFrame);
	    prop = new PrintWriter( new FileWriter(dataChooser.getSelectedFile().toString()));
	    prop.println("# "+xLab);
	    prop.println("# "+yLab);
	    for(int i=0; i<numframes; i++){
		data = new double[numprop];
		data=(double[])frames.elementAt(i);
		prop.println(data[xIndex-1]+" "+data[yIndex-1]);
	    };
	    prop.close();
	} catch (Exception ee) { System.out.println("Error writing to property file"); };
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
