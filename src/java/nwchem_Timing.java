import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

class nwchem_Timing extends JFrame implements ActionListener, ChangeListener, WindowListener, MouseListener {
    
    Font defaultFont;

    // Set up for a maximum of 5 different timings sets

    int numSets = 0;
    TimingSet[] timer = new TimingSet[5];

    Graph aPlot = new Graph();
    Graph bPlot = new Graph();
    Graph cPlot = new Graph();
    Graph dPlot = new Graph();
    
    DefaultListModel timeList = new DefaultListModel();
    JList tList = new JList(timeList);
    JScrollPane timePane = new JScrollPane(tList);
    
    JButton doneButton = new JButton("Done");
    JButton clearButton = new JButton("Clear");
    JButton sumButton = new JButton("Sum");
    JButton newButton = new JButton("New");
    JButton prevButton = new JButton("Prev");
    JButton nextButton = new JButton("Next");
    JButton orderButton = new JButton("Order");

    int current = -1;

    public nwchem_Timing(){
	
	super("Timing Analysis");
	
	defaultFont = new Font("Dialog", Font.BOLD,12);
	
	super.getContentPane().setLayout(new GridBagLayout());
	super.getContentPane().setForeground(Color.black);
	super.getContentPane().setBackground(Color.lightGray);
	super.getContentPane().setFont(defaultFont);
	super.addWindowListener(this);
	
	JPanel header = new JPanel();
	header.setLayout(new GridBagLayout());
	header.setForeground(Color.black);
	header.setBackground(Color.lightGray);
	addComponent(super.getContentPane(),header,0,0,2,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.WEST);
	
	addComponent(header,newButton,11,6,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.CENTER);
	newButton.addActionListener(this);

	addComponent(header,clearButton,12,6,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.CENTER);
	clearButton.addActionListener(this);

	addComponent(header,sumButton,13,6,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.CENTER);
	sumButton.addActionListener(this);

	addComponent(header,doneButton,14,6,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.CENTER);
	doneButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){ 
		    setVisible(false); }});
	
	addComponent(header,prevButton,11,7,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.CENTER);
	prevButton.addActionListener(this);

	addComponent(header,nextButton,12,7,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.CENTER);
	nextButton.addActionListener(this);

	addComponent(header,orderButton,13,7,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.CENTER);
	orderButton.addActionListener(this);
	
	addComponent(header,timePane,11,0,5,5,10,10,
		     GridBagConstraints.NONE,GridBagConstraints.WEST);
	tList.addMouseListener(this);
	tList.setVisibleRowCount(15);
	
	setLocation(25,225);	
	setSize(1400,800);
	setVisible(true);

	timer[0] = new TimingSet(numSets,timeList); numSets=1;

	// set up graphs

	try{
	    addComponent(header,aPlot,0,0,5,5,1,1,
			 GridBagConstraints.NONE,GridBagConstraints.CENTER);
	    addComponent(header,bPlot,5,0,5,5,1,1,
			 GridBagConstraints.NONE,GridBagConstraints.CENTER);
	    addComponent(header,cPlot,0,5,5,5,1,1,
			 GridBagConstraints.NONE,GridBagConstraints.CENTER);
	    addComponent(header,dPlot,5,5,5,5,1,1,
			 GridBagConstraints.NONE,GridBagConstraints.CENTER);
	    aPlot.init();
	    aPlot.resize(700,300);
	    aPlot.setTitle("Processor Wall Clock Time");
	    aPlot.setXLabel("Processor");
	    aPlot.setBars(1.1,0.0);
	    aPlot.setMarksStyle("none");
	    bPlot.init();
	    bPlot.resize(700,300);
	    bPlot.setTitle("Accumulated Wall Clock Time");
	    bPlot.setXLabel("Time");
	    cPlot.init();
	    cPlot.resize(700,300);
	    cPlot.setTitle("Accumulated Wall Clock Time");
	    cPlot.setXLabel("Time");
	    dPlot.init();
	    dPlot.resize(700,300);
	    dPlot.setTitle("Accumulated Wall Clock Time");
	    dPlot.setXLabel("Time");
	    validate();
	} catch(Exception e) {e.printStackTrace();};

	/*
	    br = new BufferedReader(new FileReader(chooser.getSelectedFile().toString()));
	    String card;
	    card=br.readLine();
	    numProcs1 = Integer.parseInt(card.substring(1,5).trim());
	    ntimes = Integer.parseInt(card.substring(6,10).trim());
	    numFrames1=0;
	    while((card=br.readLine()) != null){
		if(card.startsWith("timings")){
		    numFrames1++;
		};
	    };
	    br.close();
	    System.out.println("Number of frames is "+numFrames1);
	    dataSet1 = new double[numFrames1][numProcs1][ntimes];
	    timeSet1 = new double[numFrames1];
	    accuData1 = new double[numFrames1];
	    br = new BufferedReader(new FileReader(chooser.getSelectedFile().toString()));
	    card=br.readLine();
	    for(int itime=0; itime<ntimes; itime++){
		card=br.readLine(); System.out.println(card);
		timeList.addElement(card);
	    };
	    for(int iframe=0; iframe<numFrames1; iframe++){
		card=br.readLine(); card=br.readLine();
		timeSet1[iframe]=Double.valueOf(card.substring(0,12)).doubleValue();
		for(int iproc=0; iproc<numProcs1; iproc++){
		    int icount=0;
		    for(int itime=0; itime<ntimes; itime++){
			if(icount==0) card=br.readLine();
			dataSet1[iframe][iproc][itime]=Double.valueOf(card.substring(icount*7,(icount+1)*7)).doubleValue();
			icount++; if(icount>9) icount=0;
		    };
		};
	    };
	    addComponent(header,aPlot,0,2,5,5,1,1,
			 GridBagConstraints.NONE,GridBagConstraints.CENTER);
	    addComponent(header,bPlot,0,10,5,5,1,1,
			 GridBagConstraints.NONE,GridBagConstraints.CENTER);
	    aPlot.init();
	    aPlot.resize(700,300);
	    aPlot.setTitle("Processor Wall Clock Time");
	    aPlot.setXLabel("Processor");
	    aPlot.setBars(1.1,0.0);
	    aPlot.setMarksStyle("none");
	    bPlot.init();
	    bPlot.resize(700,300);
	    bPlot.setTitle("Accumulated Wall Clock Time");
	    bPlot.setXLabel("Time");
	    //      bPlot.setBars(1.1,0.0);
	    //      bPlot.setMarksStyle("none");
	    validate();
	} catch(Exception e) {e.printStackTrace();};
	
	*/

    /*
    try{
      addComponent(header,aPlot,0,2,5,5,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
      addComponent(header,bPlot,0,8,5,5,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
      aPlot.init();
      aPlot.resize(700,300);
      aPlot.setTitle("Wall Clock Time Decomposition per Processor");
      aPlot.setXLabel("Processor");
      aPlot.setBars(1.1,0.0);
      aPlot.setMarksStyle("none");
      bPlot.init();
      bPlot.resize(700,300);
      bPlot.setTitle("Wall Clock Time Decomposition Accumulated");
      bPlot.setXLabel("Time");
      //      bPlot.setBars(1.1,0.0);
      //      bPlot.setMarksStyle("none");
      validate();
    } catch(Exception e) {e.printStackTrace();};

    try{
      br = new BufferedReader(new FileReader(chooser.getSelectedFile().toString()));
      String card;
      card=br.readLine();
      np = Integer.parseInt(card.substring(1,5).trim());
      nt = Integer.parseInt(card.substring(6,10).trim());
      for(int it=0; it<nt; it++){
	  card=br.readLine();
      };
      firstElement1=true;
      ndx = new int[np];
      dataSet1 = new double[1000][np][17];
      ndata = new double[np][17];
      tdata = new double[17];
    } catch(Exception e) {e.printStackTrace();};

    for(int ip=0; ip<np; ip++){
      ndx[ip]=ip;
    };
    */

  }

    /*
  boolean readData(){
    try{
      card=br.readLine();
      card=br.readLine();
      time=Double.valueOf(card.substring(1,12)).doubleValue();
      int num = 0;
      for(int ip=0; ip<np; ip++){
	card=br.readLine();
	for(int i=0; i<17; i++){
	  ndata[ip][i]=Double.valueOf(card.substring(i*7+1,i*7+7)).doubleValue();
	};
      };
      for(int i=0; i<17; i++){
	tdata[i]=0.0;
	for(int ip=0; ip<np; ip++){
	  tdata[i]=tdata[i]+ndata[ip][i];
	};
      };
      for(int i=1; i<17; i++){
	tdata[i]=tdata[i]+tdata[i-1];
	for(int ip=0; ip<np; ip++){  
	  ndata[ip][i]=ndata[ip][i]+ndata[ip][i-1];
	};
      };
      storeData();
      return true;
    } catch (Exception e) {return false;};
  };

    void storeData(){
	for(int i=1; i<17; i++){
	    for(int ip=0; ip<np; ip++){  
		dataSet1[numData][ip][i]=ndata[ip][i];
	    };
	};
	numData++; curData=numData;
    };

    void retrieveData(){
	for(int i=1; i<17; i++){
	    for(int ip=0; ip<np; ip++){  
		ndata[ip][i]=dataSet1[curData][ip][i];
	    };
	};
    };

  void addNodeData(){
    for(int i=0; i<17; i++){
	//	aPlot.clear(i);
      try{ aPlot.removeSet(i); } catch(Exception ee){};
    };
    for(int i=0; i<17; i++){
      for(int ip=0; ip<np; ip++){
	aPlot.addData(i,ip,ndata[ndx[ip]][i],false,false);
      };
    };
  };

  void addAccuData(){
    for(int i=16; i>=0; i--){
      bPlot.addData(16-i,time,tdata[i],!firstElement1,false);
      //bPlot.addData(16-i,time,tdata[i],false,false);
    };
    firstElement1=false;
  };

    */

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
      if(e.getSource()==newButton){
	  timer[numSets]= new TimingSet(numSets,timeList); numSets++;
      };
      if(e.getSource()==sumButton){
	  if(numSets==1){
	      bPlot.clear(true); 
	      bPlot.setTitle("Accumulated Wall Clock Time Set 1");
	      bPlot.setXLabel("Time");
	      timer[0].sumPlot(bPlot);
	      bPlot.fillPlot();
	  } else {
	      aPlot.clear(true);  
	      aPlot.setTitle("Accumulated Wall Clock Time Set 1");
	      aPlot.setXLabel("Time");
	      timer[0].sumPlot(aPlot);
	      aPlot.fillPlot();
	      bPlot.clear(true);  
	      bPlot.setTitle("Accumulated Wall Clock Time Set 2");
	      bPlot.setXLabel("Time");
	      timer[1].sumPlot(bPlot);
	      bPlot.fillPlot();
	  };
      };
      if(e.getSource()==clearButton){
	  bPlot.clear(true); bPlot.fillPlot();
      };
      if(e.getSource()==nextButton){
	  current++;
	  cPlot.clear(true);
	  cPlot.setBars(1.1,0.0);
	  cPlot.setMarksStyle("none");
	  cPlot.setConnected(false);
	  timer[0].nodPlot(current,cPlot); cPlot.fillPlot();
      };
      if(e.getSource()==prevButton){
	  current--;
	  cPlot.clear(true);
	  cPlot.setBars(1.1,0.0);
	  cPlot.setMarksStyle("none");
	  cPlot.setConnected(false);
	  timer[0].nodPlot(current,cPlot); cPlot.fillPlot();
      };
      if(e.getSource()==orderButton){
	  timer[0].order=!timer[0].order;
	  cPlot.clear(true);
	  cPlot.setBars(1.1,0.0);
	  cPlot.setMarksStyle("none");
	  cPlot.setConnected(false);
	  timer[0].nodPlot(current,cPlot); cPlot.fillPlot();
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
    if(mouse.getModifiers()==MouseEvent.BUTTON1_MASK){
      if(mouse.getSource()==tList){
	  int index=tList.getSelectedIndex(); 
	  bPlot.clear(true);  
	  bPlot.setTitle(tList.getSelectedValue().toString());
	  for(int i=0; i<numSets; i++){
	      timer[i].plot(index,i,bPlot); 
	  };
	  bPlot.fillPlot(); 
	  //super.getContentPane().validate();
      };
    };
  }

  public void mouseEntered(MouseEvent mouse){
  }

  public void mouseExited(MouseEvent mouse){
  }
  
}
