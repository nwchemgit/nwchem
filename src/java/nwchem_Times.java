import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

class nwchem_Times extends JFrame implements ActionListener, ChangeListener, WindowListener {
  
  Font defaultFont;
  int setnumber=0;
  JFileChooser chooser;
  ExtensionFilter synchFilter;
  JFrame dialogFrame;

  Graph nodPlot = new Graph();
  Graph accPlot = new Graph();
  BufferedReader br;
  String card;

  int np;
  int nt;
  boolean first = true;

  int ndx[];
    double sdata[][][];
  double ndata[][];
  double tdata[];
  double time;

    int numData = 0;
    int curData = 0;

  public nwchem_Times(){

    super("Timing Analysis");

    defaultFont = new Font("Dialog", Font.BOLD,12);

    super.getContentPane().setLayout(new GridBagLayout());
    super.getContentPane().setForeground(Color.black);
    super.getContentPane().setBackground(Color.lightGray);
    super.getContentPane().setFont(defaultFont);
    super.addWindowListener(this);

    chooser = new JFileChooser("./");
    synchFilter = new ExtensionFilter(".tim");
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
    addComponent(header,systemLabel,0,0,5,1,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    systemLabel.setForeground(Color.black);
    
    JButton doneButton = new JButton("Done");
    addComponent(header,doneButton,5,1,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    doneButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	setVisible(false); }});
    
    JButton orderButton = new JButton("Order");
    addComponent(header,orderButton,3,1,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    orderButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	int k;
	for(int i=0; i<np-1; i++){
	  for(int j=i+1; j<np; j++){
	    if(ndata[ndx[i]][16]-ndata[ndx[i]][15]<ndata[ndx[j]][16]-ndata[ndx[j]][15]){
	      k=ndx[i]; ndx[i]=ndx[j]; ndx[j]=k;
	    };
	  };
	};
	addNodeData();
	nodPlot.fillPlot();
      }});
    
    JButton resetButton = new JButton("Reset");
    addComponent(header,resetButton,4,1,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    resetButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	for(int i=0; i<np; i++){ndx[i]=i;};
	addNodeData();
	nodPlot.fillPlot();
      }});
    
    JButton allButton = new JButton("All");
    addComponent(header,allButton,2,1,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    allButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	while(readData()){addAccuData();};
	addNodeData();
	nodPlot.fillPlot();
	accPlot.fillPlot();
      }});
    
    JButton nextButton = new JButton("Next");
    addComponent(header,nextButton,0,1,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    nextButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	  System.out.println("Next "+curData+" "+numData);
	  if(curData==numData) {readData(); addAccuData(); accPlot.fillPlot();} else {curData++; retrieveData();}; 
	  addNodeData(); nodPlot.fillPlot(); 
      }
    });
    
    JButton previousButton = new JButton("Previous");
    addComponent(header,previousButton,1,1,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    previousButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	  System.out.println("Prev "+curData+" "+numData);
	  if(curData==0) { readData(); } else { curData--; retrieveData();};
	  addNodeData(); nodPlot.fillPlot();
      }
    });

    setLocation(25,225);	
    setSize(900,700);
    setVisible(true);

    try{
      addComponent(header,nodPlot,0,2,5,5,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
      addComponent(header,accPlot,0,8,5,5,1,1,
    		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
      nodPlot.init();
      nodPlot.resize(700,300);
      nodPlot.setTitle("Wall Clock Time Decomposition per Processor");
      nodPlot.setXLabel("Processor");
      nodPlot.setBars(1.1,0.0);
      nodPlot.setMarksStyle("none");
      accPlot.init();
      accPlot.resize(700,300);
      accPlot.setTitle("Wall Clock Time Decomposition Accumulated");
      accPlot.setXLabel("Time");
      //      accPlot.setBars(1.1,0.0);
      //      accPlot.setMarksStyle("none");
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
      first=true;
      ndx = new int[np];
      sdata = new double[1000][np][17];
      ndata = new double[np][17];
      tdata = new double[17];
    } catch(Exception e) {e.printStackTrace();};

    for(int ip=0; ip<np; ip++){
      ndx[ip]=ip;
    };
  }

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
    } catch (Exception e) {return false;}
  };

    void storeData(){
	for(int i=1; i<17; i++){
	    for(int ip=0; ip<np; ip++){  
		sdata[numData][ip][i]=ndata[ip][i];
	    };
	};
	numData++; curData=numData;
    };

    void retrieveData(){
	for(int i=1; i<17; i++){
	    for(int ip=0; ip<np; ip++){  
		ndata[ip][i]=sdata[curData][ip][i];
	    };
	};
    };

  void addNodeData(){
    for(int i=0; i<17; i++){
	//	nodPlot.clear(i);
      try{ nodPlot.removeSet(i); } catch(Exception ee){};
    };
    for(int i=0; i<17; i++){
      for(int ip=0; ip<np; ip++){
	nodPlot.addData(i,ip,ndata[ndx[ip]][i],false,false);
      };
    };
  };

  void addAccuData(){
    for(int i=16; i>=0; i--){
      accPlot.addData(16-i,time,tdata[i],!first,false);
      //accPlot.addData(16-i,time,tdata[i],false,false);
    };
    first=false;
  };

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
