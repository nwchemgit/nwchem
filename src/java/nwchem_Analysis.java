import java.util.*;
import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

class nwchem_Analysis extends JFrame implements ActionListener, ChangeListener, WindowListener, MouseListener {
  
    Font defaultFont;
    int setnumber=0;
    JFileChooser chooser;
    ExtensionFilter analFilter;
    JFrame dialogFrame;
  
    BufferedReader br;
    String card;

    JFileChooser dataChooser;
    PrintWriter prop;
    ExtensionFilter dataFilter;

    Graph anaPlot = new Graph();

    JLabel systemLabel = new JLabel();
    JButton doneButton = new JButton("done");
    JButton clearButton = new JButton("clear");
    JButton writeButton = new JButton("write");
    JRadioButton histButton = new JRadioButton("histo");
    
    double time;
    double data[][];
    int numdat;
    int iset=0;

    int nb, nh, nd, no, number, j;
    
    DefaultListModel propList = new DefaultListModel();
    JList pList = new JList(propList);
    JScrollPane propPane = new JScrollPane(pList);

    int currentSelection = 0;

    boolean histogram = false;
    
    public nwchem_Analysis(){
	
	super("Analysis Viewer");
	
	defaultFont = new Font("Dialog", Font.BOLD,12);
	
	super.getContentPane().setLayout(new GridBagLayout());
	super.getContentPane().setForeground(Color.black);
	super.getContentPane().setBackground(Color.lightGray);
	super.getContentPane().setFont(defaultFont);
	super.addWindowListener(this);
	
	chooser = new JFileChooser("./");
	analFilter = new ExtensionFilter(".ana");
	chooser.setFileFilter(analFilter);
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
	addComponent(header,systemLabel,0,0,10,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	systemLabel.setForeground(Color.black);
	
	BufferedReader br;

	try{
	    br = new BufferedReader(new FileReader(chooser.getSelectedFile().toString()));
	    String card;
	    card=br.readLine();
	    nb = Integer.parseInt(card.substring(1,10).trim());
	    nh = Integer.parseInt(card.substring(11,20).trim());
	    nd = Integer.parseInt(card.substring(21,30).trim());
	    no = Integer.parseInt(card.substring(31,40).trim());
            number=nb+nh+nd+no+1;
	    
	    data = new double[number][25000];
	    numdat=0;
	    
	    number--;

	    for(int i=0; i<number; i++){
		card=br.readLine();
		propList.addElement(card);
	    };

	    while((card=br.readLine()) != null){
		int k=1;
		for(int l=0; l<number+1; l=l+1){
		    if(k==61){card=br.readLine(); k=1;};
		    data[l][numdat]=Double.valueOf(card.substring(k,k+11)).doubleValue();
		    k=k+12;
		};
		numdat++;
	    };
	    br.close();
	} catch(Exception ee) {ee.printStackTrace();};

	anaPlot.setSize(800,600);
        anaPlot.init();

	addComponent(header,anaPlot,0,1,20,10,10,10,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	
	addComponent(header,histButton,21,2,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.SOUTHWEST);
	histButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){ 
		    for(int i=0; i<iset; i++){anaPlot.removeSet(i);}; iset=0; anaPlot.fillPlot();
		    histogram=!histogram;
		}});
    
	addComponent(header,propPane,21,3,3,1,1,10,
		     GridBagConstraints.NONE,GridBagConstraints.WEST);
       	pList.addMouseListener(this);
	
	addComponent(header,writeButton,21,4,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	writeButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){ 
		    if(currentSelection>0){write_ana(currentSelection);}; }});

	addComponent(header,clearButton,21,5,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	clearButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){ 
		    for(int i=0; i<iset; i++){anaPlot.removeSet(i);}; iset=0; anaPlot.fillPlot(); }});

	addComponent(header,doneButton,21,6,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	doneButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){ 
		    setVisible(false); }});
	/*
        number=0;
	for(int i=0; i<nb; i++){
	    propList.addElement("Distance "+(i+1));
	    number++;
	};

	for(int i=0; i<nh; i++){
	    number++;
	    propList.addElement("Angle "+(i+1));
	};

	for(int i=0; i<nd; i++){
	    number++;
	    propList.addElement("Torsion "+(i+1));
	};

	for(int i=0; i<no; i++){
	    number++;
	    propList.addElement("Improper "+(i+1));
	};
	*/

	setLocation(25,225);	
	setSize(900,700);
	setVisible(true);

	pList.setVisibleRowCount(15);

    }

    void plot_ana(int ndx){
        if(histogram){
	    double dmax = data[ndx][0];
	    double dmin = data[ndx][0];
	    for(int i=0; i<numdat; i++){
		if(data[ndx][i]>dmax){dmax=data[ndx][i];};
		if(data[ndx][i]<dmin){dmin=data[ndx][i];};
	    };
	    int id;
	    Double dd;
	    int hist[] = new int[1001];
	    for(int i=0; i<1001; i++){hist[i]=0;};
	    for(int i=0; i<numdat; i++){
		id =(int)  Math.round(((data[ndx][i]-dmin)*1000)/(dmax-dmin));
		hist[id]=hist[id]+1;
	    };
	    boolean first=true;
	    double dstep=(dmax-dmin)/1000.0;
	    for(int i=0; i<1000; i++){
		anaPlot.addData(iset,dmin+(i+0.5)*dstep,hist[i],!first,false); first=false;
	    };
	} else {
	    boolean first=true;
	    for(int i=0; i<numdat; i++){
		if(ndx>nb && i>0 && Math.abs(data[ndx][i-1]-data[ndx][i])>3.14){first=true;};
		anaPlot.addData(iset,data[0][i],data[ndx][i],!first,false); first=false;
	    };
	};
	anaPlot.fillPlot();
	iset++;
    };

    void write_ana(int ndx){
	try{ 
	    dataChooser.showOpenDialog(dialogFrame);
	    prop = new PrintWriter( new FileWriter(dataChooser.getSelectedFile().toString()));
	    System.out.println("Writing data "+ndx);
	    for(int i=0; i<numdat; i++){
		prop.println(data[0][i]+" "+data[ndx][i]);
	    };
	    prop.close();
	} catch (Exception ee) { System.out.println("Error writing to data file"); };
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

    public void actionPerformed(ActionEvent e) {
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
	    if(mouse.getSource()==pList){
		currentSelection=pList.getSelectedIndex()+1;
		plot_ana(currentSelection);
	    };
	};
    }

    public void mouseEntered(MouseEvent mouse){
    }

    public void mouseExited(MouseEvent mouse){
    }
  
}
