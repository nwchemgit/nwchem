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

    Graph anaPlot = new Graph();


    JLabel systemLabel = new JLabel();
    JButton doneButton = new JButton("done");
    JButton clearButton = new JButton("clear");
    JButton plotButton = new JButton("plot");
    
    double time;
    double data[][];
    int numdat;
    int iset=0;

    int nb, nh, nd, no, number, j;
    
    public nwchem_Analysis(){
	
	super("Analysis Viewer");
	
	defaultFont = new Font("Dialog", Font.BOLD,12);
	
	super.getContentPane().setLayout(new GridBagLayout());
	super.getContentPane().setForeground(Color.black);
	super.getContentPane().setBackground(Color.lightGray);
	super.getContentPane().setFont(defaultFont);
	super.addWindowListener(this);
	
	chooser = new JFileChooser();
	analFilter = new ExtensionFilter(".ana");
	chooser.setFileFilter(analFilter);
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
	    //	} catch(Exception ee) {ee.printStackTrace();};

	System.out.println("Number is "+number);

	data = new double[number][5000];
	numdat=0;

	System.out.println("---");

	//	try{
	    while((card=br.readLine()) != null){
		int k=1;
		for(int l=0; l<number; l=l+1){
		    if(k==49){card=br.readLine(); k=1;};
		    //		    System.out.println(card.substring(k,k+11));
		    data[l][numdat]=Double.valueOf(card.substring(k,k+11)).doubleValue();
		    k=k+12;
		};
		numdat++;
	    };
	    br.close();
	} catch(Exception ee) {ee.printStackTrace();};

	System.out.println("Number of frames is "+numdat);
	System.out.println(nb+" "+nh+" "+nd+" "+no);
	
	JButton tButton = new JButton("t");
	addComponent(header,tButton,0,1,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);

        number=0;
	for(int i=0; i<nb; i++){
	    number++;
	    JButton bButton = new JButton(Integer.toString(number));
	    addComponent(header,bButton,1,i+1,1,1,1,1,
			 GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	    bButton.addActionListener(new ActionListener(){
		    public void actionPerformed(ActionEvent e){ 
			plot_ana(e.getActionCommand()); }});
	};

	for(int i=0; i<nh; i++){
	    number++;
	    JButton bButton = new JButton(Integer.toString(number));
	    addComponent(header,bButton,2,i+1,1,1,1,1,
			 GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	    bButton.addActionListener(new ActionListener(){
		    public void actionPerformed(ActionEvent e){ 
			plot_ana(e.getActionCommand()); }});
	};

	for(int i=0; i<nd; i++){
	    number++;
	    JButton bButton = new JButton(Integer.toString(number));
	    addComponent(header,bButton,3,i+1,1,1,1,1,
			 GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	    bButton.addActionListener(new ActionListener(){
		    public void actionPerformed(ActionEvent e){ 
			plot_ana(e.getActionCommand()); }});
	};
	for(int i=0; i<no; i++){
	    number++;
	    JButton bButton = new JButton(Integer.toString(number));
	    addComponent(header,bButton,4,i+1,1,1,1,1,
			 GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	    bButton.addActionListener(new ActionListener(){
		    public void actionPerformed(ActionEvent e){ 
			plot_ana(e.getActionCommand()); }});
	};

	j=nb; if(nh>j){j=nh;}; if(nd>j){j=nd;}; if(no>j){j=no;}; j++; j++;

	addComponent(header,plotButton,0,j,2,1,1,10,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	plotButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){ 
		    setVisible(false); }});
	addComponent(header,clearButton,2,j,2,1,1,10,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	clearButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){ 
		    setVisible(false); }});
	addComponent(header,doneButton,4,j,2,1,1,10,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	doneButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){ 
		    setVisible(false); }});

	addComponent(header,anaPlot,6,1,20,10,10,10,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	
	setLocation(25,225);	
	setSize(900,700);
	setVisible(true);

    }

    void plot_ana(String s){
	int ndx = Integer.parseInt(s);
	System.out.println("Plot index "+ndx);
	boolean first=true;
	for(int i=0; i<numdat; i++){
	    anaPlot.addData(iset,data[0][i],data[ndx][i],!first,false); first=false;
	};
	anaPlot.fillPlot();
	iset++;
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
    }

    public void mouseEntered(MouseEvent mouse){
    }

    public void mouseExited(MouseEvent mouse){
    }
  
}
