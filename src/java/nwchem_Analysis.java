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
	    
	    data = new double[number][10000];
	    numdat=0;
	    
	    while((card=br.readLine()) != null){
		int k=1;
		for(int l=0; l<number; l=l+1){
		    if(k==61){card=br.readLine(); k=1;};
		    data[l][numdat]=Double.valueOf(card.substring(k,k+11)).doubleValue();
		    k=k+12;
		};
		numdat++;
	    };
	    br.close();
	} catch(Exception ee) {ee.printStackTrace();};
	
	JButton tButton = new JButton("t");
	addComponent(header,tButton,0,1,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	
        number=0;
        int ix=0;
        int iy=2;
	for(int i=0; i<nb; i++){
	    number++;
	    JButton bButton = new JButton(Integer.toString(number));
	    addComponent(header,bButton,ix,iy,1,1,1,1,
			 GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	    bButton.addActionListener(new ActionListener(){
		    public void actionPerformed(ActionEvent e){ 
			plot_ana(e.getActionCommand()); }});
	    ix++; if(ix>11){ix=0; iy++;};
	};
	iy++;

	for(int i=0; i<nh; i++){
	    number++;
	    JButton bButton = new JButton(Integer.toString(number));
	    addComponent(header,bButton,ix,iy,1,1,1,1,
			 GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	    bButton.addActionListener(new ActionListener(){
		    public void actionPerformed(ActionEvent e){ 
			plot_ana(e.getActionCommand()); }});
	    ix++; if(ix>11){ix=0; iy++;};
	};
	iy++;

	for(int i=0; i<nd; i++){
	    number++;
	    JButton bButton = new JButton(Integer.toString(number));
	    addComponent(header,bButton,ix,iy,1,1,1,1,
			 GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	    bButton.addActionListener(new ActionListener(){
		    public void actionPerformed(ActionEvent e){ 
			plot_ana(e.getActionCommand()); }});
	    ix++; if(ix>11){ix=0; iy++;};
	};
	iy++;

	for(int i=0; i<no; i++){
	    number++;
	    JButton bButton = new JButton(Integer.toString(number));
	    addComponent(header,bButton,i,5,1,1,1,1,
			 GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	    bButton.addActionListener(new ActionListener(){
		    public void actionPerformed(ActionEvent e){ 
			plot_ana(e.getActionCommand()); }});
	    ix++; if(ix>11){ix=0; iy++;};
	};
	iy++;

	addComponent(header,plotButton,0,iy,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	plotButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){ 
		    anaPlot.fillPlot(); }});
	anaPlot.setSize(700,400);
	addComponent(header,anaPlot,1,iy,20,10,10,10,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	iy++;

	addComponent(header,clearButton,0,iy,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	clearButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){ 
		    for(int i=0; i<iset; i++){anaPlot.removeSet(i);}; iset=0; anaPlot.fillPlot(); }});
	iy++;

	addComponent(header,doneButton,0,iy,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	doneButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){ 
		    setVisible(false); }});
	
	setLocation(25,225);	
	setSize(900,700);
	setVisible(true);

    }

    void plot_ana(String s){
	int ndx = Integer.parseInt(s);
	boolean first=true;
	for(int i=0; i<numdat; i++){
	    if(ndx>nb && i>0 && Math.abs(data[ndx][i-1]-data[ndx][i])>3.14){first=true;};
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
