import java.util.*;
import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

class nwchem_RMS extends JFrame implements ActionListener, ChangeListener, WindowListener, MouseListener {
  
    Font defaultFont;
    int setnumber=0;
    JFileChooser chooser;
    ExtensionFilter rmsFilter;
    JFrame dialogFrame;
  
    BufferedReader br;
    String card;

    Graph rmsPlot = new Graph();
    Graph rmsaPlot = new Graph();
    Graph rmsrPlot = new Graph();
    Graph bfacaPlot = new Graph();
    Graph bfacrPlot = new Graph();

    JLabel systemLabel = new JLabel();
    JButton doneButton = new JButton("done");
    
    double time,rms1,rms2;
    
    public nwchem_RMS(){
	
	super("RMS Viewer 2");
	
	defaultFont = new Font("Dialog", Font.BOLD,12);
	
	super.getContentPane().setLayout(new GridBagLayout());
	super.getContentPane().setForeground(Color.black);
	super.getContentPane().setBackground(Color.lightGray);
	super.getContentPane().setFont(defaultFont);
	super.addWindowListener(this);
	
	chooser = new JFileChooser("./");
	rmsFilter = new ExtensionFilter(".rms");
	chooser.setFileFilter(rmsFilter);
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
	addComponent(header,systemLabel,2,0,10,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	systemLabel.setForeground(Color.black);
	
	addComponent(header,doneButton,0,0,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	doneButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){ 
		    setVisible(false); }});

	rmsPlot.init();
	rmsaPlot.init();
	rmsrPlot.init();
	bfacaPlot.init();
	bfacrPlot.init();
	rmsPlot.setTitle("RMS Deviation vs Time");
	rmsaPlot.setTitle("Atomic RMS Deviation");
	rmsrPlot.setTitle("Segment RMS Deviation");
	bfacaPlot.setTitle("Atomic B factor");
	bfacrPlot.setTitle("Segment B Factor");
	//	rmsrPlot.setBars(1.0,0.0);
	rmsPlot.setSize(500,300);
	rmsaPlot.setSize(350,300);
	rmsrPlot.setSize(350,300);
	bfacaPlot.setSize(350,300);
	bfacrPlot.setSize(350,300);
	addComponent(header,rmsPlot,0,1,20,10,10,10,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	addComponent(header,rmsaPlot,0,11,10,10,10,10,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	addComponent(header,rmsrPlot,11,11,10,10,10,10,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	addComponent(header,bfacaPlot,0,21,10,10,10,10,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	addComponent(header,bfacrPlot,11,21,10,10,10,10,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	try{
	    BufferedReader br = new BufferedReader(new FileReader(chooser.getSelectedFile().toString()));
	    String card;
	    card=br.readLine();
	    int numt=0;
	    boolean first=true;
	    while(!card.startsWith("analysis")){
		time=Double.valueOf(card.substring(1,12)).doubleValue();
		rms1=Double.valueOf(card.substring(13,24)).doubleValue();
		rms2=Double.valueOf(card.substring(25,36)).doubleValue();
		rmsPlot.addData(0,time,rms1,!first,false);
		rmsPlot.addData(1,time,rms2,!first,false); first=false;
		card=br.readLine();
	    };
	    rmsPlot.fillPlot();
	    card=br.readLine();
	    int numa=0;
	    first=true;
	    while(!card.startsWith("analysis")){
		rms1=Double.valueOf(card.substring(32,43)).doubleValue();
		rms2=Double.valueOf(card.substring(44,55)).doubleValue();
		if(!card.substring(30,31).equals("0")){
		    rmsaPlot.addData(0,numa,rms1,!first,false);
		    bfacaPlot.addData(0,numa,rms2,!first,false); first=false;
		    numa++;
		};
		card=br.readLine();
	    };
	    numa=0;
            int n;
	    first=true;
	    while((card=br.readLine()) != null){
		rms1=Double.valueOf(card.substring(12,23)).doubleValue();
		rms2=Double.valueOf(card.substring(24,35)).doubleValue();
		rmsrPlot.addData(0,numa,rms1,!first,false);
		bfacrPlot.addData(0,numa,rms2,!first,false); first=false; 
		numa++;
	    };
	    rmsaPlot.fillPlot();
	    rmsrPlot.fillPlot();
	    bfacaPlot.fillPlot();
	    bfacrPlot.fillPlot();
	    br.close();
	    
	} catch(Exception ee) {ee.printStackTrace();};

	
	setLocation(25,225);	
	setSize(900,700);
	setVisible(true);

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
