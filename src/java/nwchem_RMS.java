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
    Graph smrPlot = new Graph();

    JLabel systemLabel = new JLabel();
    JButton doneButton = new JButton("done");
    
    double time,rms1,rms2;
    int ndx[] = new int [1000];
    
    public nwchem_RMS(){
	
	super("RMS Viewer");
	
	defaultFont = new Font("Dialog", Font.BOLD,12);
	
	super.getContentPane().setLayout(new GridBagLayout());
	super.getContentPane().setForeground(Color.black);
	super.getContentPane().setBackground(Color.lightGray);
	super.getContentPane().setFont(defaultFont);
	super.addWindowListener(this);
	
	chooser = new JFileChooser();
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
	addComponent(header,systemLabel,0,0,10,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	systemLabel.setForeground(Color.black);
	
	addComponent(header,doneButton,0,0,2,1,1,10,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	doneButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){ 
		    setVisible(false); }});

	addComponent(header,rmsPlot,1,0,20,10,10,10,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	addComponent(header,smrPlot,21,0,20,10,10,10,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	try{
	    BufferedReader br = new BufferedReader(new FileReader(chooser.getSelectedFile().toString()));
	    String card;
	    card=br.readLine();
	    int numt=0;
	    boolean first=true;
	    while(!card.startsWith("Analysis")){
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
	    int m,n;
            int ndxs=-1;
	    int mprev=-1;
	    first=true;
	    while(!card.startsWith("Analysis")){
		rms1=Double.valueOf(card.substring(27,38)).doubleValue();
		System.out.println(card); 
		System.out.println(card.substring(11,16));
		System.out.println(card.substring(17,21));
		m=Integer.valueOf(card.substring(12,16)).intValue();
                if(m!=mprev){mprev=m; ndxs++;}; 
		n=Integer.valueOf(card.substring(18,21)).intValue();
		ndx[ndxs]=n;
		numa++;
		smrPlot.addData(0,numa,rms1,!first,false); first=false;
		card=br.readLine();
	    };
	    card=br.readLine();
	    numa=0;
	    first=true;
	    while(!card.startsWith("Analysis")){
		rms1=Double.valueOf(card.substring(12,23)).doubleValue();
		if(numa==0){
		    smrPlot.addData(0,0,rms1,!first,false); first=false;
		} else {
		    smrPlot.addData(0,ndx[numa-1],rms1,!first,false);
		}; 
		smrPlot.addData(0,ndx[numa],rms1,!first,false); numa++;
		card=br.readLine();
	    };
	    smrPlot.fillPlot();
	    br.close();
	    
	    System.out.println("Number is "+numt+","+numa);

	    System.out.println("---");

	    br = new BufferedReader(new FileReader(chooser.getSelectedFile().toString()));

	    card=br.readLine();

	    br.close();

	} catch(Exception ee) {ee.printStackTrace();};

	
	setLocation(25,225);	
	setSize(900,700);
	setVisible(true);

    }

    //    void plot_rms(String s){
    //	int ndx = Integer.parseInt(s);
    //	System.out.println("Plot index "+ndx);
    //	boolean first=true;
    //	for(int i=0; i<numdat; i++){
    //	    rmsPlot.addData(iset,data[0][i],data[ndx][i],!first,false); first=false;
    //	};
    //	rmsPlot.fillPlot();
    //	iset++;
    // };

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
