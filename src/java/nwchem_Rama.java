import java.util.*;
import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

class nwchem_Rama extends JFrame implements ActionListener, ChangeListener, WindowListener, MouseListener {
  
    Font defaultFont;
    int setnumber=0;
    JFileChooser chooser;
    ExtensionFilter ramFilter;
    JFrame dialogFrame;
  
    BufferedReader br;
    String card;

    Graph ramPlot = new Graph();

    JLabel systemLabel = new JLabel();
    JButton doneButton = new JButton("done");
    
    double time,rms1,rms2;
    
    public nwchem_Rama(){
	
	super("Ramachandran Viewer");
	
	defaultFont = new Font("Dialog", Font.BOLD,12);
	
	super.getContentPane().setLayout(new GridBagLayout());
	super.getContentPane().setForeground(Color.black);
	super.getContentPane().setBackground(Color.lightGray);
	super.getContentPane().setFont(defaultFont);
	super.addWindowListener(this);
	
	chooser = new JFileChooser("./");
	ramFilter = new ExtensionFilter(".ram");
	chooser.setFileFilter(ramFilter);
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

	ramPlot.setTitle("Ramachandran");
	ramPlot.setSize(700,700);
	addComponent(header,ramPlot,0,1,20,10,10,10,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	try{
	    BufferedReader br = new BufferedReader(new FileReader(chooser.getSelectedFile().toString()));
	    String card;
            int nump=0;
	    boolean first=true;
	    double phi, psi;
	    ramPlot.setMarksStyle("points");
	    while((card=br.readLine()) != null){
		time=Double.valueOf(card.substring(1,12)).doubleValue();
		card=br.readLine();
                first=true;
		while(!card.startsWith("      0")){
		    phi=Double.valueOf(card.substring(11,22)).doubleValue();
		    psi=Double.valueOf(card.substring(23,34)).doubleValue();
		    ramPlot.addData(nump,phi,psi,false,false); first=false;
		    card=br.readLine();
                    first=false;
		};
	    ramPlot.fillPlot();
	    nump++;
	    };
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
