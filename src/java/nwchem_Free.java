import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

class nwchem_Free extends JFrame implements ActionListener, ChangeListener, WindowListener {
    
    Font defaultFont;
    int setnumber=0;
    JFileChooser chooser;
    ExtensionFilter gibFilter;
    JFrame dialogFrame;
    int type = 0;
    double lambda;
    double dlambda;
    double deriv[] = new double[24];
    double free;
    double dfree;
    Graph gibPlot = new Graph();
    Graph cnvPlot = new Graph();
    
    public nwchem_Free(){
	
	super("Free Energy");
	
	defaultFont = new Font("Dialog", Font.BOLD,12);
	
	super.getContentPane().setLayout(new GridBagLayout());
	super.getContentPane().setForeground(Color.black);
	super.getContentPane().setBackground(Color.lightGray);
	super.getContentPane().setFont(defaultFont);
	super.addWindowListener(this);
	
	chooser = new JFileChooser();
	gibFilter = new ExtensionFilter(".gib");
	chooser.setFileFilter(gibFilter);
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
		    /*
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
		      gibPlot.addData(setnumber,time,synt,!first,true);
		      cnvPlot.addData(setnumber,time,stpt,!first,true); first=false;
		      if(type==1){card=br.readLine();};
		      for(int i=0; i<np; i++) {card=br.readLine();};
		      };
		      gibPlot.fillPlot();
		      cnvPlot.fillPlot();
		      br.close();
		      } catch(Exception ee) {ee.printStackTrace();};
		    */
		}});
	
	setLocation(25,225);	
	setSize(900,700);
	setVisible(true);
	
	try{
	    addComponent(header,gibPlot,0,1,5,5,1,1,
			 GridBagConstraints.NONE,GridBagConstraints.CENTER);
	    addComponent(header,cnvPlot,0,7,5,5,1,1,
			 GridBagConstraints.NONE,GridBagConstraints.CENTER);
	    gibPlot.init();
	    gibPlot.resize(500,300);
	    gibPlot.setNumSets(5);
	    gibPlot.setTitle("Free Energy");
	    gibPlot.setXLabel("Lambda");
	    cnvPlot.init();
	    cnvPlot.resize(500,300);
	    cnvPlot.setNumSets(5);
	    cnvPlot.setTitle("Free Energy");
	    cnvPlot.setXLabel("Time");
	    validate();
	} catch(Exception e) {e.printStackTrace();};
	
	try{
	    BufferedReader br = new BufferedReader(new FileReader(chooser.getSelectedFile().toString()));
	    String card;
	    free=0.0;
	    boolean first=true;
            int j,ndata;
	    double cnv[] = new double[10000];
            int mdata=10000;
            int numdat=0;
	    while((card=br.readLine()) != null){
		int nderiv = Integer.parseInt(card.substring(1,7).trim());
		ndata = Integer.parseInt(card.substring(8,14).trim());
                numdat=ndata;
                if(numdat>mdata){numdat=mdata;};
		if(first){ 
		    for(int i=0; i<numdat; i++){cnv[i]=0.0;};
		};
		lambda=Double.valueOf(card.substring(15,26)).doubleValue();
		dlambda=Double.valueOf(card.substring(27,38)).doubleValue();
		if(first){gibPlot.addData(0,lambda,free,!first,true);};
		first=false;
		for(int i=0; i<6; i++){
		    card=br.readLine();
		    deriv[4*i]=Double.valueOf(card.substring(1,20)).doubleValue();
		    deriv[4*i+1]=Double.valueOf(card.substring(21,40)).doubleValue();
		    deriv[4*i+2]=Double.valueOf(card.substring(41,60)).doubleValue();
		    deriv[4*i+3]=Double.valueOf(card.substring(61,80)).doubleValue();
		};
                j=0;
		for(int i=0; i<ndata; i=i+4){
		    card=br.readLine();
                    if(j<numdat){cnv[j]=cnv[j]+dlambda*Double.valueOf(card.substring(1,20)).doubleValue(); j++;};
                    if(j<numdat){cnv[j]=cnv[j]+dlambda*Double.valueOf(card.substring(21,40)).doubleValue(); j++;};
                    if(j<numdat){cnv[j]=cnv[j]+dlambda*Double.valueOf(card.substring(41,60)).doubleValue(); j++;};
                    if(j<numdat){cnv[j]=cnv[j]+dlambda*Double.valueOf(card.substring(61,80)).doubleValue(); j++;};
		};
		dfree=0.0;
		for(int i=0; i<24; i++){
		    dfree=dfree+deriv[i];
		};
		dfree=dfree/nderiv;
		free=free+dfree*dlambda;
		lambda=lambda+dlambda;
		gibPlot.addData(0,lambda,free,!first,false);
		card=br.readLine();
	    };
	    gibPlot.fillPlot();
	    br.close();
	    for(int i=1; i<numdat; i++){ cnv[i]=cnv[i]+cnv[i-1]; };
	    for(int i=1; i<numdat; i++){ cnv[i]=cnv[i]/(i+1); };
            first=true;
            for(int i=0; i<numdat; i++){ cnvPlot.addData(0,i+1,cnv[i],!first,false); first=false;};
	    cnvPlot.fillPlot();
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
