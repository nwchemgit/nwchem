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
    double free, freem, freep, freepm, ep2, ep3, tmp, ep2m, ep3m;
    double dfree, freeb, dfbias, dfreem;
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
	
	chooser = new JFileChooser("./");
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
	    gibPlot.setTitle("Free Energy vs Lambda");
	    gibPlot.setXLabel("Lambda");
	    gibPlot.setMarksStyle("dots",2);
	    cnvPlot.init();
	    cnvPlot.resize(500,300);
	    cnvPlot.setTitle("Free Energy Convergence");
	    cnvPlot.setXLabel("Time");
	    validate();
	} catch(Exception e) {e.printStackTrace();};
	
	try{
	    BufferedReader br = new BufferedReader(new FileReader(chooser.getSelectedFile().toString()));
	    String card;
	    free=0.0;
	    freem=0.0;
            freep=0.0;
	    freepm=0.0;
            freeb=0.0;
            lambda=0.0;
	    boolean first=true;
            int j,ndata;
            double lam;
            double freem=0.0;
	    double cnv[] = new double[10000];
	    double cnvm[] = new double[10000];
            int mdata=10000;
            int numdat=0;
            int ndec=0;
            int nsa=0;
	    while((card=br.readLine()) != null){
		int nderiv = Integer.parseInt(card.substring(0,8).trim());
		ndata = Integer.parseInt(card.substring(8,15).trim());
                numdat=ndata;
                if(numdat>mdata){numdat=mdata;};
		if(first){ 
		    for(int i=0; i<numdat; i++){cnv[i]=0.0;cnvm[i]=0.0;};
		};
		lam=Double.valueOf(card.substring(15,27)).doubleValue();
		dlambda=Double.valueOf(card.substring(27,39)).doubleValue();
                ndec=Integer.parseInt(card.substring(39,47).trim());
                nsa=Integer.parseInt(card.substring(47,54).trim());
		if(first){
		    gibPlot.addData(0,lambda,free,!first,true);
		    gibPlot.addData(1,lambda,freep,!first,true);
		    gibPlot.addData(2,lambda,freeb,!first,true);
		    lambda=lam;
		};
		first=false;
		for(int i=0; i<6; i++){
		    card=br.readLine();
		    deriv[4*i]=Double.valueOf(card.substring(0,20)).doubleValue();
		    deriv[4*i+1]=Double.valueOf(card.substring(20,40)).doubleValue();
		    deriv[4*i+2]=Double.valueOf(card.substring(40,60)).doubleValue();
		    deriv[4*i+3]=Double.valueOf(card.substring(60,80)).doubleValue();
		};
                j=0;
		for(int i=0; i<ndata; i=i+4){
		    card=br.readLine();
                    if(j<numdat){cnv[j]=cnv[j]+dlambda*Double.valueOf(card.substring(0,20)).doubleValue(); j++;};
                    if(j<numdat){cnv[j]=cnv[j]+dlambda*Double.valueOf(card.substring(20,40)).doubleValue(); j++;};
                    if(j<numdat){cnv[j]=cnv[j]+dlambda*Double.valueOf(card.substring(40,60)).doubleValue(); j++;};
                    if(j<numdat){cnv[j]=cnv[j]+dlambda*Double.valueOf(card.substring(60,80)).doubleValue(); j++;};
		};
                j=0;
		for(int i=0; i<ndata; i=i+4){
		    card=br.readLine();
                    if(j<numdat){cnvm[j]=cnvm[j]+dlambda*Double.valueOf(card.substring(0,20)).doubleValue(); j++;};
                    if(j<numdat){cnvm[j]=cnvm[j]+dlambda*Double.valueOf(card.substring(20,40)).doubleValue(); j++;};
                    if(j<numdat){cnvm[j]=cnvm[j]+dlambda*Double.valueOf(card.substring(40,60)).doubleValue(); j++;};
                    if(j<numdat){cnvm[j]=cnvm[j]+dlambda*Double.valueOf(card.substring(60,80)).doubleValue(); j++;};
		};
		dfree=0.0;
		for(int i=0; i<24; i++){
		    dfree=dfree+deriv[i];
		};
		dfreem=dfree-deriv[0]-deriv[12];
		dfree=dfree/nderiv;
		dfreem=dfreem/nderiv;
		free=free+dfree*dlambda;
		freem=freem+dfreem*dlambda;
		lambda=lambda+dlambda;
		gibPlot.addData(0,lambda,free,!first,false);
		gibPlot.addData(1,lambda,freem,!first,false);
		card=br.readLine();
		card=br.readLine();
		tmp=Double.valueOf(card.substring(0,20)).doubleValue();
		ep2=Double.valueOf(card.substring(20,40)).doubleValue();
		ep3=Double.valueOf(card.substring(40,60)).doubleValue();
		dfbias=Double.valueOf(card.substring(60,80)).doubleValue();
		card=br.readLine();
		ep2m=Double.valueOf(card.substring(0,20)).doubleValue();
		ep3m=Double.valueOf(card.substring(20,40)).doubleValue();
		freeb=freeb+dfbias*dlambda;
		if(freeb!=0.0) gibPlot.addData(2,lambda,freeb,!first,false);
		freep=freep+0.00831151*tmp*(Math.log(ep2)-Math.log(ep3));
		gibPlot.addData(3,lambda,freep,!first,false);
		freepm=freepm+0.00831151*tmp*(Math.log(ep2m)-Math.log(ep3m));
		gibPlot.addData(4,lambda,freepm,!first,false);
		if(ndec>0) {for(int i=0; i<6*nsa; i=i+4){ card=br.readLine(); };};
	    };
	    gibPlot.fillPlot();
	    br.close();
	    for(int i=1; i<numdat; i++){ cnv[i]=cnv[i]+cnv[i-1];cnvm[i]=cnvm[i]+cnvm[i-1]; };
	    for(int i=1; i<numdat; i++){ cnv[i]=cnv[i]/(i+1);cnvm[i]=cnvm[i]/(i+1); };
            first=true;
            for(int i=0; i<numdat; i++){ 
		cnvPlot.addData(0,i+1,cnv[i],!first,false);
		cnvPlot.addData(1,i+1,cnvm[i],!first,false); first=false;
	    };
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
