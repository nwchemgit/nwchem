import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

class nwchem_Main extends JFrame implements ActionListener, WindowListener {

  // Menubar

  JMenuBar menubar = new JMenuBar();
  
  JMenuItem quit;
  JMenuItem nwchem, rasmol, gopmol;
  JMenuItem proper, anal, rms, free, synch, times, timing, rama;
  JMenuItem frgmnt, segmnt, param, seqnce, topol;

  // Status TextArea

  JTextArea statusarea;
  JScrollPane statuspane;

  nwchem_Main(){

    super("NWChem Computational Chemistry Software");
    super.setSize(400,200);
    super.setLocation(100,100);
    super.setBackground(Color.lightGray);
    super.setForeground(Color.blue);
    super.addWindowListener(this);
    super.setFont(new Font("Serif",Font.BOLD,10));

    GridBagLayout nwchemlayout = new GridBagLayout();
    super.getContentPane().setLayout(nwchemlayout);
    GridBagConstraints nwchemconstraints = new GridBagConstraints();

    // Setup MenuBar

    this.setJMenuBar(menubar);

    // Setup File Menu

    JMenu select = new JMenu("File");
    menubar.add(select);
      
    select.add(quit   = new JMenuItem("Quit"));
    quit.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	setVisible(false); System.exit(0); }});

    // Setup Run Menu

    JMenu run = new JMenu("Run");
    menubar.add(run);

    run.add(nwchem = new JMenuItem("NWChem"));
    nwchem.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	System.out.println("New NWChem Job");
	nwchem_Job job = new nwchem_Job();
	job.setVisible(true); 
      }});

    JMenu view = new JMenu("View");
    menubar.add(view);
    view.add(rasmol = new JMenuItem("Rasmol"));
    rasmol.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	System.out.println("RASMOL Coordinate Viewer");
	nwchem_Rasmol rasmol_View = new nwchem_Rasmol();
	rasmol_View.setVisible(true);
      }});
    view.add(gopmol = new JMenuItem("gOpenMol"));
    view.add(proper = new JMenuItem("property"));
    proper.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	System.out.println("Property Viewer");
	nwchem_Property proper_View = new nwchem_Property();
	proper_View.setVisible(true);
      }});
    view.add(anal = new JMenuItem("analysis"));
    anal.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	System.out.println("Analysis Viewer");
	nwchem_Analysis anal_View = new nwchem_Analysis();
	anal_View.setVisible(true);
      }});
    view.add(rms = new JMenuItem("rms"));
    rms.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	System.out.println("RMS Viewer");
	nwchem_RMS rms_View = new nwchem_RMS();
	rms_View.setVisible(true);
      }});
    view.add(rama = new JMenuItem("ramachandran"));
    rama.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	System.out.println("Ramachandran Viewer");
	nwchem_Rama rama_View = new nwchem_Rama();
	rama_View.setVisible(true);
      }});
    view.add(free   = new JMenuItem("free energy"));
    free.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	System.out.println("Free Energy Viewer");
	nwchem_Free free_View = new nwchem_Free();
	free_View.setVisible(true); 
      }});
    view.add(synch  = new JMenuItem("synchronization"));
    synch.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	System.out.println("Synchronization Viewer");
	nwchem_Synch synch_View = new nwchem_Synch();
	synch_View.setVisible(true); 
      }});
    view.add(times = new JMenuItem("times"));
    times.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	System.out.println("Timing Analysis Viewer");
	nwchem_Times times_View = new nwchem_Times();
	times_View.setVisible(true); 
      }});
    view.add(timing = new JMenuItem("timings"));
    timing.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	System.out.println("Timing Analysis Viewer");
	nwchem_Timing timing_View = new nwchem_Timing();
	timing_View.setVisible(true); 
      }});

    JMenu edit = new JMenu("Edit");
    menubar.add(edit);
    edit.add(frgmnt = new JMenuItem("fragment"));
    frgmnt.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	System.out.println("Fragment Editor");
	nwchem_Fragment frgmnt_View = new nwchem_Fragment();
	frgmnt_View.setVisible(true); 
      }});
    edit.add(segmnt = new JMenuItem("segment"));
    segmnt.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	System.out.println("Segment Editor");
	nwchem_Segment segmnt_View = new nwchem_Segment();
	segmnt_View.setVisible(true); 
      }});
    edit.add(param  = new JMenuItem("parameter"));
    param.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	System.out.println("Parameter Database Editor");
	nwchem_Param param_View = new nwchem_Param();
	param_View.setVisible(true); 
      }});
    edit.add(seqnce = new JMenuItem("sequence"));
    edit.add(topol  = new JMenuItem("topology"));

    buildConstraints(nwchemconstraints,0,0,1,1,1,1);
    nwchemconstraints.fill = GridBagConstraints.BOTH;
    nwchemconstraints.anchor = GridBagConstraints.CENTER;
    statusarea = new JTextArea("Status Window");
    statusarea.setEditable(false);
    statusarea.setBackground(Color.white);
    statuspane = new JScrollPane(statusarea);
    nwchemlayout.setConstraints(statuspane,nwchemconstraints);
    super.getContentPane().add(statuspane);

    setVisible(true);

  }
 
  void buildConstraints(GridBagConstraints gbc, int gx, int gy, 
			int gw,int gh, int wx, int wy){
    gbc.gridx = gx;
    gbc.gridy = gy;
    gbc.gridwidth = gw;
    gbc.gridheight = gh;
    gbc.weightx = wx;
    gbc.weighty = wy;}

  public void actionPerformed(ActionEvent e){}

  public void windowClosing(WindowEvent event){ setVisible(false); System.exit(0);}

  public void windowClosed(WindowEvent event) { System.exit(0); }
  
  public void windowDeiconified(WindowEvent event) {}
  
  public void windowIconified(WindowEvent event) {}
  
  public void windowActivated(WindowEvent event) {}
  
  public void windowDeactivated(WindowEvent event) {}
  
  public void windowOpened(WindowEvent event) {}

}
