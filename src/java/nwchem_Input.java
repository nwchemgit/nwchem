import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

class nwchem_Input extends JFrame implements ActionListener, ChangeListener, WindowListener {
  
  Font defaultFont;

  JTabbedPane inputpane;
  TextField systemName;
  JComboBox operation;
  nwchem_Task tsk;
  JLabel taskLabel = new JLabel();
  JPanel body = new JPanel();

  // Menubar

  JMenuBar menubar = new JMenuBar();
  JMenuItem quit, prepar, md, analyz, scf, dft, esp;
  JMenuItem energy, optimize, dynamics, thermodynamics;

  nwchem_MD md_Input;
  nwchem_Prepare prepare_Input;

  public nwchem_Input(nwchem_Task task){

    super("NWChem Task");

    defaultFont = new Font("Serif",Font.BOLD,12);

    super.getContentPane().setLayout(new GridBagLayout());
    super.getContentPane().setForeground(Color.black);
    super.getContentPane().setBackground(Color.lightGray);
    super.getContentPane().setFont(defaultFont);
    super.addWindowListener(this);

    tsk=task;

    // Setup MenuBar

    super.setJMenuBar(menubar);

    addComponent(super.getContentPane(),taskLabel,0,0,5,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);

    addComponent(super.getContentPane(),body,0,1,5,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.NORTHWEST);
    body.setLayout(new GridBagLayout());
    body.setBackground(Color.yellow);

    JMenu select = new JMenu("File");
    menubar.add(select);
    select.add(quit = new JMenuItem("Quit"));
    quit.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	setVisible(false);
      }});

    JMenu theory = new JMenu("Theory");
    menubar.add(theory);
    theory.add(prepar = new JMenuItem("prepare"));
    prepar.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	tsk.theory="prepare";
	setLabels();
	body.removeAll();
	prepare_Input = new nwchem_Prepare(tsk);
	addComponent(body,prepare_Input,0,0,1,1,1,1,
		     GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
	validate();
      }});
    theory.add(md     = new JMenuItem("md"));
    md.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	tsk.theory="md";
	setLabels();
	body.removeAll();
    	md_Input = new nwchem_MD(tsk);
	addComponent(body,md_Input,0,0,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
	validate();
      }});
    theory.add(analyz = new JMenuItem("analyze"));
    analyz.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	tsk.theory="analyze";
	setLabels();
	body.removeAll();
	validate();
      }});
    theory.add(scf    = new JMenuItem("scf"));
    scf.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	tsk.theory="scf";
	setLabels();
	body.removeAll();
	validate();
      }});
    theory.add(dft    = new JMenuItem("dft"));
    scf.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	tsk.theory="dft";
	setLabels();
	body.removeAll();
	validate();
      }});
    theory.add(esp    = new JMenuItem("esp"));
    esp.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	tsk.theory="esp";
	setLabels();
	body.removeAll();
	validate();
      }});

    JMenu operation = new JMenu("Operation");
    menubar.add(operation);
    operation.add(energy = new JMenuItem("energy"));
    energy.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	tsk.operation="energy";
	setLabels();
	validate();
      }});
    operation.add(optimize = new JMenuItem("optimize"));
    optimize.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	tsk.operation="optimize";
	setLabels();
	validate();
      }});
    operation.add(dynamics = new JMenuItem("dynamics"));
    dynamics.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	tsk.operation="dynamics";
	setLabels();
	validate();
      }});
    operation.add(thermodynamics = new JMenuItem("thermodynamics"));
    thermodynamics.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	tsk.operation="thermodynamics";
	setLabels();
	validate();
      }});

    setLocation(100,400);	
    setSize(800,400);

    body.removeAll();

    if(task.theory.equals("md")){
      if(tsk.operation.equals("")){tsk.operation="energy";};
      nwchem_MD md_Input = new nwchem_MD(tsk);
      addComponent(body,md_Input,0,0,1,1,1,1,
		   GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
      validate();
    };
    setLabels();

  }	

  void setLabels(){
    if(tsk.theory.equals("md")){ 
      taskLabel.setText(tsk.theory+" "+tsk.operation);
    } else if(tsk.theory.equals("scf") || tsk.theory.equals("dft")){
      if(tsk.operation.equals("thermodynamics")){tsk.operation="dynamics";};
      taskLabel.setText(tsk.theory+" "+tsk.operation);
    } else {
      tsk.operation=""; taskLabel.setText(tsk.theory);
    };
    taskLabel.setForeground(Color.black);
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
