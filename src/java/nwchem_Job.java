import java.util.*;
import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.lang.reflect.*;

class nwchem_Job extends JFrame implements ActionListener, ChangeListener, WindowListener{
  
  Font defaultFont;
  
  JPanel panel = new JPanel();

  int heap = 1;
  int stack = 48;
  int global = 24;
  boolean verify = false;

  private Vector tasks;
  int numberTasks = 0;
  int currentTask = -1;
  DefaultListModel taskModel = new DefaultListModel();
  JList taskList = new JList(taskModel);
  JScrollPane taskPane = new JScrollPane(taskList);

  String title = new String("Default Title For NWChem Job");
  String root = new String("test");
  TextField titlefield;
  TextField rtdbfield;
  IntegerField heapfield;
  IntegerField stackfield;
  IntegerField globalfield;
  JRadioButton verifybutton;
  JButton runJobButton;

  JTextArea taskarea;

  // Menubar

  JMenuBar menubar = new JMenuBar();
  JMenuItem test, quit;

  nwchem_Task task;
  nwchem_Input task_Input;

  PrintWriter nw;

  public nwchem_Job(){

    super("NWChem Job");

    defaultFont = new Font("Dialog", Font.BOLD,12);

    super.getContentPane().setLayout(new GridBagLayout());
    super.getContentPane().setForeground(Color.black);
    super.getContentPane().setBackground(Color.lightGray);
    super.getContentPane().setFont(defaultFont);
    super.addWindowListener(this);

    tasks = new Vector();

    // Setup MenuBar

    this.setJMenuBar(menubar);

    JMenu fileItem = new JMenu("File");
    menubar.add(fileItem);
    fileItem.add(quit  = new JMenuItem("Quit"));
    quit.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	setVisible(false);
      }});

    fileItem.add(test  = new JMenuItem("Test"));
    test.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	list_Task();
      }});
      
    addComponent(super.getContentPane(),panel,0,0,1,1,1.0,1.0,
		 GridBagConstraints.BOTH,GridBagConstraints.NORTHWEST);
    panel.setLayout(new GridBagLayout());

    JLabel titlelabel = new JLabel("Title ");
    addComponent(panel,titlelabel,0,0,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    titlefield = new TextField(title,30);
    addComponent(panel,titlefield,1,0,7,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.WEST);
    titlefield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){title=titlefield.getText();};
      public void focusGained(FocusEvent e){}
    });
    JLabel rootlabel = new JLabel("RTDB ");
    addComponent(panel,rootlabel,0,1,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    rtdbfield = new TextField(root,10);
    addComponent(panel,rtdbfield,1,1,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.WEST);
    rtdbfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){root=rtdbfield.getText();};
      public void focusGained(FocusEvent e){}
    });
    JLabel memorylabel = new JLabel("Memory:  ");
    addComponent(panel,memorylabel,2,1,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.EAST);
    JLabel heaplabel = new JLabel("heap ");
    addComponent(panel,heaplabel,3,1,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    heapfield = new IntegerField(heap,1,5);
    addComponent(panel,heapfield,4,1,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    heapfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){heap=heapfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    JLabel mb1label = new JLabel(" mb");
    addComponent(panel,mb1label,5,1,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    JLabel stacklabel = new JLabel("stack ");
    addComponent(panel,stacklabel,3,2,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    stackfield = new IntegerField(stack,24,5);
    addComponent(panel,stackfield,4,2,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    stackfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){stack=stackfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    JLabel mb2label = new JLabel(" mb");
    addComponent(panel,mb2label,5,2,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    JLabel globallabel = new JLabel("global ");
    addComponent(panel,globallabel,3,3,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    globalfield = new IntegerField(global,24,5);
    addComponent(panel,globalfield,4,3,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    globalfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){global=globalfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    JLabel mb3label = new JLabel(" mb");
    addComponent(panel,mb3label,5,3,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    verifybutton = new JRadioButton("verify");
    addComponent(panel,verifybutton,2,2,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.EAST);
    verifybutton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	verify=verifybutton.isSelected();
      }});
    JButton addTaskButton = new JButton("Add Task");
    addComponent(panel,addTaskButton,4,5,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.NORTH);
    addTaskButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	currentTask++; task=new nwchem_Task(); 
	tasks.addElement(task); nwchemInput();
      }});
    JButton insTaskButton = new JButton("Insert Task");
    addComponent(panel,insTaskButton,5,5,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.NORTH);
    insTaskButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	if(currentTask<0) currentTask=0; task=new nwchem_Task(); 
	tasks.insertElementAt(task,currentTask); nwchemInput();
      }});
    JButton writeJobButton = new JButton("Write Job");
    addComponent(panel,writeJobButton,4,6,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.NORTH);
    writeJobButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	writeJob(); 
      }});

    runJobButton = new JButton("Run Job");
    addComponent(panel,runJobButton,5,6,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.NORTH);
    runJobButton.addActionListener(this);

    addComponent(panel,taskPane,0,5,4,5,20,20,
		 GridBagConstraints.BOTH,GridBagConstraints.CENTER);
    JButton editTaskButton = new JButton("Edit");
    addComponent(panel,editTaskButton,4,7,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.NORTH);
    editTaskButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	currentTask=taskList.getSelectedIndex();
	task=(nwchem_Task)tasks.elementAt(currentTask);
	nwchemInput();
      }});
    JButton delTaskButton = new JButton("Delete");
    addComponent(panel,delTaskButton,5,7,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.NORTH);
    delTaskButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	try{task_Input.setVisible(false);} catch(Exception ee){};
	tasks.removeElementAt(taskList.getSelectedIndex());
	taskModel.removeAllElements();
	updateTaskList();
      }});
    JButton writeTaskButton = new JButton("Write");
    addComponent(panel,writeTaskButton,4,8,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.NORTH);
    JButton runTaskButton = new JButton("Run");
    addComponent(panel,runTaskButton,5,8,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.NORTH);
    JButton upTaskButton = new JButton("Move up");
    addComponent(panel,upTaskButton,4,9,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.NORTH);
    JButton downTaskButton = new JButton("Move down");
    addComponent(panel,downTaskButton,5,9,1,1,1,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.NORTH);

    setLocation(500,1);	
    setSize(800,400);

  }

  void list_Task(){
    System.out.println("Number of stored tasks is "+tasks.size());
    nwchem_Task t_temp = new nwchem_Task();
    for(int i=0; i<tasks.size(); i++){
      t_temp=(nwchem_Task)tasks.elementAt(i);
      System.out.println("Listed at "+i+" "+t_temp.theory+" "+t_temp.operation);
    };
    System.out.println("List done ");
  }

  public int num_Tasks(){ return tasks.size(); }

  public void updateTaskList(){
    nwchem_Task t_temp;
    taskModel.setSize(tasks.size());
    for(int i=0; i<tasks.size(); i++){
      //      nwchem_Task t_temp = new nwchem_Task();
      t_temp=(nwchem_Task)tasks.elementAt(i);
      System.out.println("task "+(i+1)+": "+t_temp.theory+" "+t_temp.operation);
      taskModel.setElementAt("task "+(i+1)+": "+t_temp.theory+" "+t_temp.operation,i);
    };
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
  public void actionPerformed(ActionEvent e){
    if(e.getSource()==runJobButton){
      nwchem_NWChem nn = new nwchem_NWChem();
    };
  }

  public void stateChanged(ChangeEvent e) {}

  public void windowClosing(WindowEvent event)
  { setVisible(false); System.exit(0);}

  public void windowClosed(WindowEvent event) { System.exit(0); }
  
  public void windowDeiconified(WindowEvent event) {}
  
  public void windowIconified(WindowEvent event) {}
  
  public void windowActivated(WindowEvent event) {}
  
  public void windowDeactivated(WindowEvent event) {}
  
  public void windowOpened(WindowEvent event) {}

  void nwchemInput(){
    try{task_Input.setVisible(false);} catch(Exception ee){};
    task_Input = new nwchem_Input(task);
    task_Input.setVisible(true);
    updateTaskList();
    task_Input.quit.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent ee){ updateTaskList(); }});
    task_Input.prepar.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent ee){ 
	task.theory="prepare"; task.operation=""; updateTaskList(); }});
    task_Input.md.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent ee){  
	task.theory="md"; if(task.operation.equals(""))
	  {task.operation="energy";}; 
	updateTaskList(); }});
    task_Input.analyz.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent ee){  
	task.theory="analyze"; task.operation=""; updateTaskList(); }});
    task_Input.scf.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent ee){   
	task.theory="scf"; if(task.operation.equals(""))
	  {task.operation="energy";}; 
	updateTaskList(); }});
    task_Input.dft.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent ee){   
	task.theory="dft"; if(task.operation.equals(""))
	  {task.operation="energy";}; 
	updateTaskList(); }});
    task_Input.esp.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent ee){   
	task.theory="esp"; task.operation="";
	updateTaskList(); }});
    task_Input.energy.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent ee){ 
	task.operation="energy"; 
	if(task.theory.equals("prepare") || task.theory.equals("analyze") 
	   || task.theory.equals("esp")) {task.operation="";};
	updateTaskList(); }});
    task_Input.optimize.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent ee){  
	task.operation="optimize";
	if(task.theory.equals("prepare") || task.theory.equals("analyze") 
	   || task.theory.equals("esp")) {task.operation="";}; 
	updateTaskList(); }});
    task_Input.dynamics.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent ee){  
	task.operation="dynamics";
	if(task.theory.equals("prepare") || task.theory.equals("analyze") 
	   || task.theory.equals("esp")) {task.operation="";}; 
	updateTaskList(); }});
    task_Input.thermodynamics.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent ee){   
	task.operation="thermodynamics";
	if(task.theory.equals("scf") || task.theory.equals("dft")) 
	  {task.operation="energy";}; 
	if(task.theory.equals("prepare") || task.theory.equals("analyze") 
	   || task.theory.equals("esp")) {task.operation="";}; 
	updateTaskList(); }});
  }

  void writeJob(){
    try{
      nw = new PrintWriter( new FileWriter(root+".nw"));
      nw.println("title; "+title+"\n");
      nw.print("memory "); 
      if(verify){nw.print("verify ");} else {nw.print("noverify ");};
      nw.println("heap "+heap+" mb stack "+stack+" mb global "+global+" mb\n");
      nw.println("start "+root+"\n");
      for(int i=0; i<tasks.size(); i++){
	task=(nwchem_Task)tasks.elementAt(i);
	if(task.theory.equals("md")){writeTaskMD();};
      };
      nw.close();
    } catch (Exception e) { System.out.println("Error writing to file"); };
  };

  void runJob(){
    writeJob();
    System.out.println("Running Job");
  };

  void writeTaskMD(){
    nw.println("md");

    nw.print(" set "+task.set);
    if(task.pset2||task.pset3){
      nw.print(" pset"); 
      if(task.pset2) nw.print(" 2");
      if(task.pset3) nw.print(" 3");
    }; nw.println();

    if(task.twin)
      nw.print(" cutoff short "+task.r_short+" long "+task.r_short);
    else nw.print(" cutoff "+task.r_short);
    if(task.qmmm) nw.println(" qmmm "+task.r_qmmm); else nw.println();

    if(task.interaction==1) nw.println(" polar first");
    if(task.interaction==2) nw.println(" polar scf "+task.polmax+" "+task.poltol);
    if(task.interaction==3) {
      nw.print(" pme grid ");
      if(task.gridy==task.gridz){
	if(task.gridx==task.gridy){ nw.print(task.gridx);} 
	else { nw.print(task.gridx+" "+task.gridy);};
      } else { nw.print(task.gridx+" "+task.gridy+" "+task.gridz);};
      nw.print(" order "+task.order);
      if(task.alpha>0.0) nw.print(" alpha "+task.alpha);
      if(task.nodpme>0) nw.print(" nodpme "+task.nodpme);
      if(task.fft>1) nw.print(" fft "+task.fft);
      nw.println();
    };
    
    if(task.operation.equals("optimize")){
      nw.println(" sd "+task.sd_iter+" init "+task.sd_init+" min "+task.sd_xmin+" max "+task.sd_xmax);
      nw.println(" cg "+task.cg_iter+" init "+task.cg_init+" min "+task.cg_xmin+" cy "+task.cg_cycle);
    };
    
    if(task.operation.equals("dynamics") || task.operation.equals("thermodynamics")){
      nw.print(" "+task.integration+" equil "+task.equi+" data "+task.data);
      if(task.operation.equals("thermodynamics")) nw.print(" over "+task.ti_over);
      nw.println(" time "+task.time+" step "+task.timestep);
      if(task.isotherm){
	nw.print(" isotherm "+task.temp+" trelax "+task.tmprlxw);
	if(task.isotherm_both) nw.print(" "+task.tmprlxs); nw.println();
      };
      if(task.isobar) nw.println(" isobar "+task.press+" trelax "+task.prsrlx+" compr "+task.compr);
      if(task.freq_reas>0) nw.println(" vreass "+task.freq_reas+" "+task.temp_reas);
    };

    if(task.operation.equals("thermodynamics")){
      nw.print(" "+task.ti_start+" "+task.ti_direct+" "+task.ti_steps+" of "+task.ti_max+
	       " error "+task.ti_error+" drift "+task.ti_drift+" factor "+task.ti_factor);
      if(task.ti_sss) nw.print(" sss delta "+task.ti_delta);
      if(task.ti_decomp) nw.print(" decomp"); nw.println();
    };

    if(task.shake_w||task.shake_s) nw.println(" shake "+task.shake_w_iter+" "+task.shake_s_iter+" "+
					    task.shake_w_toler+" "+task.shake_s_toler);
    if(!(task.shake_w&&task.shake_s)){
      nw.print(" noshake");
      if(!task.shake_w) nw.print(" solvent");
      if(!task.shake_s) nw.print(" solute"); nw.println();
    };

    if(task.fix||task.unfix) {
      if(task.fix_none||task.fix_all||task.fix_w||task.fix_s||task.fix_X){
	if(task.fix) { nw.print(" fix"); } else { nw.print(" unfix"); }
	if(task.fix_none) nw.print(" none");
	if(task.fix_all) { nw.print(" all"); } else {
	  if(task.fix_w) nw.print(" solvent");
	  if(task.fix_s) nw.print(" solute");
	  if(task.fix_X) nw.print(" non-H");
	};
	nw.println();
      };
    };

    if(task.prt_top||task.prt_step>0||task.prt_stat>0||task.prt_expect>0){
      nw.print(" print");
      if(task.prt_top){
	nw.print(" topol");
	if(task.prt_topn) nw.print(" nonbond");
	if(task.prt_topw) nw.print(" solvent");
	if(task.prt_tops) nw.print(" solute");
      };
      if(task.prt_step>0){
	nw.print(" step "+task.prt_step);
	if(task.prt_extra) nw.print(" extra");
	if(task.prt_energy) nw.print(" energy");
      };
      if(task.prt_stat>0) nw.print(" stat "+task.prt_stat);
      if(task.prt_expect>0) nw.print(" expect "+task.prt_expect);
      nw.println();
    };

    nw.print(" update pairs "+task.freq_pair);
    if(task.twin) nw.print(" long "+task.freq_long);
    nw.println(" center "+task.freq_center+" motion "+task.freq_motion);
    if(task.freq_rdf>0) nw.println(" update rdf "+task.freq_rdf+" range "+task.r_rdf+" bins "+task.num_rdf); 
    nw.print(" record rest "+task.freq_rst);  if(task.keep) nw.print(" keep");
    nw.print(" coord "+task.freq_coo+" scoor "+task.freq_sco+" veloc "+task.freq_vel+" svelo "+task.freq_svl);
    nw.print(" prop "+task.freq_prp+" sync "+task.freq_syn+" times "+task.freq_tim);
    if(task.operation.equals("optimize")) nw.print(" mind "+task.freq_mnd);
    if(task.operation.equals("thermodynamics")) {
      nw.print(" free "+task.freq_gib);
      if(task.cnv) nw.print(" cnv");
      if(task.fet) nw.print(" fet");
      if(task.acf) nw.print(" acf");
    };
    nw.println();

    if(task.distar){
      nw.print(" distar");
      if(task.draver) nw.print(" draver");
      if(task.drsscl!=1.0) nw.print(" scale "+task.drsscl);
      if(task.nfdrss!=0) nw.print(" after "+task.nfdrss);
      nw.println();
    };

    if(task.load_none||task.load_pair||task.load_size||task.load_reset){
      nw.print(" load");
      if(task.load_none) nw.print(" none");
      if(task.load_reset) nw.print(" reset");
      if(task.load_pair) nw.print(" pairs");
      if(task.load_size) {
	if(task.load_pair) nw.print(" "+task.load_num_pair);
	nw.print(" size "+task.load_factor);
      };
      nw.println();
    };

    if(task.npx>0) nw.print(" nodes "+task.npx+" "+task.npy+" "+task.npz);
    if(task.nbx>0) nw.print(" boxes "+task.nbx+" "+task.nby+" "+task.nbz);
    if(task.mad>6) nw.print(" extra "+task.mad);
    if(task.mwm>0) nw.print(" mwm "+task.mwm);
    if(task.msa>0) nw.print(" msa "+task.msa);
    if(task.npx>0||task.nbx>0||task.mad>6||task.mwm>0||task.msa>0) nw.println();

    if(task.Limit>0) nw.print(" limit "+task.Limit);
    if(task.test>0) nw.print(" test "+task.test);
    if(task.debug>0) nw.print(" debug "+task.debug);
    if(task.Limit>0||task.test>0||task.debug>0) nw.println();

    nw.println("end\n");
    nw.println("task md "+task.operation+"\n");
  }

  public void task_Fields(){
    try{
      Class cl = Class.forName("nwchem_Task");
      Field f[] = cl.getDeclaredFields();
      for(int j=0; j<f.length; j++){
	System.out.println("Field "+j+" = "+f[j].getName()+" ");
      };
    } catch(Exception ej){};
  }

}
