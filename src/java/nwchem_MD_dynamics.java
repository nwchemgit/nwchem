import java.applet.Applet;
import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

public class nwchem_MD_dynamics extends JPanel implements ActionListener {
    
    ////////////////////////////////////////////////////////////////////////
    //	Dynamics panel variables
    ////////////////////////////////////////////////////////////////////////
    
  JLabel temperlabel;
  JLabel integrationlabel;
  JLabel equillabel;
  JLabel datalabel;
  JLabel timelabel;
  JLabel tsteplabel;
  JLabel ensemblelabel;
  JLabel vellabel;
  JLabel recordlabel;
  JLabel coordlabel;
  JLabel proplabel;
  JLabel rstlabel;
  JLabel synclabel;
  JLabel veloclabel;
  JLabel solvlabel;
  JLabel solxlabel;
  JLabel presslabel;
  JLabel relaxlabel;
  JLabel compresslabel;
  JLabel templabel;
  JLabel fillerlabel;
  JLabel updateslabel;
  JLabel pairslabel;
  JLabel centerlabel;
  JLabel motionlabel;
  JLabel rdflabel;
  
  IntegerField equilfield;
  IntegerField datafield;
  DoubleField timefield; 
  DoubleField tstepfield;
  IntegerField velfield;
  DoubleField temperfield;
  DoubleField pressfield;
  DoubleField relaxfield;
  DoubleField compressfield;
  DoubleField tempfield;
  DoubleField solventfield;
  DoubleField solutefield;
  IntegerField coordfield;
  IntegerField solxfield;
  IntegerField solvfield;
  IntegerField propfield;
  IntegerField rstfield;
  IntegerField syncfield;
  IntegerField velocfield;
  IntegerField centerfield;
  IntegerField motionfield;
  IntegerField rdffield;
  IntegerField pairsfield;
  
  JComboBox ensemblebox;
  JComboBox integrationbox;
  JComboBox eccebox;
  JComboBox relaxbox;

  String[] intstrings = {"leapfrog", "velocity Verlet"};
  String[] integrationstrings = {"leapfrog", "vverlet"};
  
  GridBagLayout dynamicslayout;
  GridBagConstraints dynamicsconstraints;
  Font defaultfont;
  
  private nwchem_Task task = null;
  
  
  public nwchem_MD_dynamics(nwchem_Task t) {
    
    defaultfont = new Font("Dialog", Font.BOLD,12);
    dynamicslayout = new GridBagLayout();
    setLayout(dynamicslayout);
    dynamicsconstraints = new GridBagConstraints();
    
    task = t;
    
    ////////////////////////////////////////////////////////////////////////////
    //	Setting up dynamics panel
    ///////////////////////////////////////////////////////////////////////////
    //JLabel integrationlabel
    
    setForeground(Color.black);
    setFont(defaultfont);
    
    integrationlabel = new JLabel("Integration: ");	
    integrationlabel.setFont(defaultfont);
    integrationlabel.setForeground(Color.black);
    addComponent(this,integrationlabel,0,0,2,1,15,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);


    integrationbox = new JComboBox(intstrings);
    integrationbox.setSelectedIndex(0);
    if(task.integration.equals("vverlet")) integrationbox.setSelectedIndex(1);
    addComponent(this,integrationbox,2,0,2,1,15,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    integrationbox.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	task.integration=integrationstrings[integrationbox.getSelectedIndex()];
      }});

    equillabel = new JLabel("equilibration:");
    equillabel.setFont(defaultfont);
    equillabel.setForeground(Color.black);
    addComponent(this,equillabel,4,0,2,1,15,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);

    equilfield = new IntegerField(task.equi,0,5);
    addComponent(this,equilfield,5,0,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    equilfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.equi=equilfield.getValue();};
      public void focusGained(FocusEvent e){}
    });

    datalabel = new JLabel("data:");
    datalabel.setFont(defaultfont);
    datalabel.setForeground(Color.black);
    addComponent(this,datalabel,4,1,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
	
    datafield = new IntegerField(task.data,1000, 5);
    addComponent(this,datafield,5,1,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    datafield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.data=datafield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    timelabel = new JLabel("time:");
    timelabel.setFont(defaultfont);
    timelabel.setForeground(Color.black);
    addComponent(this,timelabel,6,0,1,1,2,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    timefield = new DoubleField(task.time,0.0,5);
    addComponent(this,timefield,7,0,1,1,12,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    timefield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.time=timefield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    tsteplabel = new JLabel("time step:");
    tsteplabel.setFont(defaultfont);
    tsteplabel.setForeground(Color.black);
    addComponent(this,tsteplabel,6,1,1,1,14,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    tstepfield = new DoubleField(task.timestep,0.001,7);
    addComponent(this,tstepfield,7,1,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    tstepfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.timestep=tstepfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    ensemblelabel = new JLabel("Ensemble:");
    ensemblelabel.setFont(defaultfont);
    ensemblelabel.setForeground(Color.black);
    addComponent(this,ensemblelabel,0,2,2,1,15,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    String[] ensemstring = {"microcanonical", "canonical",
			    "isobaric", "isothermal-isobaric"};
    ensemblebox = new JComboBox(ensemstring);
    if(task.isotherm){
      if(task.isobar){ ensemblebox.setSelectedIndex(4); } 
      else { ensemblebox.setSelectedIndex(2); };
    } else {
      if(task.isobar){ ensemblebox.setSelectedIndex(3); } 
      else { ensemblebox.setSelectedIndex(0); };
    };
    addComponent(this,ensemblebox,2,2,2,1,15,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    ensemblebox.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	int i = ensemblebox.getSelectedIndex();
	task.isotherm=(i==1||i==3);
	task.isobar=(i==2||i==3);
	enableSelections();
      }});
    
    vellabel = new JLabel("Velocity reassignment frequency:");
    vellabel.setFont(defaultfont);
    vellabel.setForeground(Color.black);
    addComponent(this,vellabel,2,6,2,1,28,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    velfield = new IntegerField(task.freq_reas,0,8);
    addComponent(this,velfield,4,6,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    velfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.freq_reas=velfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    temperlabel = new JLabel("temp:");
    temperlabel.setFont(defaultfont);
    temperlabel.setForeground(Color.black);
    addComponent(this,temperlabel,5,6,2,1,14,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    temperfield = new DoubleField(task.temp_reas,298.15,7);
    addComponent(this,temperfield,6,6,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    temperfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.temp_reas=temperfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    presslabel = new JLabel("Pressure:");
    presslabel.setFont(defaultfont);
    presslabel.setForeground(Color.gray);
    addComponent(this,presslabel,3,5,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    pressfield = new DoubleField(task.press,102500.0,9);
    addComponent(this,pressfield,4,5,1,1,14,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    pressfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.press=pressfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    //    pressfield.setEnabled(false);
    
    relaxlabel = new JLabel("rlx time");
    relaxlabel.setFont(defaultfont);
    relaxlabel.setForeground(Color.gray);
    addComponent(this,relaxlabel,5,5,2,1,14,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    relaxfield = new DoubleField(task.prsrlx,0.4,5);
    addComponent(this,relaxfield,6,5,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    relaxfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.prsrlx=relaxfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
      //    relaxfield.setEnabled(false);
    
    compresslabel = new JLabel("compr:");
    compresslabel.setFont(defaultfont);
    compresslabel.setForeground(Color.gray);
    addComponent(this,compresslabel,7,5,1,1,14,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    compressfield = new DoubleField(task.compr,4.53E-10,6);
    addComponent(this,compressfield,9,5,1,1,1,1,
    		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    compressfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.compr=compressfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    dynamicslayout.setConstraints(compressfield,dynamicsconstraints);
    
    templabel = new JLabel("Temperature");
    templabel.setFont(defaultfont);
    templabel.setForeground(Color.gray);
    addComponent(this,templabel,3,4,2,1,14,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    tempfield = new DoubleField(task.temp,298.15,7);
    addComponent(this,tempfield,4,4,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    tempfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.temp=tempfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    String[] relaxstring = {"System rlx time", 
			    "Solvent/solute rlx times"};
    relaxbox = new JComboBox(relaxstring);
    if(task.isotherm_both){
      relaxbox.setSelectedIndex(1);
    } else {
      relaxbox.setSelectedIndex(0);
    };
    addComponent(this,relaxbox,5,4,2,1,21,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    relaxbox.addActionListener(this);
    relaxbox.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
	task.isotherm_both=(relaxbox.getSelectedIndex()==1);
	enableSelections();
      }});
    
    solventfield = new DoubleField(task.tmprlxw,0.1,5);
    addComponent(this,solventfield,7,4,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    solventfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.tmprlxw=solventfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    solutefield = new DoubleField(task.tmprlxs,0.1,5);
    addComponent(this,solutefield,8,4,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    solutefield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.tmprlxs=solutefield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    updateslabel = new JLabel("Updates:");
    updateslabel.setFont(defaultfont);
    updateslabel.setForeground(Color.black);
    addComponent(this,updateslabel,0,7,2,1,15,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    pairslabel = new JLabel("pairs:");
    pairslabel.setFont(defaultfont);
    pairslabel.setForeground(Color.black);
    addComponent(this,pairslabel,2,7,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    pairsfield = new IntegerField(task.freq_pair,1, 5);
    addComponent(this,pairsfield,3,7,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    pairsfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.freq_pair=pairsfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    centerlabel = new JLabel("center:");
    centerlabel.setFont(defaultfont);
    centerlabel.setForeground(Color.black);
    addComponent(this,centerlabel,4,7,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    centerfield = new IntegerField(task.freq_center,0, 5);
    addComponent(this,centerfield,5,7,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    centerfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.freq_center=centerfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    motionlabel = new JLabel("motion:");
    motionlabel.setFont(defaultfont);
    motionlabel.setForeground(Color.black);
    addComponent(this,motionlabel,6,7,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    motionfield = new IntegerField(task.freq_motion,0,5);
    addComponent(this,motionfield,7,7,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    motionfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.freq_motion=motionfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    rdflabel = new JLabel("rdf:");
    rdflabel.setFont(defaultfont);
    rdflabel.setForeground(Color.black);
    addComponent(this,rdflabel,8,7,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    rdffield = new IntegerField(task.freq_rdf,0,5);
    addComponent(this,rdffield,9,7,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    rdffield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.freq_rdf=rdffield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    recordlabel = new JLabel("Recording:");
    recordlabel.setFont(defaultfont);
    recordlabel.setForeground(Color.black);
    addComponent(this,recordlabel,0,8,2,1,15,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    coordlabel = new JLabel("coord:");
    coordlabel.setFont(defaultfont);
    coordlabel.setForeground(Color.black);
    addComponent(this,coordlabel,2,8,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    coordfield = new IntegerField(task.freq_coo,0,5);
    addComponent(this,coordfield,3,8,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    coordfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.freq_coo=coordfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    solxlabel = new JLabel("solute x");
    solxlabel.setFont(defaultfont);
    solxlabel.setForeground(Color.black);
    addComponent(this,solxlabel,4,8,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    solxfield = new IntegerField(task.freq_sco,0,5);
    addComponent(this,solxfield,5,8,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    solxfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.freq_sco=solxfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    rstlabel = new JLabel("restart");
    rstlabel.setFont(defaultfont);
    rstlabel.setForeground(Color.black);
    addComponent(this,rstlabel,6,8,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    rstfield = new IntegerField(task.freq_rst,0,5);
    addComponent(this,rstfield,7,8,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    rstfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.freq_rst=rstfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    synclabel = new JLabel("sync:");
    synclabel.setFont(defaultfont);
    synclabel.setForeground(Color.black);
    addComponent(this,synclabel,8,8,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    syncfield = new IntegerField(task.freq_syn,0,5);
    addComponent(this,syncfield,9,8,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    syncfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.freq_syn=syncfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    veloclabel = new JLabel("veloc:");
    veloclabel.setFont(defaultfont);
    veloclabel.setForeground(Color.black);
    addComponent(this,veloclabel,2,9,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    velocfield = new IntegerField(task.freq_vel,0,5);
    addComponent(this,velocfield,3,9,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    velocfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.freq_vel=velocfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    solvlabel = new JLabel("solute v");
    solvlabel.setFont(defaultfont);
    solvlabel.setForeground(Color.black);
    addComponent(this,solvlabel,4,9,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    solvfield = new IntegerField(task.freq_svl,0,5);
    addComponent(this,solvfield,5,9,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    solvfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.freq_svl=solvfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    proplabel = new JLabel("property");
    proplabel.setFont(defaultfont);
    proplabel.setForeground(Color.black);
    addComponent(this,proplabel,6,9,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    propfield = new IntegerField(task.freq_prp,0,5);
    addComponent(this,propfield,7,9,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    propfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.freq_prp=propfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    String[] eccestring = {"ecce"};
    eccebox = new JComboBox(eccestring);
    addComponent(this,eccebox,9,9,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    fillerlabel = new JLabel(" ");
    addComponent(this,fillerlabel,10,10,1,1,35,35,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    enableSelections();

  }
  
  void enableSelections(){
    if (ensemblebox.getSelectedIndex() == 0) {
      templabel.setForeground(Color.gray);
      tempfield.setEnabled(false);
      relaxbox.setSelectedIndex(0);
      relaxbox.setEnabled(false);
      solventfield.setEnabled(false);
      solutefield.setEnabled(false);
      presslabel.setForeground(Color.gray);
      pressfield.setEnabled(false);
      relaxlabel.setForeground(Color.gray);
      relaxfield.setEnabled(false);
      compressfield.setEnabled(false);
      compresslabel.setForeground(Color.gray);
    };
    if (ensemblebox.getSelectedIndex() == 1) {
      templabel.setForeground(Color.black);
      tempfield.setEnabled(true); 
      relaxbox.setEnabled(true); 
      solventfield.setEnabled(true); 
      presslabel.setForeground(Color.gray);
      pressfield.setEnabled(false);
      relaxlabel.setForeground(Color.gray);
      relaxfield.setEnabled(false);
      compressfield.setEnabled(false);
      compresslabel.setForeground(Color.gray);
    };
    if (ensemblebox.getSelectedIndex() == 2) {
      templabel.setForeground(Color.gray);
      tempfield.setEnabled(false);
      relaxbox.setSelectedIndex(0); 
      relaxbox.setEnabled(false);
      solventfield.setEnabled(false);
      solutefield.setEnabled(false);
      presslabel.setForeground(Color.black);
      pressfield.setEnabled(true);
      relaxlabel.setForeground(Color.black);
      relaxfield.setEnabled(true);
      compressfield.setEnabled(true);
      compresslabel.setForeground(Color.black);
    };
    if (ensemblebox.getSelectedIndex() == 3) {
      templabel.setForeground(Color.black);
      tempfield.setEnabled(true); 
      relaxbox.setEnabled(true); 
      solventfield.setEnabled(true); 
      presslabel.setForeground(Color.black);
      pressfield.setEnabled(true);
      relaxlabel.setForeground(Color.black);
      relaxfield.setEnabled(true);
      compressfield.setEnabled(true);
      compresslabel.setForeground(Color.black);
    };
    if(relaxbox.getSelectedIndex() == 0) {
      solutefield.setEnabled(false);
    };
    if(relaxbox.getSelectedIndex() == 1) {
      solutefield.setEnabled(true);
    };
  };

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
  
  public void actionPerformed(ActionEvent e){}
  
}
