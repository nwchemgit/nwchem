import java.applet.Applet;
import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.lang.*;

class nwchem_MD_common extends JPanel {
  
  ///////////////////////////////////////////////////////////////////////
  //	Variables that all the panels share.
  //////////////////////////////////////////////////////////////////////
  Font defaultfont;
  Font disabledfont;
  
  JComboBox interactbox;
  JComboBox cutbox;
  JComboBox fixbox;
  
  JLabel interactlabel;
  JLabel setlabel;
  JLabel maxlabel;
  JLabel tolerlabel;
  JLabel gridlabel;
  JLabel alphalabel;
  JLabel orderlabel;
  JLabel nodeslabel;
  JLabel fftlabel;
  JLabel fillerlabel;
  JLabel cutlabel;
  JLabel shortlabel;
  JLabel longlabel;
  JLabel freqlabel;
  JLabel shakelabel;
  JLabel solventlabel;
  JLabel solutelabel;
  JLabel maxit1label;
  JLabel maxit2label;
  JLabel tol1label;
  JLabel tol2label;
  JLabel distlabel;
  JLabel scale1label;
  JLabel scale2label;
  JLabel fixlabel;
  JLabel printlabel;
  JLabel steplabel;
  JLabel statlabel;
  JLabel expectlabel;
  JLabel calculationlabel;
  
  ButtonGroup  setgroup;
  ButtonGroup  setpgroup;
  ButtonGroup  pickcalcgroup;	
  JRadioButton classicalradio;
  JRadioButton qmmmradio;
  JRadioButton quantumradio;
  JRadioButton set1button;
  JRadioButton set2button;	
  JRadioButton set3button;
  JCheckBox setp1button;
  JCheckBox setp2button;
  JRadioButton soluteradio;
  JRadioButton solventradio;
  JRadioButton averageradio;
  JRadioButton applyradio;
  JRadioButton noneradio;
  JRadioButton allradio;
  JRadioButton solventrbutton;
  JRadioButton soluterbutton;
  JRadioButton nonHradio;
  JRadioButton topradio;
  JRadioButton nonbradio;
  JRadioButton solvradio;
  JRadioButton soluradio;
  JRadioButton resetradio;
  JRadioButton lnoneradio;
  JRadioButton pairsradio;
  JRadioButton boxsizeradio;
  JRadioButton verifyradio;
  JRadioButton extraradio;
  JRadioButton energyradio;
  
  IntegerField maxfield;
  DoubleField tolerfield;
  IntegerField grid1field;
  IntegerField grid2field;
  IntegerField grid3field;
  DoubleField alphafield;
  IntegerField orderfield;
  IntegerField nodesfield;
  IntegerField fftfield;
  DoubleField shortfield;
  DoubleField longfield;
  IntegerField freqfield;
  IntegerField maxit1field;
  IntegerField maxit2field;
  DoubleField tol1field;
  DoubleField tol2field;
  DoubleField scale1field;
  DoubleField scale2field;
  IntegerField stepfield;
  IntegerField statfield;

  private nwchem_Task parentFrame = null;
    
////////////////////////////////////////////////////////////////////////
//	Building constraints for panel with settings common between
//	all calculations
///////////////////////////////////////////////////////////////////////

  public nwchem_MD_common(nwchem_Task task) {
    
    defaultfont = new Font("Dialog", Font.BOLD,12);
    disabledfont = new Font("Dialog", Font.PLAIN,12);
    parentFrame = task;
    
    GridBagLayout commonlayout = new GridBagLayout();
    GridBagConstraints commonconstraints = new GridBagConstraints();
    setLayout(commonlayout);
    setgroup = new ButtonGroup();
    
    interactlabel = new JLabel("Interaction:");
    interactlabel.setFont(defaultfont);
    interactlabel.setForeground(Color.black);
    addComponent(this,interactlabel,0,0,4,1,20,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    String[] interactstring = {"effective pair","first order polarization",
			       "SCF polarization", "particle mesh Ewald"};
    interactbox = new JComboBox(interactstring);
    addComponent(this,interactbox,4,0,2,1,20,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    interactbox.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	parentFrame.interaction=interactbox.getSelectedIndex();
	enableSelections(parentFrame);
      }});
    
    JPanel setbuttons = new JPanel();
    setbuttons.setLayout(new GridBagLayout());
    addComponent(this,setbuttons,6,0,3,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);

    setlabel = new JLabel(" set: ");
    setlabel.setFont(defaultfont);
    setlabel.setForeground(Color.black);
    addComponent(setbuttons,setlabel,1,0,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    set1button = new JRadioButton("1");
    set2button = new JRadioButton("2");
    set3button = new JRadioButton("3");
    setgroup.add(set1button);
    setgroup.add(set2button);
    setgroup.add(set3button);
    addComponent(setbuttons,set1button,2,0,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    addComponent(setbuttons,set2button,3,0,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    addComponent(setbuttons,set3button,4,0,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    set1button.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	parentFrame.set=1;
	enableSelections(parentFrame);
      }});
    set2button.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	parentFrame.set=2;
	enableSelections(parentFrame);
      }});
    set3button.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	parentFrame.set=3;
	enableSelections(parentFrame);
      }});
    
    setp1button = new JCheckBox("P1");
    setp2button = new JCheckBox("P2");
    addComponent(setbuttons,setp1button,5,0,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    addComponent(setbuttons,setp2button,6,0,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    setp1button.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	parentFrame.pset2=setp1button.isSelected();
      }});
    setp2button.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	parentFrame.pset3=setp2button.isSelected();
      }});
    
    maxlabel = new JLabel("max iters:");
    maxlabel.setFont(defaultfont);
    maxlabel.setForeground(Color.gray);
    addComponent(this,maxlabel,4,1,4,1,20,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    maxfield = new IntegerField(parentFrame.polmax,14,5);
    addComponent(this,maxfield,5,1,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    //    maxfield.setEnabled(false);
    maxfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ parentFrame.polmax=maxfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    tolerlabel = new JLabel("tolerance:");
    tolerlabel.setFont(defaultfont);
    tolerlabel.setForeground(Color.gray);
    addComponent(this,tolerlabel,6,1,1,1,20,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    tolerfield = new DoubleField(parentFrame.poltol,0.0010,10);
    addComponent(this,tolerfield,7,1,1,1,20,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    tolerfield.setEnabled(false);
    tolerfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ parentFrame.poltol=tolerfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    gridlabel = new JLabel("grid");
    gridlabel.setFont(defaultfont);
    gridlabel.setForeground(Color.gray);
    addComponent(this,gridlabel,4,2,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    grid1field = new IntegerField(parentFrame.gridx,8,5);
    grid2field = new IntegerField(parentFrame.gridy,8,5);
    grid3field = new IntegerField(parentFrame.gridz,8,5);
    addComponent(this,grid1field,5,2,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    addComponent(this,grid2field,6,2,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    addComponent(this,grid3field,7,2,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    grid1field.setEnabled(false);
    grid2field.setEnabled(false);
    grid3field.setEnabled(false);
    grid1field.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ parentFrame.gridx=grid1field.getValue();};
      public void focusGained(FocusEvent e){}
    });
    grid2field.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ parentFrame.gridy=grid2field.getValue();};
      public void focusGained(FocusEvent e){}
    });
    grid3field.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ parentFrame.gridz=grid3field.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    alphalabel = new JLabel("alpha:");
    alphalabel.setFont(defaultfont);
    alphalabel.setForeground(Color.gray);
    addComponent(this,alphalabel,4,3,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    alphafield = new DoubleField(parentFrame.alpha,0.0,10);
    addComponent(this,alphafield,5,3,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    alphafield.setEnabled(false);
    alphafield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ parentFrame.alpha=alphafield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    orderlabel = new JLabel("order:");
    orderlabel.setFont(defaultfont);
    orderlabel.setForeground(Color.gray);
    addComponent(this,orderlabel,6,3,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    orderfield = new IntegerField(parentFrame.order,4,5);
    addComponent(this,orderfield,7,3,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    orderfield.setEnabled(false);
    orderfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ parentFrame.order=orderfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    nodeslabel = new JLabel("nodes:");
    nodeslabel.setFont(defaultfont);
    nodeslabel.setForeground(Color.gray);
    addComponent(this,nodeslabel,8,2,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    nodesfield = new IntegerField(parentFrame.nodpme,0,5);
    addComponent(this,nodesfield,9,2,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    nodesfield.setEnabled(false);
    nodesfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ parentFrame.nodpme=nodesfield.getValue();};
      public void focusGained(FocusEvent e){}
    });

    fftlabel = new JLabel("fft:");
    fftlabel.setFont(defaultfont);
    fftlabel.setForeground(Color.gray);
    addComponent(this,fftlabel,8,3,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    fftfield = new IntegerField(parentFrame.fft,1,5);
    addComponent(this,fftfield,9,3,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    fftfield.setEnabled(false);
    fftfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ parentFrame.fft=fftfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    cutlabel = new JLabel("Cut off:");
    cutlabel.setFont(defaultfont);
    cutlabel.setForeground(Color.black);
    addComponent(this,cutlabel,0,3,4,1,20,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    String[] cutstrings = {"single range", "twin range"};
    cutbox = new JComboBox(cutstrings);
    addComponent(this,cutbox,4,4,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    cutbox.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	parentFrame.twin=cutbox.getSelectedIndex()==1;
	enableSelections(parentFrame);
      }});
    
    shortlabel = new JLabel("short:");
    shortlabel.setFont(defaultfont);
    shortlabel.setForeground(Color.gray);
    addComponent(this,shortlabel,5,4,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    shortfield = new DoubleField(parentFrame.r_short,0.9,5);
    addComponent(this,shortfield,6,4,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    shortfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ parentFrame.r_short=shortfield.getValue();};
      public void focusGained(FocusEvent e){}
    });

    longlabel = new JLabel("long:");
    longlabel.setFont(defaultfont);
    longlabel.setForeground(Color.gray);
    addComponent(this,longlabel,5,5,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    longfield = new DoubleField(parentFrame.r_long,1.2,5);
    addComponent(this,longfield,6,5,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    longfield.setEnabled(false);
    longfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ parentFrame.r_long=longfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    freqlabel = new JLabel("frequency:  ");
    freqlabel.setFont(defaultfont);
    freqlabel.setForeground(Color.gray);
    addComponent(this,freqlabel,7,4,2,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    freqfield = new IntegerField(parentFrame.freq_long,5,5);
    addComponent(this,freqfield,8,4,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    freqfield.setEnabled(false);
    freqfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ parentFrame.freq_long=freqfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    shakelabel = new JLabel("Shake:");
    shakelabel.setFont(defaultfont);
    shakelabel.setForeground(Color.black);
    addComponent(this,shakelabel,0,6,4,1,20,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    solventradio = new JRadioButton("solvent");
    addComponent(this,solventradio,4,6,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    solventradio.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	parentFrame.shake_w=solventradio.isSelected();
	enableSelections(parentFrame);
      }});
    
    maxit1label = new JLabel("max iters:");
    maxit1label.setFont(defaultfont);
    maxit1label.setForeground(Color.gray);
    addComponent(this,maxit1label,5,6,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    maxit1field = new IntegerField(parentFrame.shake_w_iter,100, 5);
    addComponent(this,maxit1field,6,6,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);	
    maxit1field.setEnabled(false);
    maxit1field.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ parentFrame.shake_w_iter=maxit1field.getValue();};
      public void focusGained(FocusEvent e){}
    });

    tol1label = new JLabel("tolerance:");
    tol1label.setFont(defaultfont);
    tol1label.setForeground(Color.gray);
    addComponent(this,tol1label,7,6,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);	
    
    tol1field = new DoubleField(parentFrame.shake_w_toler,0.0010, 5);
    addComponent(this,tol1field,8,6,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    tol1field.setEnabled(false);
    tol1field.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ parentFrame.shake_w_toler=tol1field.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    soluteradio = new JRadioButton("solute");
    addComponent(this,soluteradio,4,7,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    soluteradio.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	parentFrame.shake_s=soluteradio.isSelected();
	enableSelections(parentFrame);
      }});
    
    maxit2label = new JLabel("max iters:");
    maxit2label.setFont(defaultfont);
    maxit2label.setForeground(Color.gray);
    addComponent(this,maxit2label,5,7,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    maxit2field = new IntegerField(parentFrame.shake_s_iter,100,5);	
    addComponent(this,maxit2field,6,7,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    maxit2field.setEnabled(false);
    maxit2field.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ parentFrame.shake_s_iter=maxit2field.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    tol2label = new JLabel("tolerance:");
    tol2label.setFont(defaultfont);
    tol2label.setForeground(Color.gray);	
    addComponent(this,tol2label,7,7,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    tol2field = new DoubleField(parentFrame.shake_s_toler,0.0010,5);
    addComponent(this,tol2field,8,7,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    tol2field.setEnabled(false);
    tol2field.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){parentFrame.shake_s_toler=tol2field.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    distlabel = new JLabel("dist. rest.:");
    distlabel.setFont(defaultfont);
    distlabel.setForeground(Color.black);
    addComponent(this,distlabel,0,8,4,1,20,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    applyradio = new JRadioButton("apply");
    addComponent(this,applyradio,4,8,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    applyradio.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	parentFrame.distar=applyradio.isSelected();
	enableSelections(parentFrame);
      }});
    
    averageradio = new JRadioButton("average");
    addComponent(this,averageradio,5,8,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    averageradio.setEnabled(false);
    
    scale1label = new JLabel("scale:");
    scale1label.setFont(defaultfont);
    scale1label.setForeground(Color.gray);
    addComponent(this,scale1label,6,8,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    scale1field = new DoubleField(parentFrame.drsscl,1.0,5);
    addComponent(this,scale1field,7,8,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    scale1field.setEnabled(false);
    scale1field.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){parentFrame.drsscl=scale1field.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    scale2label = new JLabel("scale:");
    scale2label.setFont(defaultfont);
    scale2label.setForeground(Color.gray);
    addComponent(this,scale2label,8,8,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    scale2field = new DoubleField(parentFrame.drsscl,1.0,5);
    addComponent(this,scale2field,9,8,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    scale2field.setEnabled(false);
    scale2field.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){parentFrame.drsscl=scale2field.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    fixlabel = new JLabel("fix:");
    fixlabel.setFont(defaultfont);
    fixlabel.setForeground(Color.black);
    addComponent(this,fixlabel,0,9,4,1,20,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    String[] fixstrings = {"fix", "unfix"};
    fixbox = new JComboBox(fixstrings);
    addComponent(this,fixbox,4,9,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    noneradio = new JRadioButton("none");
    addComponent(this,noneradio,5,9,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    allradio = new JRadioButton("all");
    addComponent(this,allradio,6,9,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    solventrbutton = new JRadioButton("solvent");
    addComponent(this,solventrbutton,7,9,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    soluterbutton = new JRadioButton("solute");
    addComponent(this,soluterbutton,8,9,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    nonHradio = new JRadioButton("non hydrogen");
    addComponent(this,nonHradio,9,9,2,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    printlabel = new JLabel("Print:");
    printlabel.setFont(defaultfont);
    printlabel.setForeground(Color.black);
    addComponent(this,printlabel,0,10,4,1,20,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    topradio = new JRadioButton("topology");
    addComponent(this,topradio,4,10,2,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    topradio.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	parentFrame.prt_top=topradio.isSelected();
	enableSelections(parentFrame);
      }});
    
    nonbradio = new JRadioButton("nonb");
    addComponent(this,nonbradio,6,10,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    nonbradio.setEnabled(false);
    
    solvradio = new JRadioButton("solvent");
    addComponent(this,solvradio,7,10,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    solvradio.setEnabled(false);
    
    soluradio = new JRadioButton("solute");
    addComponent(this,soluradio,8,10,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    soluradio.setEnabled(false);
    
    steplabel = new JLabel("step:");
    steplabel.setFont(defaultfont);
    steplabel.setForeground(Color.black);
    addComponent(this,steplabel,4,11,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);

    stepfield = new IntegerField(parentFrame.prt_step,0,5);
    addComponent(this,stepfield,5,11,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    stepfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){parentFrame.prt_step=stepfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    extraradio = new JRadioButton("extra");
    addComponent(this,extraradio,6,11,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    energyradio = new JRadioButton("energy");
    addComponent(this,energyradio,7,11,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    statlabel = new JLabel("stat:");
    statlabel.setFont(defaultfont);
    statlabel.setForeground(Color.black);
    addComponent(this,statlabel,4,12,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    statfield = new IntegerField(parentFrame.prt_stat,0,5);
    addComponent(this,statfield,5,12,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    statfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){parentFrame.prt_stat=statfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    fillerlabel = new JLabel("                 ");
    addComponent(this,fillerlabel,10,18,3,3,35,35,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    setVisible(true);

	enableSelections(parentFrame);
  }
 
  void enableSelections(nwchem_Task tsk){
     
    interactbox.setSelectedIndex(tsk.interaction);
    if(tsk.set == 1){
      set1button.setSelected(true);
      setp1button.setEnabled(true);
      setp2button.setEnabled(true);
    } else {
      if(tsk.set == 2){ set2button.setSelected(true); };
      if(tsk.set == 3){ set3button.setSelected(true); };
      setp1button.setEnabled(false);
      setp2button.setEnabled(false);
    };
    setp1button.setSelected(tsk.pset2);
    setp2button.setSelected(tsk.pset3);
    if(tsk.interaction == 0 || tsk.interaction == 1){
      maxlabel.setForeground(Color.gray);
      maxfield.setEnabled(false);
      tolerlabel.setForeground(Color.gray);
      tolerfield.setEnabled(false);
      gridlabel.setForeground(Color.gray);
      grid1field.setEnabled(false);
      grid2field.setEnabled(false);
      grid3field.setEnabled(false);	
      alphalabel.setForeground(Color.gray);
      alphafield.setEnabled(false);
      orderlabel.setForeground(Color.gray);
      orderlabel.setEnabled(false);
      nodeslabel.setForeground(Color.gray);
      nodesfield.setEnabled(false);
      fftlabel.setForeground(Color.gray);
      fftfield.setEnabled(false);
    };
    if(tsk.interaction == 2){
      maxlabel.setForeground(Color.black);
      maxfield.setEnabled(true);
      tolerlabel.setForeground(Color.black);
      tolerfield.setEnabled(true);
      gridlabel.setForeground(Color.gray);
      grid1field.setEnabled(false);
      grid2field.setEnabled(false);
      grid3field.setEnabled(false);	
      alphalabel.setForeground(Color.gray);
      alphafield.setEnabled(false);
      orderlabel.setForeground(Color.gray);
      orderlabel.setEnabled(false);
      nodeslabel.setForeground(Color.gray);
      nodesfield.setEnabled(false);
      fftlabel.setForeground(Color.gray);
      fftfield.setEnabled(false);
    };
    if(tsk.interaction == 3){
      maxlabel.setForeground(Color.gray);
      maxfield.setEnabled(false);
      tolerlabel.setForeground(Color.gray); 
      tolerfield.setEnabled(false);
      gridlabel.setForeground(Color.black);
      grid1field.setEnabled(true); 
      grid2field.setEnabled(true); 
      grid3field.setEnabled(true); 	
      alphalabel.setForeground(Color.black);
      alphafield.setEnabled(true);
      orderlabel.setForeground(Color.black);
      orderlabel.setEnabled(true); 
      nodeslabel.setForeground(Color.black);
      nodesfield.setEnabled(true); 
      fftlabel.setForeground(Color.black);
      fftfield.setEnabled(true); 
    };
    if(tsk.prt_top){
      nonbradio.setEnabled(true); 
      solvradio.setEnabled(true); 
      soluradio.setEnabled(true);
    } else {
      nonbradio.setSelected(false);
      solvradio.setSelected(false);
      soluradio.setSelected(false);
      nonbradio.setEnabled(false);
      solvradio.setEnabled(false);
      soluradio.setEnabled(false);
    };
    if(tsk.distar){
      averageradio.setEnabled(true); 
      scale1label.setForeground(Color.black);
      scale1field.setEnabled(true);
      scale2label.setForeground(Color.black);
      scale2field.setEnabled(true);
    } else {
      averageradio.setEnabled(false);
      scale1label.setForeground(Color.gray);
      scale1field.setEnabled(false);
      scale2label.setForeground(Color.gray);
      scale2field.setEnabled(false);
    };
    solventradio.setSelected(tsk.shake_w);
    if(tsk.shake_w){
      maxit1label.setForeground(Color.black);
      maxit1field.setEnabled(true);  
      tol1label.setForeground(Color.black);
      tol1field.setEnabled(true);
    } else {
      maxit1label.setForeground(Color.gray);
      maxit1field.setEnabled(false);
      tol1label.setForeground(Color.gray);
      tol1field.setEnabled(false);
    };
    soluteradio.setSelected(tsk.shake_s);
    if(tsk.shake_s){
      maxit2label.setForeground(Color.black);
      maxit2field.setEnabled(true);  
      tol2label.setForeground(Color.black);
      tol2field.setEnabled(true);
    } else {
      maxit2label.setForeground(Color.gray);
      maxit2field.setEnabled(false);
      tol2label.setForeground(Color.gray);
      tol2field.setEnabled(false);
    };
    if(tsk.twin){
      cutbox.setSelectedIndex(1);
      shortlabel.setForeground(Color.black);
      longlabel.setForeground(Color.black);
      longfield.setEnabled(true); 
      freqlabel.setForeground(Color.black);
      freqfield.setEnabled(true);
    } else {
      cutbox.setSelectedIndex(0);
      shortlabel.setForeground(Color.gray);
      longlabel.setForeground(Color.gray);
      longfield.setEnabled(false);
      freqlabel.setForeground(Color.gray);
      freqfield.setEnabled(false);
    };
  }

///////////////////////////////////////////////////////////////////////////
//	Helper method for setting GridBagLayout
//////////////////////////////////////////////////////////////////////////
   
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
  
  public void setInitialVariables() {

    grid1field.setText(new Integer(parentFrame.gridx).toString());
    grid2field.setText(new Integer(parentFrame.gridy).toString());
    grid3field.setText(new Integer(parentFrame.gridz).toString());
    alphafield.setText(new Double(parentFrame.alpha).toString());
    orderfield.setText(new Integer(parentFrame.order).toString());
    nodesfield.setText(new Integer(parentFrame.nodpme).toString());
    fftfield.setText(new Integer(parentFrame.fft).toString());
    if (parentFrame.twin) cutbox.setSelectedIndex(0);
    else cutbox.setSelectedIndex(1);
    shortfield.setText(new Double(parentFrame.r_short).toString());
    longfield.setText(new Double(parentFrame.r_long).toString());
    if (parentFrame.shake_w) solventradio.setSelected(true);
    if (parentFrame.shake_s) soluteradio.setSelected(true);
    maxit1field.setText(new Integer(parentFrame.shake_w_iter).toString());
    maxit2field.setText(new Integer(parentFrame.shake_s_iter).toString());
    tol1field.setText(new Double(parentFrame.shake_w_toler).toString());
    tol2field.setText(new Double(parentFrame.shake_s_toler).toString());
    if (parentFrame.distar) applyradio.setSelected(true);
    if (parentFrame.draver) averageradio.setSelected(true);
    scale1field.setText(new Double(parentFrame.drsscl).toString());
    scale2field.setText(new Double(parentFrame.nfdrss).toString());
    if (parentFrame.fix) fixbox.setSelectedIndex(0);
    if (parentFrame.unfix) fixbox.setSelectedIndex(1);
    if (parentFrame.fix_none) noneradio.setSelected(true);
    if (parentFrame.fix_all) allradio.setSelected(true);
    if (parentFrame.fix_w) solventrbutton.setSelected(true);
    if (parentFrame.fix_s) soluterbutton.setSelected(true);
    if (parentFrame.fix_X) nonHradio.setSelected(true);
    if (parentFrame.prt_top) topradio.setSelected(true);
    if (parentFrame.prt_topn) nonbradio.setSelected(true);
    if (parentFrame.prt_topw) solvradio.setSelected(true);
    if (parentFrame.prt_tops) soluradio.setSelected(true);
    stepfield.setText(new Integer(parentFrame.prt_step).toString());
    statfield.setText(new Integer(parentFrame.prt_stat).toString());
    if (parentFrame.prt_extra) extraradio.setSelected(true);
    if (parentFrame.prt_energy) energyradio.setSelected(true);
    expectfield.setText(new Integer(parentFrame.prt_expect).toString());
    if (parentFrame.load_pair) pairsradio.setSelected(true);
    maxupfield.setText(new Integer(parentFrame.load_num_pair).toString());
    if (parentFrame.load_reset) resetradio.setSelected(true);
    if (parentFrame.load_size) boxsizeradio.setSelected(true);
    
  }
  
}
