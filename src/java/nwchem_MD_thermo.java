import java.applet.Applet;
import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

class nwchem_MD_thermo extends JPanel implements ActionListener {

  JLabel freeElabel;
  JLabel oflabel;
  JLabel overlabel;
  JLabel errorlabel;
  JLabel driftlabel;
  JLabel freelabel;
  JLabel recordlabel;
  JLabel deltalabel;
  JLabel factorlabel;
  JLabel fillerlabel;
  
  DoubleField factorfield;
  IntegerField stepfield;
  IntegerField offield;
  IntegerField overfield; 
  DoubleField errorfield;
  DoubleField driftfield;
  DoubleField deltafield;
  IntegerField freefield; 

  String[] reversestring = {"forward","reverse"};
  String[] extendstring = {"new", "renew", "extend"};

  JComboBox energybox;
  JComboBox reversebox;
  
  JRadioButton sssbutton;
  JRadioButton decomposebutton;
  JRadioButton cnvbutton;
  JRadioButton fetbutton;
  JRadioButton acfbutton; 
  
  Font defaultfont;
  private nwchem_Task task = null;
  
  public nwchem_MD_thermo(nwchem_Task t) {
    
    defaultfont = new Font("Dialog", Font.BOLD,12);
    
    GridBagLayout thermolayout = new GridBagLayout();
    setLayout(thermolayout);
    GridBagConstraints thermoconstraints = new GridBagConstraints();
    
    task = t;

    freeElabel = new JLabel("Free energy:");
    freeElabel.setFont(defaultfont);
    freeElabel.setForeground(Color.black);
    addComponent(this,freeElabel,0,0,2,1,15,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    energybox = new JComboBox(extendstring);
    energybox.setSelectedItem(task.ti_start);
    addComponent(this,energybox,2,0,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    energybox.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	task.ti_direct=extendstring[energybox.getSelectedIndex()];
      }});
    
    reversebox = new JComboBox(reversestring);
    reversebox.setSelectedItem(task.ti_direct);
    addComponent(this,reversebox,3,0,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    reversebox.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	task.ti_direct=reversestring[reversebox.getSelectedIndex()];
      }});
    
    stepfield = new IntegerField(task.ti_steps,21,5);
    addComponent(this,stepfield,4,0,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    stepfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.ti_steps=stepfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    oflabel = new JLabel("of");
    oflabel.setFont(defaultfont);
    oflabel.setForeground(Color.black);
    addComponent(this,oflabel,5,0,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    offield = new IntegerField(task.ti_max,21,5);
    addComponent(this,offield,6,0,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    offield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.ti_max=offield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    overlabel = new JLabel("over");
    overlabel.setFont(defaultfont);
    overlabel.setForeground(Color.black);
    addComponent(this,overlabel,7,0,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    overfield = new IntegerField(task.ti_over,1000,5);
    addComponent(this,overfield,8,0,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    overfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.ti_over=overfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    errorlabel = new JLabel("error:");
    errorlabel.setFont(defaultfont);
    errorlabel.setForeground(Color.black);
    addComponent(this,errorlabel,2,1,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    errorfield = new DoubleField(task.ti_error,5.0,5);
    addComponent(this,errorfield,3,1,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    errorfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.ti_error=errorfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    driftlabel = new JLabel("drift:");
    driftlabel.setFont(defaultfont);
    driftlabel.setForeground(Color.black);
    addComponent(this,driftlabel,2,2,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    driftfield = new DoubleField(task.ti_drift,5.0,5);
    addComponent(this,driftfield,3,2,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    driftfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.ti_drift=driftfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    factorlabel = new JLabel("factor:");
    factorlabel.setFont(defaultfont);
    factorlabel.setForeground(Color.black);
    addComponent(this,factorlabel,2,3,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    factorfield = new DoubleField(task.ti_factor,0.75,5);
    addComponent(this,factorfield,3,3,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    factorfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.ti_factor=factorfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    sssbutton = new JRadioButton("sss");
    sssbutton.setSelected(task.ti_sss);
    addComponent(this,sssbutton,4,1,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    sssbutton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	task.ti_sss=sssbutton.isSelected();
	enableSelections();}
    });
    
    deltalabel = new JLabel("delta:");
    deltalabel.setFont(defaultfont);
    deltalabel.setForeground(Color.gray);
    addComponent(this,deltalabel,5,1,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    deltafield = new DoubleField(task.ti_delta,0.075,5);
    addComponent(this,deltafield,6,1,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    deltafield.setEnabled(false);
    deltafield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.ti_delta=deltafield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    decomposebutton = new JRadioButton("decompose");
    decomposebutton.setSelected(task.ti_decomp);
    addComponent(this,decomposebutton,4,2,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    decomposebutton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	task.ti_decomp=decomposebutton.isSelected(); }
    });
    
    recordlabel = new JLabel("Recording:");
    recordlabel.setFont(defaultfont);
    recordlabel.setForeground(Color.black);
    addComponent(this,recordlabel,0,4,2,1,15,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    freelabel = new JLabel("free:");
    freelabel.setFont(defaultfont);
    freelabel.setForeground(Color.black);
    addComponent(this,freelabel,2,4,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    freefield = new IntegerField(task.freq_gib,0,5);
    addComponent(this,freefield,3,4,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    freefield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.freq_gib=freefield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    cnvbutton = new JRadioButton("cnv");
    cnvbutton.setSelected(task.cnv);
    addComponent(this,cnvbutton,4,4,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    cnvbutton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	task.cnv=cnvbutton.isSelected(); }
    });
    
    fetbutton = new JRadioButton("fet");
    fetbutton.setSelected(task.fet);
    addComponent(this,fetbutton,5,4,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    fetbutton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	task.fet=fetbutton.isSelected(); }
    });
    
    acfbutton = new JRadioButton("acf");
    acfbutton.setSelected(task.acf);
    addComponent(this,acfbutton,6,4,1,1,7,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    acfbutton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	task.acf=acfbutton.isSelected(); }
    });
    
    fillerlabel = new JLabel(" ");
    addComponent(this,fillerlabel,7,5,1,1,7,35,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);

    setVisible(true);
    enableSelections();
  }
  
  void enableSelections(){
    if(task.ti_sss) {
      deltalabel.setForeground(Color.black);
      deltafield.setEnabled(true);
    } else {
      deltalabel.setForeground(Color.gray);
      deltafield.setEnabled(false);
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
  ///////////////////////////////////////////////////////////////////////////
  //	Helper method for setting GridBagLayout
  //////////////////////////////////////////////////////////////////////////
  
  void buildConstraints(GridBagConstraints gbc, int gx, int gy, int gw, int gh, 
			int wx, int wy){
    
    gbc.gridx = gx;
    gbc.gridy = gy;
    gbc.gridwidth = gw;
    gbc.gridheight = gh;
    gbc.weightx = wx;
    gbc.weighty = wy;
  }
  
  public void actionPerformed(ActionEvent e){
    
    ///////////////////////////////////////////////////////////////////////////
    //	The following actions are called such that certain textfields, 
    //	buttons, and comboboxes are enabled or disabled when another 
    //	designated button or textfield has been chosen, in order that only
    //	necessary input is taken.
    //////////////////////////////////////////////////////////////////////////
    
    if (e.getSource() == sssbutton) {
      if (!sssbutton.isSelected()) {
	deltalabel.setForeground(Color.gray);
	deltafield.setEnabled(false);
      }
      
      if (sssbutton.isSelected()) {
	deltalabel.setForeground(Color.black);
	deltafield.setEnabled(true); 
      }
    }
  }
}
