import java.applet.Applet;
import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.lang.*;

class nwchem_MD_control extends JPanel {
  
  ///////////////////////////////////////////////////////////////////////
  //	Variables that all the panels share.
  //////////////////////////////////////////////////////////////////////
  Font defaultfont;
  Font disabledfont;
  JLabel fixlabel;
  JLabel printlabel;
  JLabel steplabel;
  JLabel statlabel;
  JLabel expectlabel;
  JLabel controllabel;
  JLabel loadlabel;
  JLabel cnodeslabel;
  JLabel mwmlabel;
  JLabel testlabel;
  JLabel boxeslabel;
  JLabel extralabel;
  JLabel msalabel;
  JLabel debuglabel;
  JLabel memorylabel;
  JLabel stacklabel;
  JLabel heaplabel;
  JLabel globallabel;
  JLabel limitlabel;
  JLabel maxuplabel;
  JLabel scalelabel;
  JLabel calculationlabel;
  JLabel fillerlabel;
  
  JRadioButton resetradio;
  JRadioButton lnoneradio;
  JRadioButton pairsradio;
  JRadioButton boxsizeradio;
  JRadioButton verifyradio;
  JRadioButton extraradio;
  JRadioButton energyradio;
  
  IntegerField expectfield;
  IntegerField nodes1field;
  IntegerField nodes2field;
  IntegerField nodes3field;
  IntegerField mwmfield;
  IntegerField testfield;
  IntegerField boxes1field;
  IntegerField boxes2field;
  IntegerField boxes3field;
  IntegerField extrafield;
  IntegerField msafield;
  IntegerField debugfield;
  IntegerField limitfield;
  IntegerField maxupfield;
  DoubleField scalefield;

  private nwchem_Task task = null;
    
////////////////////////////////////////////////////////////////////////
//	Building constraints for panel with settings common between
//	all calculations
///////////////////////////////////////////////////////////////////////

  public nwchem_MD_control(nwchem_Task t) {
    
    defaultfont = new Font("Dialog", Font.BOLD,12);
    disabledfont = new Font("Dialog", Font.PLAIN,12);
    task = t;
    
    GridBagLayout commonlayout = new GridBagLayout();
    GridBagConstraints commonconstraints = new GridBagConstraints();
    setLayout(commonlayout);
    
    printlabel = new JLabel("Print:");
    printlabel.setFont(defaultfont);
    printlabel.setForeground(Color.black);
    addComponent(this,printlabel,0,0,4,1,20,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    expectlabel = new JLabel("expect:");
    expectlabel.setFont(defaultfont);
    expectlabel.setForeground(Color.black);
    addComponent(this,expectlabel,4,0,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    expectfield = new IntegerField(task.prt_expect,1,5);
    addComponent(this,expectfield,5,0,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    expectfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){task.prt_expect=expectfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    loadlabel = new JLabel("Load balance:");
    loadlabel.setFont(defaultfont);
    loadlabel.setForeground(Color.black);
    addComponent(this,loadlabel,0,1,2,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    lnoneradio = new JRadioButton("none");
    lnoneradio.setSelected(task.load_none);
    addComponent(this,lnoneradio,4,1,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    lnoneradio.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	task.load_none=lnoneradio.isSelected();
	enableSelections(task);
      }});
    
    resetradio = new JRadioButton("reset");
    resetradio.setSelected(task.load_reset);
    addComponent(this,resetradio,4,2,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    resetradio.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	task.load_reset=resetradio.isSelected();
	enableSelections(task);
      }});
    
    pairsradio = new JRadioButton("pairs:");
    pairsradio.setSelected(task.load_pair);
    addComponent(this,pairsradio,5,1,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    pairsradio.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	task.load_pair=pairsradio.isSelected();
	enableSelections(task);
      }});
    
    boxsizeradio = new JRadioButton("boxsize");
    boxsizeradio.setSelected(task.load_size);
    addComponent(this,boxsizeradio,5,2,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    boxsizeradio.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	task.load_size=boxsizeradio.isSelected();
	enableSelections(task);
      }});
    
    scalelabel = new JLabel("scale:");
    scalelabel.setFont(defaultfont);
    scalelabel.setForeground(Color.gray);
    addComponent(this,scalelabel,6,2,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    scalefield = new DoubleField(task.load_factor,0.75,5);
    addComponent(this,scalefield,7,2,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    scalefield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){task.load_factor=scalefield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    maxuplabel = new JLabel("maxup:");
    maxuplabel.setFont(defaultfont);
    maxuplabel.setForeground(Color.gray);
    addComponent(this,maxuplabel,6,1,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    maxupfield = new IntegerField(task.load_num_pair,5,5);
    addComponent(this,maxupfield,7,1,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    maxupfield.setEnabled(false);
    maxupfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){task.load_num_pair=maxupfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    

    controllabel = new JLabel("Decomposition:");
    controllabel.setFont(defaultfont);
    controllabel.setForeground(Color.black);
    addComponent(this,controllabel,0,3,4,1,20,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);

    cnodeslabel = new JLabel("nodes:");
    cnodeslabel.setFont(defaultfont);
    cnodeslabel.setForeground(Color.black);
    addComponent(this,cnodeslabel,4,3,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    nodes1field = new IntegerField(task.npx,0,5);
    nodes2field = new IntegerField(task.npy,0,5);
    nodes3field = new IntegerField(task.npz,0,5);
    addComponent(this,nodes1field,5,3,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    addComponent(this,nodes2field,6,3,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    addComponent(this,nodes3field,7,3,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    nodes1field.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){task.npx=nodes1field.getValue();};
      public void focusGained(FocusEvent e){}
    });
    nodes2field.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){task.npy=nodes2field.getValue();};
      public void focusGained(FocusEvent e){}
    });
    nodes3field.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){task.npz=nodes3field.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    mwmlabel = new JLabel("mwm:");
    mwmlabel.setFont(defaultfont);
    mwmlabel.setForeground(Color.black);
    addComponent(this,mwmlabel,4,5,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    mwmfield = new IntegerField(task.mwm,0,5);
    addComponent(this,mwmfield,5,5,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    mwmfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){task.mwm=mwmfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    testlabel = new JLabel("test:");
    testlabel.setFont(defaultfont);
    testlabel.setForeground(Color.red);
    addComponent(this,testlabel,8,3,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    testfield = new IntegerField(task.test,0,5);
    addComponent(this,testfield,9,3,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    testfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){task.test=testfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    boxeslabel = new JLabel("boxes");
    boxeslabel.setFont(defaultfont);
    boxeslabel.setForeground(Color.black);
    addComponent(this,boxeslabel,4,4,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    boxes1field = new IntegerField(task.nbx,0,5);
    boxes2field = new IntegerField(task.nby,0,5);
    boxes3field = new IntegerField(task.nbz,0,5);
    addComponent(this,boxes1field,5,4,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    addComponent(this,boxes2field,6,4,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    addComponent(this,boxes3field,7,4,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    boxes1field.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){task.nbx=boxes1field.getValue();};
	  public void focusGained(FocusEvent e){}
    });
    boxes2field.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){task.nby=boxes2field.getValue();};
      public void focusGained(FocusEvent e){}
    });
    boxes3field.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){task.nbz=boxes3field.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    msalabel = new JLabel("msa:");
    msalabel.setFont(defaultfont);
    msalabel.setForeground(Color.black);
    addComponent(this,msalabel,6,5,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    msafield = new IntegerField(task.msa,0,5);
    addComponent(this,msafield,7,5,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    msafield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){task.msa=msafield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    extralabel = new JLabel("extra");
    extralabel.setFont(defaultfont);
    extralabel.setForeground(Color.red);
    addComponent(this,extralabel,8,4,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    extrafield = new IntegerField(task.mad,6,5);
    addComponent(this,extrafield,9,4,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    extrafield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){task.mad=extrafield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    debuglabel = new JLabel("debug:");
    debuglabel.setFont(defaultfont);
    debuglabel.setForeground(Color.red);
    addComponent(this,debuglabel,8,5,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    debugfield = new IntegerField(task.debug,0,5);
    addComponent(this,debugfield,9,5,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    debugfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){task.debug=debugfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    memorylabel = new JLabel("Memory:");
    memorylabel.setFont(defaultfont);
    memorylabel.setForeground(Color.black);
    addComponent(this,memorylabel,0,6,4,1,20,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    limitlabel = new JLabel("limit:");
    limitlabel.setFont(defaultfont);
    limitlabel.setForeground(Color.black);
    addComponent(this,limitlabel,4,6,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    limitfield = new IntegerField(task.Limit,0,5);
    addComponent(this,limitfield,5,6,1,1,5,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    limitfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){task.Limit=limitfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    fillerlabel = new JLabel("                 ");
    addComponent(this,fillerlabel,10,7,3,3,35,35,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    setVisible(true);

    enableSelections(task);
  }
 
  void enableSelections(nwchem_Task tsk){
     
    if(tsk.load_none){
      pairsradio.setEnabled(false);
      boxsizeradio.setEnabled(false);
      scalelabel.setForeground(Color.gray);
      scalefield.setEnabled(false);
      maxuplabel.setForeground(Color.gray);
      maxupfield.setEnabled(false);
    } else {
      pairsradio.setEnabled(true);
      boxsizeradio.setEnabled(true);
      if(tsk.load_size){
	scalelabel.setForeground(Color.black);
	scalefield.setEnabled(true);
	if(tsk.load_pair){
	  maxuplabel.setForeground(Color.black);
	  maxupfield.setEnabled(true);
	};
      } else {
	scalelabel.setForeground(Color.gray);
	scalefield.setEnabled(false);
	maxuplabel.setForeground(Color.gray);
	maxupfield.setEnabled(false);
      };
      if(tsk.load_pair){
	if(tsk.load_size){
	  maxuplabel.setForeground(Color.black);
	  maxupfield.setEnabled(true);
	};
      } else {
	maxuplabel.setForeground(Color.gray);
	maxupfield.setEnabled(false);
      };
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
  
}
