import java.applet.Applet;
import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

class nwchem_MD_optimize extends JPanel implements ActionListener {
    
  GridBagLayout optimizelayout;
  GridBagConstraints optimizeconstraints;
  Font defaultfont;
  
  JLabel steeplabel;
  JLabel steepmaxlabel;
  JLabel conjlabel;
  JLabel conjmaxlabel;
  JLabel initlabel;
  JLabel tolinitlabel;
  JLabel msteplabel;
  JLabel intlabel;
  JLabel tolintlabel;
  JLabel refrlabel;
  JLabel updateslabel;
  JLabel pairslabel;
  JLabel recordlabel;
  JLabel minxlabel;
  JLabel proplabel;
  JLabel fillerlabel;
  
  IntegerField sditer;
  DoubleField sdinit;
  DoubleField sdxmin;
  DoubleField sdxmax;
  IntegerField cgiter;
  DoubleField cginit;
  DoubleField cgxmin;
  IntegerField cgcycle;
  IntegerField pairsfield;
  IntegerField minxfield;
  IntegerField propfield;
  private nwchem_Task task = null;
  

public nwchem_MD_optimize(nwchem_Task t) {

    optimizelayout = new GridBagLayout();
    setLayout(optimizelayout);
    optimizeconstraints = new GridBagConstraints();
    
    defaultfont = new Font("Dialog", Font.BOLD,12);
    task = t;

    steeplabel = new JLabel("Steepest descent");
    steeplabel.setFont(defaultfont);
    steeplabel.setForeground(Color.black);
    addComponent(this,steeplabel,0,0,2,1,20,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    steepmaxlabel = new JLabel("max steps:");
    steepmaxlabel.setFont(defaultfont);
    steepmaxlabel.setForeground(Color.black);
    addComponent(this,steepmaxlabel,2,0,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    sditer = new IntegerField(task.sd_iter,100,5);
    addComponent(this,sditer,3,0,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    sditer.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.sd_iter=sditer.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    initlabel = new JLabel("initial step");
    initlabel.setFont(defaultfont);
    initlabel.setForeground(Color.black);
    addComponent(this,initlabel,4,0,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    sdinit = new DoubleField(task.sd_init,0.01, 5);
    addComponent(this,sdinit,5,0,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    sdinit.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.sd_init=sdinit.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    tolinitlabel = new JLabel("tolerance:");
    tolinitlabel.setForeground(Color.black);
    tolinitlabel.setFont(defaultfont);
    addComponent(this,tolinitlabel,6,0,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    sdxmin = new DoubleField(task.sd_xmin,1.0E-4,5);
    addComponent(this,sdxmin,7,0,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    sdxmin.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.sd_xmin=sdxmin.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    msteplabel = new JLabel("maximum step:");
    msteplabel.setFont(defaultfont);
    msteplabel.setForeground(Color.black);
    addComponent(this,msteplabel,8,0,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);

    sdxmax = new DoubleField(task.sd_xmax,0.1,5);
    addComponent(this,sdxmax,9,0,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    sdxmax.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.sd_xmax=sdxmax.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    conjlabel = new JLabel("Conjugate gradient");
    conjlabel.setFont(defaultfont);
    conjlabel.setForeground(Color.black);
    addComponent(this,conjlabel,0,1,2,1,20,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);

    conjmaxlabel = new JLabel("max steps:");
    conjmaxlabel.setFont(defaultfont);
    conjmaxlabel.setForeground(Color.black);
    addComponent(this,conjmaxlabel,2,1,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    cgiter = new IntegerField(task.cg_iter,0,5);
    addComponent(this,cgiter,3,1,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    cgiter.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.cg_iter=cgiter.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    intlabel = new JLabel("interval:");
    intlabel.setFont(defaultfont);
    intlabel.setForeground(Color.black);
    addComponent(this,intlabel,4,1,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);

    cginit = new DoubleField(task.cg_init,0.3,5);
    addComponent(this,cginit,5,1,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    cginit.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.cg_init=cginit.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    tolintlabel = new JLabel("tolerance:");
    tolintlabel.setForeground(Color.black);
    tolintlabel.setFont(defaultfont);
    addComponent(this,tolintlabel,6,1,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    cgxmin = new DoubleField(task.cg_xmin,1.0E-4,5);
    addComponent(this,cgxmin,7,1,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    cgxmin.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.cg_xmin=cgxmin.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    refrlabel = new JLabel("refresh cycle:");
    refrlabel.setForeground(Color.black);
    refrlabel.setFont(defaultfont);
    addComponent(this,refrlabel,8,1,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    cgcycle = new IntegerField(task.cg_cycle,10,5);
    addComponent(this,cgcycle,9,1,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    cgcycle.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.cg_cycle=cgcycle.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    updateslabel = new JLabel("Updates:");
    updateslabel.setFont(defaultfont);
    updateslabel.setForeground(Color.black);
    addComponent(this,updateslabel,0,2,2,1,20,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    pairslabel = new JLabel("pairs:");
    pairslabel.setFont(defaultfont);
    pairslabel.setForeground(Color.black);
    addComponent(this,pairslabel,2,2,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    pairsfield = new IntegerField(task.freq_pair,1,5);
    addComponent(this,pairsfield,3,2,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    pairsfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.freq_pair=pairsfield.getValue();};
      public void focusGained(FocusEvent e){}
    });

    recordlabel = new JLabel("Recording:");
    recordlabel.setFont(defaultfont);
    recordlabel.setForeground(Color.black);
    addComponent(this,recordlabel,0,3,2,1,20,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    minxlabel = new JLabel("minim x:");
    minxlabel.setFont(defaultfont);
    minxlabel.setForeground(Color.black);
    addComponent(this,minxlabel,2,3,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    minxfield = new IntegerField(task.freq_mnd,0,5);
    addComponent(this,minxfield,3,3,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    minxfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.freq_mnd=minxfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    proplabel = new JLabel("property:");
    proplabel.setFont(defaultfont);
    proplabel.setForeground(Color.black);
    addComponent(this,proplabel,4,3,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    
    propfield = new IntegerField(task.freq_prp,0,5);
    addComponent(this,propfield,5,3,1,1,10,1,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
    propfield.addFocusListener(new FocusListener(){
      public void focusLost(FocusEvent e){ task.freq_prp=propfield.getValue();};
      public void focusGained(FocusEvent e){}
    });
    
    fillerlabel = new JLabel(" ");
    addComponent(this,fillerlabel,10,4,1,1,35,35,
		 GridBagConstraints.HORIZONTAL,GridBagConstraints.CENTER);
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

}


