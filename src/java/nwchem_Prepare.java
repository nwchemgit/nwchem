import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

class nwchem_Prepare extends JPanel implements ActionListener, ChangeListener {

  Font defaultFont;

  JTabbedPane inputpane;
  TextField systemName;
  nwchem_Task tsk;

  public nwchem_Prepare(nwchem_Task task){

    setLayout(new GridBagLayout());
    setForeground(Color.black);
    setBackground(Color.lightGray);
    setFont(defaultFont);

    tsk=task;

    System.out.println("nwchem_Prepare constructor");

    JLabel systemLabel = new JLabel("system: ");
    addComponent(this,systemLabel,2,0,1,1,1,1,
    		 GridBagConstraints.EAST,GridBagConstraints.EAST);
    systemLabel.setForeground(Color.black);
    
    System.out.println(tsk.system);
    systemName = new TextField(tsk.system,25);
    addComponent(this,systemName,3,0,3,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.WEST);

    inputpane = new JTabbedPane();
    inputpane.setForeground(Color.black);
    inputpane.setBackground(Color.gray);
    inputpane.setFont(defaultFont);
    inputpane.addChangeListener(this);
    addComponent(this,inputpane,0,1,6,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.CENTER);

    nwchem_Prepare_topology topologypanel = new nwchem_Prepare_topology(task);
    inputpane.insertTab("Topology",null,topologypanel,"Topology",0);

    nwchem_Prepare_restart restartpanel = new nwchem_Prepare_restart(task);
    inputpane.insertTab("Restart",null,restartpanel,"Restart",1);

    validate();
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

  public void mouseClicked(MouseEvent mouse) {}
  
}
