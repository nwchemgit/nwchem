import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

class nwchem_Param extends JFrame implements ActionListener, ChangeListener, WindowListener {
  
  Font defaultFont;

    String ffield = "unknown";
    String file_s = " ";
    String file_x = " ";
    String file_u = " ";
    String file_t = " ";

    int mAtoms=100, mCross=100, mBonds=250, mAngles=250, mTorsions=200, mImpropers=100, mRules=100;

    String Atoms[][];
    String Cross[];
    String Bonds[] = new String[500];
    String Angles[];
    String Torsions[];
    String Impropers[];
    String Rules[][];

    int idAtoms[];
    int idCross[];
    int idBonds[] = new int[500];
    int idAngles[];
    int idTorsions[];
    int idImpropers[];
    int idRules[];

    int nAtoms=0, nCross=0, nBonds=0, nAngles=0, nTorsions=0, nImpropers=0, nRules=0;

  public nwchem_Param(){

    super("Parameters");

    defaultFont = new Font("Dialog", Font.BOLD,12);

    super.getContentPane().setLayout(new GridBagLayout());
    super.getContentPane().setForeground(Color.black);
    super.getContentPane().setBackground(Color.lightGray);
    super.getContentPane().setFont(defaultFont);
    super.addWindowListener(this);

    JPanel header = new JPanel();
    header.setLayout(new GridBagLayout());
    header.setForeground(Color.black);
    header.setBackground(Color.lightGray);
    addComponent(super.getContentPane(),header,0,0,2,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.WEST);

    JButton doneButton = new JButton("Done");
    addComponent(header,doneButton,5,0,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    doneButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	setVisible(false); }});

    setLocation(25,225);	
    setSize(900,700);
    setVisible(true);

    Atoms = new String[mAtoms][2];
    Cross = new String[mCross];
    //    Bonds = new String[mBonds];
    Angles = new String[mAngles];
    Torsions = new String[mTorsions];
    Impropers = new String[mImpropers];
    Rules = new String[mRules][3];

    idAtoms = new int[mAtoms];
    idCross = new int[mCross];
    //    idBonds = new int[mBonds];
    idAngles = new int[mAngles];
    idTorsions = new int[mTorsions];
    idImpropers = new int[mImpropers];
    idRules = new int[mRules];

    try{
      BufferedReader br = new BufferedReader(new FileReader("/home/strtsm/.nwchemrc"));
      String card;
      while((card=br.readLine()) != null){
	  if(card.startsWith("ffield")) {ffield=card.substring(7,12).trim();};
	  if(card.startsWith(ffield+"_s")) {file_s=card.substring(card.indexOf(" ")+1,card.length())+ffield+".par";};
	  if(card.startsWith(ffield+"_x")) {file_x=card.substring(card.indexOf(" ")+1,card.length())+ffield+".par";};
	  if(card.startsWith(ffield+"_u")) {file_u=card.substring(card.indexOf(" ")+1,card.length())+ffield+".par";};
	  if(card.startsWith(ffield+"_t")) {file_t=card.substring(card.indexOf(" ")+1,card.length())+ffield+".par";};
      };
      br.close();
    } catch(Exception e) {e.printStackTrace();};
    if(file_s!=" ") { param_Read(file_s,1);};
    if(file_x!=" ") { param_Read(file_x,2);};
    if(file_u!=" ") { param_Read(file_u,3);};
    if(file_t!=" ") { param_Read(file_t,4);};
    System.out.println(nAtoms+" "+nCross+" "+nBonds+" "+nAngles+" "+nTorsions+" "+nImpropers+" "+nRules);
  }	

    void param_Read(String file, int ftype){
	try {
	    BufferedReader br = new BufferedReader(new FileReader(file));
	    String card;
	    int dataType=0;
	    while((card=br.readLine()) != null){
		if(card.startsWith("Atoms")){
		    dataType=1;
		} else if(card.startsWith("Cross")){
		    dataType=2;
		} else if(card.startsWith("Bonds")){
		    dataType=3;
		} else if(card.startsWith("Angles")){
		    dataType=4;
		} else if(card.startsWith("Proper dihedrals")){
		    dataType=5;
		} else if(card.startsWith("Improper dihedrals")){
		    dataType=6;
		} else if(card.startsWith("Atom types")){
		    dataType=7;
		} else if(card.startsWith("#")){
		} else {
		    if(dataType==0){
		    } else if(dataType==1){
			Atoms[nAtoms][0]=card; card=br.readLine(); Atoms[nAtoms][1]=card; idAtoms[nAtoms]=ftype; nAtoms++;
			if(nAtoms>=mAtoms) { System.out.println("Increase mAtoms"); setVisible(false); };
		    } else if(dataType==2){
			Atoms[nCross][0]=card; idCross[nCross]=ftype; nCross++; 
			if(nCross>=mCross) { System.out.println("Increase mCross"); setVisible(false); };
		    } else if(dataType==3){
			Atoms[nBonds][0]=card; idBonds[nBonds]=ftype; nBonds++; 
			if(nBonds>=mBonds) { System.out.println("Increase mBonds"); setVisible(false); };
		    } else if(dataType==4){
			Atoms[nAngles][0]=card; idAngles[nAngles]=ftype; nAngles++; 
			if(nAngles>=mAngles) { System.out.println("Increase mAngles"); setVisible(false); };
		    } else if(dataType==5){
			Atoms[nTorsions][0]=card; idTorsions[nTorsions]=ftype; nTorsions++; 
			if(nTorsions>=mTorsions) { System.out.println("Increase mTorsions"); setVisible(false); };
		    } else if(dataType==6){
			Atoms[nImpropers][0]=card; idImpropers[nImpropers]=ftype; nImpropers++; 
			if(nImpropers>=mImpropers) { System.out.println("Increase mImpropers"); setVisible(false); };
		    } else if(dataType==7){
			Atoms[nRules][0]=card; idRules[nRules]=ftype; nRules++; 
			if(nRules>=mRules) { System.out.println("Increase mRules"); setVisible(false); };
		    } else {
		    };
		};
    System.out.println(nAtoms+" "+nCross+" "+nBonds+" "+nAngles+" "+nTorsions+" "+nImpropers+" "+nRules+" "+card);
	    };
	    br.close();
	} catch(Exception ee) {ee.printStackTrace();};
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
