import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

class nwchem_Param extends JFrame implements ActionListener, ChangeListener, WindowListener, MouseListener {
  
  Font defaultFont;

    String ffield = "unknown";
    String file_s = " ";
    String file_x = " ";
    String file_u = " ";
    String file_t = " ";

    int mAtoms=100, mCross=100, mBonds=250, mAngles=500, mTorsions=200, mImpropers=100, mRules=100;

    Atom[] Atoms = new Atom[mAtoms];
    Cros[] Cross = new Cros[mCross];
    Bond[] Bonds = new Bond[mBonds];
    Angle[] Angles = new Angle[mAngles];
    Torsion[] Torsions = new Torsion[mTorsions];
    Improper[] Impropers = new Improper[mImpropers];
    Rule[] Rules = new Rule[mRules];

    int nAtoms=0, nCross=0, nBonds=0, nAngles=0, nTorsions=0, nImpropers=0, nRules=0;

    DefaultListModel atomList = new DefaultListModel();
    JList aList = new JList(atomList);
    JScrollPane atomPane = new JScrollPane(aList);
    DefaultListModel crossList = new DefaultListModel();
    JList cList = new JList(crossList);
    JScrollPane crossPane = new JScrollPane(cList);
    DefaultListModel bondList = new DefaultListModel();
    JList bList = new JList(bondList);
    JScrollPane bondPane = new JScrollPane(bList);
    DefaultListModel angleList = new DefaultListModel();
    JList hList = new JList(angleList);
    JScrollPane anglePane = new JScrollPane(hList);
    DefaultListModel torsionList = new DefaultListModel();
    JList tList = new JList(torsionList);
    JScrollPane torsionPane = new JScrollPane(tList);
    DefaultListModel improperList = new DefaultListModel();
    JList iList = new JList(improperList);
    JScrollPane improperPane = new JScrollPane(iList);
    DefaultListModel ruleList = new DefaultListModel();
    JList rList = new JList(ruleList);
    JScrollPane rulePane = new JScrollPane(rList);

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

    addComponent(header,atomPane,0,0,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    addComponent(header,crossPane,1,0,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    addComponent(header,bondPane,2,0,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    addComponent(header,anglePane,3,0,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    addComponent(header,torsionPane,4,0,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    addComponent(header,improperPane,5,0,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    addComponent(header,rulePane,6,0,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.WEST);
    aList.addMouseListener(this);
    cList.addMouseListener(this);
    bList.addMouseListener(this);
    hList.addMouseListener(this);
    tList.addMouseListener(this);
    iList.addMouseListener(this);
    rList.addMouseListener(this);
    aList.setVisibleRowCount(15);
    cList.setVisibleRowCount(15);
    bList.setVisibleRowCount(15);
    hList.setVisibleRowCount(15);
    tList.setVisibleRowCount(15);
    iList.setVisibleRowCount(15);
    rList.setVisibleRowCount(15);

    JButton doneButton = new JButton("Done");
    addComponent(header,doneButton,0,3,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    doneButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	setVisible(false); }});

    try{
      BufferedReader br = new BufferedReader(new FileReader("/usr/people/d3j191/.nwchemrc"));
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
    checkRedefinitions();
    displayEntries();

    setLocation(25,225);	
    setSize(900,700);
    setVisible(true);

  }	

    void param_Read(String file, int ftype){
	try {
	    BufferedReader br = new BufferedReader(new FileReader(file));
	    String card;
	    String card2;
	    String card3;
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
			card2=br.readLine(); Atoms[nAtoms] = new Atom(card,card2,ftype); nAtoms++;
			if(nAtoms>=mAtoms) { System.out.println("Increase mAtoms"); setVisible(false); };
		    } else if(dataType==2){
			card2=br.readLine(); Cross[nCross] = new Cros(card,card2,ftype); nCross++;
			if(nCross>=mCross) { System.out.println("Increase mCross"); setVisible(false); };
		    } else if(dataType==3){
			Bonds[nBonds] = new Bond(card,ftype); nBonds++;
			if(nBonds>=mBonds) { System.out.println("Increase mBonds"); setVisible(false); };
		    } else if(dataType==4){
			Angles[nAngles] = new Angle(card,ftype); nAngles++;
			if(nAngles>=mAngles) { System.out.println("Increase mAngles"); setVisible(false); };
		    } else if(dataType==5){
			Torsions[nTorsions] = new Torsion(card,ftype); nTorsions++;
			if(nTorsions>=mTorsions) { System.out.println("Increase mTorsions"); setVisible(false); };
		    } else if(dataType==6){
			Impropers[nImpropers] = new Improper(card,ftype); nImpropers++;
			if(nImpropers>=mImpropers) { System.out.println("Increase mImpropers"); setVisible(false); };
		    } else if(dataType==7){
			card2=br.readLine(); card3=br.readLine(); Rules[nRules] = new Rule(card,card2,card3,ftype); nRules++;
			if(nRules>=mRules) { System.out.println("Increase mRules"); setVisible(false); };
		    } else {
		    };
		};
	    };
	    br.close();
	} catch(Exception ee) {ee.printStackTrace();};
    }

    void checkRedefinitions(){
	if(nAtoms>1){
	    for(int i=1; i<nAtoms; i++){
		for(int j=0; j<i; j++){
		    if(Atoms[i].type.equals(Atoms[j].type)) Atoms[j].redefined=true;
		};
	    }
	};
	if(nBonds>1){
	    for(int i=1; i<nBonds; i++){
		for(int j=0; j<i; j++){
		    if(Bonds[i].type1.equals(Bonds[j].type1) && Bonds[i].type2.equals(Bonds[j].type2)) {
			Bonds[j].redefined=true;
		    } else if(Bonds[i].type1.equals(Bonds[j].type2) && Bonds[i].type2.equals(Bonds[j].type1)) {
			Bonds[j].redefined=true;
		    };
		};
	    }
	};
	if(nAngles>1){
	    for(int i=1; i<nAngles; i++){
		for(int j=0; j<i; j++){
		    if(Angles[i].type1.equals(Angles[j].type1) && Angles[i].type2.equals(Angles[j].type2) && Angles[i].type3.equals(Angles[j].type3)) {
			Angles[j].redefined=true;
		    } else if(Angles[i].type1.equals(Angles[j].type3) && Angles[i].type2.equals(Angles[j].type2) && Angles[i].type3.equals(Angles[j].type1)) {
			Angles[j].redefined=true;
		    };
		};
	    }
	};
	if(nTorsions>1){
	    for(int i=1; i<nTorsions; i++){
		for(int j=0; j<i; j++){
		    if(Torsions[i].type2.equals(Torsions[j].type2) && Torsions[i].type3.equals(Torsions[j].type3)) {
			if(Torsions[i].type1.equals(Torsions[j].type1) && Torsions[i].type4.equals(Torsions[j].type4)) Torsions[j].redefined=true;
			if(Torsions[i].type1.equals("  ") && Torsions[i].type4.equals(Torsions[j].type4)) Torsions[j].redefined=true;
			if(Torsions[i].type1.equals(Torsions[j].type1) && Torsions[i].type4.equals("  ")) Torsions[j].redefined=true;
			if(Torsions[i].type1.equals("  ") && Torsions[i].type4.equals("  ")) Torsions[j].redefined=true;
		    } else if(Torsions[i].type2.equals(Torsions[j].type3) && Torsions[i].type3.equals(Torsions[j].type2)) {
			if(Torsions[i].type1.equals(Torsions[j].type4) && Torsions[i].type4.equals(Torsions[j].type1)) Torsions[j].redefined=true;
			if(Torsions[i].type1.equals("  ") && Torsions[i].type4.equals(Torsions[j].type1)) Torsions[j].redefined=true;
			if(Torsions[i].type1.equals(Torsions[j].type4) && Torsions[i].type4.equals("  ")) Torsions[j].redefined=true;
			if(Torsions[i].type1.equals("  ") && Torsions[i].type4.equals("  ")) Torsions[j].redefined=true;
		    };
		};
	    }
	};
	if(nImpropers>1){
	    for(int i=1; i<nImpropers; i++){
		for(int j=0; j<i; j++){
		    if(Impropers[i].type3.equals(Impropers[j].type3) && Impropers[i].type4.equals(Impropers[j].type4)) {
			if(Impropers[i].type1.equals(Impropers[j].type1) && Impropers[i].type2.equals(Impropers[j].type2)) Impropers[j].redefined=true;
			if(Impropers[i].type1.equals("  ") && Impropers[i].type2.equals(Impropers[j].type2)) Impropers[j].redefined=true;
			if(Impropers[i].type1.equals(Impropers[j].type1) && Impropers[i].type2.equals("  ")) Impropers[j].redefined=true;
			if(Impropers[i].type1.equals("  ") && Impropers[i].type2.equals("  ")) Impropers[j].redefined=true;
		    };
		};
	    }
	};
    }

    void displayEntries(){
	if(nAtoms>1){
	    for(int i=0; i<nAtoms; i++){
		atomList.addElement(Atoms[i].type);
		if(Atoms[i].redefined){
		    System.out.println("Atom "+i+" "+aList.getSelectedIndex());
		};
	    };
	};
	if(nCross>1){
	    for(int i=0; i<nCross; i++){
		crossList.addElement(Cross[i].type1+"-"+Cross[i].type2);
		if(Cross[i].redefined){
		    System.out.println("Cross "+i+" "+cList.getSelectedIndex());
		};
	    }
	};
	if(nBonds>1){
	    for(int i=0; i<nBonds; i++){
		bondList.addElement(Bonds[i].type1+"-"+Bonds[i].type2);
		if(Bonds[i].redefined){
		    //		    JList.AccessibleJList.AccessibleJListChild(aList,i).setBackground(Color.blue);
		    System.out.println("Bond "+i+" "+bList.getSelectedIndex());
		};
	    }
	};
	if(nAngles>1){
	    for(int i=0; i<nAngles; i++){
		angleList.addElement(Angles[i].type1+"-"+Angles[i].type2+"-"+Angles[i].type3);
		if(Angles[i].redefined){
		    System.out.println("Angle "+i+" "+hList.getSelectedIndex());
		};
	    }
	};
	if(nTorsions>1){
	    for(int i=0; i<nTorsions; i++){
		torsionList.addElement(Torsions[i].type1+"-"+Torsions[i].type2+"-"+Torsions[i].type3+"-"+Torsions[i].type4);
		if(Torsions[i].redefined){
		    System.out.println("Torsion  "+i+" "+tList.getSelectedIndex());
		};
	    }
	};
	if(nImpropers>1){
	    for(int i=0; i<nImpropers; i++){
		improperList.addElement(Impropers[i].type1+"-"+Impropers[i].type2+"-"+Impropers[i].type3+"-"+Impropers[i].type4);
		if(Impropers[i].redefined){
		    System.out.println("Improper "+i+" "+iList.getSelectedIndex());
		};
	    }
	};
	if(nRules>1){
	    for(int i=0; i<nRules; i++){
		ruleList.addElement(Rules[i].type);
		if(Rules[i].redefined){
		    System.out.println("Rule "+i+" "+rList.getSelectedIndex());
		};
	    }
	};
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

    public void mousePressed(MouseEvent mouse){}

    public void mouseReleased(MouseEvent mouse){
    }

    public void mouseEntered(MouseEvent mouse){
    }

    public void mouseExited(MouseEvent mouse){
    }
  
}
