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

    char[] file_type = { '?', 's', 'x', 'u', 't'};

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

    JPanel header = new JPanel();
    JPanel dataPane = new JPanel();

    JLabel data00Label = new JLabel();
    JLabel data11Label = new JLabel();
    JLabel data21Label = new JLabel();
    JLabel data31Label = new JLabel();
    JLabel data41Label = new JLabel();
    JLabel data12Label = new JLabel();
    JLabel data22Label = new JLabel();
    JLabel data32Label = new JLabel();
    JLabel data42Label = new JLabel();
    JLabel data13Label = new JLabel();
    JLabel data23Label = new JLabel();
    JLabel data33Label = new JLabel();
    JLabel data43Label = new JLabel();
    JLabel data14Label = new JLabel();
    JLabel data24Label = new JLabel();
    JLabel data34Label = new JLabel();
    JLabel data44Label = new JLabel();
    JLabel data15Label = new JLabel();
    JLabel data25Label = new JLabel();
    JLabel data35Label = new JLabel();
    JLabel data45Label = new JLabel();

  public nwchem_Param(){

    super("Parameters");

    defaultFont = new Font("Dialog", Font.BOLD,12);

    super.getContentPane().setLayout(new GridBagLayout());
    super.getContentPane().setForeground(Color.black);
    super.getContentPane().setBackground(Color.lightGray);
    super.getContentPane().setFont(defaultFont);
    super.addWindowListener(this);

    header.setLayout(new GridBagLayout());
    header.setForeground(Color.black);
    header.setBackground(Color.lightGray);
    addComponent(super.getContentPane(),header,0,0,2,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.WEST);

    dataPane.setLayout(new GridBagLayout());
    dataPane.setForeground(Color.black);
    dataPane.setBackground(Color.lightGray);

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
    addComponent(header,dataPane,0,1,6,1,1,1,
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
    addComponent(header,doneButton,7,1,1,1,1,1,
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
    setSize(1000,600);
    setVisible(true);

  }	

    void param_Read(String file, int ftype){
	boolean repeat = false;
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
			Torsions[nTorsions] = new Torsion(card,ftype);
			while(Torsions[nTorsions].repeat){
			    card=br.readLine(); Torsions[nTorsions].add(card);
			}; 
			nTorsions++;
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
		    if(Atoms[i].type.equals(Atoms[j].type)) {
			Atoms[i].redefining=true; Atoms[j].redefined=true;
		    };
		};
	    }
	};
	if(nBonds>1){
	    for(int i=1; i<nBonds; i++){
		for(int j=0; j<i; j++){
		    if(Bonds[i].type1.equals(Bonds[j].type1) && Bonds[i].type2.equals(Bonds[j].type2)) {
			Bonds[i].redefining=true; Bonds[j].redefined=true;
		    } else if(Bonds[i].type1.equals(Bonds[j].type2) && Bonds[i].type2.equals(Bonds[j].type1)) {
			Bonds[i].redefining=true; Bonds[j].redefined=true;
		    };
		};
	    }
	};
	if(nAngles>1){
	    for(int i=1; i<nAngles; i++){
		for(int j=0; j<i; j++){
		    if(Angles[i].type1.equals(Angles[j].type1) && Angles[i].type2.equals(Angles[j].type2) && Angles[i].type3.equals(Angles[j].type3)) {
			Angles[i].redefining=true; Angles[j].redefined=true;
		    } else if(Angles[i].type1.equals(Angles[j].type3) && Angles[i].type2.equals(Angles[j].type2) && Angles[i].type3.equals(Angles[j].type1)) {
			Angles[i].redefining=true;Angles[j].redefined=true;
		    };
		};
	    }
	};
	if(nTorsions>1){
	    for(int i=1; i<nTorsions; i++){
		for(int j=0; j<i; j++){
		    if(Torsions[i].type2.equals(Torsions[j].type2) && Torsions[i].type3.equals(Torsions[j].type3)) {
			if(Torsions[i].type1.equals(Torsions[j].type1) && Torsions[i].type4.equals(Torsions[j].type4)) {
			    Torsions[i].redefining=true; Torsions[j].redefined=true;
			};
			if(Torsions[i].type1.equals("  ") && Torsions[i].type4.equals(Torsions[j].type4)) {
			    Torsions[i].redefining=true; Torsions[j].redefined=true;
			};
			if(Torsions[i].type1.equals(Torsions[j].type1) && Torsions[i].type4.equals("  ")) {
			    Torsions[i].redefining=true; Torsions[j].redefined=true;
			};
			if(Torsions[i].type1.equals("  ") && Torsions[i].type4.equals("  ")) {
			    Torsions[i].redefining=true; Torsions[j].redefined=true;
			};
		    } else if(Torsions[i].type2.equals(Torsions[j].type3) && Torsions[i].type3.equals(Torsions[j].type2)) {
			if(Torsions[i].type1.equals(Torsions[j].type4) && Torsions[i].type4.equals(Torsions[j].type1)) {
			    Torsions[i].redefining=true; Torsions[j].redefined=true;
			};
			if(Torsions[i].type1.equals("  ") && Torsions[i].type4.equals(Torsions[j].type1)) {
			    Torsions[i].redefining=true; Torsions[j].redefined=true;
			};
			if(Torsions[i].type1.equals(Torsions[j].type4) && Torsions[i].type4.equals("  ")) {
			    Torsions[i].redefining=true; Torsions[j].redefined=true;
			};
			if(Torsions[i].type1.equals("  ") && Torsions[i].type4.equals("  ")) {
			    Torsions[i].redefining=true; Torsions[j].redefined=true;
			};
		    };
		};
	    }
	};
	if(nImpropers>1){
	    for(int i=1; i<nImpropers; i++){
		for(int j=0; j<i; j++){
		    if(Impropers[i].type3.equals(Impropers[j].type3) && Impropers[i].type4.equals(Impropers[j].type4)) {
			if(Impropers[i].type1.equals(Impropers[j].type1) && Impropers[i].type2.equals(Impropers[j].type2)) {
			    Impropers[i].redefining=true; Impropers[j].redefined=true;
			};
			if(Impropers[i].type1.equals("  ") && Impropers[i].type2.equals(Impropers[j].type2)) {
			    Impropers[i].redefining=true; Impropers[j].redefined=true;
			};
			if(Impropers[i].type1.equals(Impropers[j].type1) && Impropers[i].type2.equals("  ")) {
			    Impropers[i].redefining=true; Impropers[j].redefined=true;
			};
			if(Impropers[i].type1.equals("  ") && Impropers[i].type2.equals("  ")) {
			    Impropers[i].redefining=true; Impropers[j].redefined=true;
			};
		    };
		};
	    }
	};
	if(nRules>1){
	    for(int i=1; i<nRules; i++){
		for(int j=0; j<i; j++){
		    if(Rules[i].type.equals(Rules[j].type)) {
			Rules[i].redefining=true; Rules[j].redefined=true; 
		    };
		};
	    }
	};
    }

    void displayEntries(){
	if(nAtoms>1){
	    cellRender aRender = new cellRender();
	    aList.setCellRenderer(aRender);
	    for(int i=0; i<nAtoms; i++){
		if(Atoms[i].redefined){aRender.setInactive(atomList.size());};
		if(Atoms[i].redefining){aRender.setReplacing(atomList.size());};
		atomList.addElement(file_type[Atoms[i].source]+": "+Atoms[i].type);
	    };
	};
	if(nCross>1){
	    cellRender cRender = new cellRender();
	    cList.setCellRenderer(cRender);
	    for(int i=0; i<nCross; i++){
		if(Cross[i].redefined){cRender.setInactive(crossList.size());};
		if(Cross[i].redefining){cRender.setReplacing(crossList.size());};
		crossList.addElement(file_type[Cross[i].source]+": "+Cross[i].type1+"-"+Cross[i].type2);
	    }
	};
	if(nBonds>1){
	    cellRender bRender = new cellRender();
	    bList.setCellRenderer(bRender);
	    for(int i=0; i<nBonds; i++){
		if(Bonds[i].redefined){bRender.setInactive(bondList.size());};
		if(Bonds[i].redefining){bRender.setReplacing(bondList.size());};
		bondList.addElement(file_type[Bonds[i].source]+": "+Bonds[i].type1+"-"+Bonds[i].type2);
	    }
	};
	if(nAngles>1){
	    cellRender hRender = new cellRender();
	    hList.setCellRenderer(hRender);
	    for(int i=0; i<nAngles; i++){
		if(Angles[i].redefined){hRender.setInactive(angleList.size());};
		if(Angles[i].redefining){hRender.setReplacing(angleList.size());};
		angleList.addElement(file_type[Angles[i].source]+": "+Angles[i].type1+"-"+Angles[i].type2+"-"+Angles[i].type3);
	    }
	};
	if(nTorsions>1){
	    cellRender tRender = new cellRender();
	    tList.setCellRenderer(tRender);
	    for(int i=0; i<nTorsions; i++){
		if(Torsions[i].redefined){tRender.setInactive(torsionList.size());};
		if(Torsions[i].redefining){tRender.setReplacing(torsionList.size());};
		torsionList.addElement(file_type[Torsions[i].source]+": "+Torsions[i].type1+"-"+Torsions[i].type2+"-"+Torsions[i].type3+"-"+Torsions[i].type4);
	    }
	};
	if(nImpropers>1){
	    cellRender iRender = new cellRender();
	    iList.setCellRenderer(iRender);
	    for(int i=0; i<nImpropers; i++){
		if(Impropers[i].redefined){iRender.setInactive(improperList.size());};
		if(Impropers[i].redefining){iRender.setReplacing(improperList.size());};
		improperList.addElement(file_type[Impropers[i].source]+": "+Impropers[i].type1+"-"+Impropers[i].type2+"-"+Impropers[i].type3+"-"+Impropers[i].type4);
	    }
	};
	if(nRules>1){
	    cellRender rRender = new cellRender();
	    rList.setCellRenderer(rRender);
	    for(int i=0; i<nRules; i++){
		if(Rules[i].redefined){rRender.setInactive(ruleList.size());};
		if(Rules[i].redefining){rRender.setReplacing(ruleList.size());};
		ruleList.addElement(file_type[Rules[i].source]+": "+Rules[i].type);
	    }
	};
    }

    class cellRender extends JLabel implements ListCellRenderer{

	boolean[] inactive = new boolean[1000];
	boolean[] replacing = new boolean[1000];

	public cellRender(){
	    for(int i=0; i<1000; i++) { inactive[i]=false; replacing[i]=false; };
	}

	public Component getListCellRendererComponent(JList list, Object value, int index, boolean isSelected, boolean cellHasFocus){
	    if(value!=null){ String text=value.toString(); setText(text); };
	    if(isSelected){
		if(inactive[index]){
		    setForeground(Color.red);
		    setBackground(list.getSelectionBackground());
		} else if(replacing[index]) {
		    setForeground(Color.blue);
		    setBackground(list.getSelectionBackground());
		} else {
		    setForeground(list.getSelectionForeground());
		    setBackground(list.getSelectionBackground());
		};
	    } else {
		if(inactive[index]){
		    setForeground(Color.red);
		    setBackground(list.getBackground());
		} else if(replacing[index]) {
		    setForeground(Color.blue);
		    setBackground(list.getBackground());
		} else {
		    setForeground(list.getForeground());
		    setBackground(list.getBackground());
		};
	    };
	    return this;
	}

	public void setInactive(int num){
	    inactive[num]=true;
	}

	public void setReplacing(int num){
	    replacing[num]=true;
	}
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
	if(mouse.getModifiers()==MouseEvent.BUTTON1_MASK){
	    if(mouse.getSource()==aList){
		data00Label.setText("Atom ");
		data11Label.setText("Type "+file_type[Atoms[aList.getSelectedIndex()].source]+": "+Atoms[aList.getSelectedIndex()].type);
		data12Label.setText("Atom number "+Atoms[aList.getSelectedIndex()].number);
		data21Label.setText("Mass");
		data22Label.setText(String.valueOf(Atoms[aList.getSelectedIndex()].mass));
		data31Label.setText("Epsilon");
		data32Label.setText(String.valueOf(Atoms[aList.getSelectedIndex()].epsilon));
		data33Label.setText(",  "+String.valueOf(Atoms[aList.getSelectedIndex()].epsilon3));
		data41Label.setText("R*");
		data42Label.setText(String.valueOf(Atoms[aList.getSelectedIndex()].rstar));
		data43Label.setText(",  "+String.valueOf(Atoms[aList.getSelectedIndex()].rstar3));
		dataPane.removeAll();
		addComponent(dataPane,data00Label,0,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data11Label,1,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data12Label,1,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data21Label,1,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data22Label,2,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data31Label,1,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data32Label,2,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data33Label,3,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data41Label,1,4,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data42Label,2,4,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data43Label,3,4,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		dataPane.validate();
	    };
	    if(mouse.getSource()==bList){
		data00Label.setText("Bond ");
		data11Label.setText("Type "+file_type[Bonds[bList.getSelectedIndex()].source]+": "+Bonds[bList.getSelectedIndex()].type1+"-"+Bonds[bList.getSelectedIndex()].type2);
		data21Label.setText("Length");
		data22Label.setText(String.valueOf(Bonds[bList.getSelectedIndex()].value));
		data31Label.setText("Force constant");
		data32Label.setText(String.valueOf(Bonds[bList.getSelectedIndex()].force));
		data41Label.setText("Charge");
		data42Label.setText(String.valueOf(Bonds[bList.getSelectedIndex()].charge));
		dataPane.removeAll();
		addComponent(dataPane,data00Label,0,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data11Label,1,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data21Label,1,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data22Label,2,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data31Label,1,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data32Label,2,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data41Label,1,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data42Label,2,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		dataPane.validate();
	    };
	    if(mouse.getSource()==hList){
		data00Label.setText("Angle ");
		data11Label.setText("Type "+file_type[Angles[hList.getSelectedIndex()].source]+": "+Angles[hList.getSelectedIndex()].type1+"-"+Angles[hList.getSelectedIndex()].type2+"-"+Angles[hList.getSelectedIndex()].type3);
		data21Label.setText("Angle");
		data22Label.setText(String.valueOf(Angles[hList.getSelectedIndex()].value));
		data31Label.setText("Force constant");
		data32Label.setText(String.valueOf(Angles[hList.getSelectedIndex()].force));
		dataPane.removeAll();
		addComponent(dataPane,data00Label,0,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data11Label,1,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data21Label,1,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data22Label,2,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data31Label,1,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data32Label,2,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		dataPane.validate();
	    };
	    if(mouse.getSource()==tList){
		data00Label.setText("Torsion ");
		data11Label.setText("Type "+file_type[Torsions[tList.getSelectedIndex()].source]+": "+Torsions[tList.getSelectedIndex()].type1+"-"+Torsions[tList.getSelectedIndex()].type2+"-"+Torsions[tList.getSelectedIndex()].type3+"-"+Torsions[tList.getSelectedIndex()].type4);
		data21Label.setText("Angle");
		data22Label.setText(String.valueOf(Torsions[tList.getSelectedIndex()].value[0]));
		data31Label.setText("Force constant");
		data32Label.setText(String.valueOf(Torsions[tList.getSelectedIndex()].force[0]));
		data41Label.setText("Multiplicity");
		data42Label.setText(String.valueOf(Torsions[tList.getSelectedIndex()].multiplicity[0]));
		dataPane.removeAll();
		addComponent(dataPane,data00Label,0,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data11Label,1,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data21Label,1,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data22Label,2,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data31Label,1,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data32Label,2,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data41Label,1,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data42Label,2,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		if(Torsions[tList.getSelectedIndex()].number>0){
		    data23Label.setText(", "+String.valueOf(Torsions[tList.getSelectedIndex()].value[1]));
		    data33Label.setText(", "+String.valueOf(Torsions[tList.getSelectedIndex()].force[1]));
		    data43Label.setText(", "+String.valueOf(Torsions[tList.getSelectedIndex()].multiplicity[1]));
		    addComponent(dataPane,data23Label,3,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		    addComponent(dataPane,data33Label,3,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		    addComponent(dataPane,data43Label,3,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		};
		if(Torsions[tList.getSelectedIndex()].number>1){
		    data24Label.setText(", "+String.valueOf(Torsions[tList.getSelectedIndex()].value[2]));
		    data34Label.setText(", "+String.valueOf(Torsions[tList.getSelectedIndex()].force[2]));
		    data44Label.setText(", "+String.valueOf(Torsions[tList.getSelectedIndex()].multiplicity[2]));
		    addComponent(dataPane,data24Label,4,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		    addComponent(dataPane,data34Label,4,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		    addComponent(dataPane,data44Label,4,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		};
		if(Torsions[tList.getSelectedIndex()].number>2){
		    data25Label.setText(", "+String.valueOf(Torsions[tList.getSelectedIndex()].value[3]));
		    data35Label.setText(", "+String.valueOf(Torsions[tList.getSelectedIndex()].force[3]));
		    data45Label.setText(", "+String.valueOf(Torsions[tList.getSelectedIndex()].multiplicity[3]));
		    addComponent(dataPane,data25Label,5,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		    addComponent(dataPane,data35Label,5,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		    addComponent(dataPane,data45Label,5,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		};
		dataPane.validate();
	    };
	    if(mouse.getSource()==iList){
		data00Label.setText("Improper ");
		data11Label.setText("Type "+file_type[Impropers[iList.getSelectedIndex()].source]+": "+Impropers[iList.getSelectedIndex()].type1+"-"+Impropers[iList.getSelectedIndex()].type2+"-"+Impropers[iList.getSelectedIndex()].type3+"-"+Impropers[iList.getSelectedIndex()].type4);
		data21Label.setText("Angle");
		data22Label.setText(String.valueOf(Impropers[iList.getSelectedIndex()].value));
		data31Label.setText("Force constant");
		data32Label.setText(String.valueOf(Impropers[iList.getSelectedIndex()].force));
		dataPane.removeAll();
		addComponent(dataPane,data00Label,0,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data11Label,1,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data21Label,1,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data22Label,2,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data31Label,1,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data32Label,2,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		dataPane.validate();
	    };
	};
    }
    
    public void mouseEntered(MouseEvent mouse){
    }

    public void mouseExited(MouseEvent mouse){
    }
  
}
