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

    int[] indexA = new int[mAtoms];
    int[] indexC = new int[mCross];
    int[] indexB = new int[mBonds];
    int[] indexH = new int[mAngles];
    int[] indexT = new int[mTorsions];
    int[] indexI = new int[mImpropers];
    int[] indexR = new int[mRules];

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

    JLabel header1 = new JLabel(" Atom Type ");
    JLabel header2 = new JLabel("   Cross   ");
    JLabel header3 = new JLabel("   Bond    ");
    JLabel header4 = new JLabel("   Angle   ");
    JLabel header5 = new JLabel("  Torsion  ");
    JLabel header6 = new JLabel("  Improper ");
    JLabel header7 = new JLabel(" Type Rule ");

    JLabel data00Label = new JLabel();
    JLabel data11Label = new JLabel();
    JLabel data21Label = new JLabel();
    JLabel data31Label = new JLabel();
    JLabel data41Label = new JLabel();
    JLabel data51Label = new JLabel();
    JLabel data12Label = new JLabel();
    JLabel data22Label = new JLabel();
    JLabel data32Label = new JLabel();
    JLabel data42Label = new JLabel();
    JLabel data52Label = new JLabel();
    JLabel data13Label = new JLabel();
    JLabel data23Label = new JLabel();
    JLabel data33Label = new JLabel();
    JLabel data43Label = new JLabel();
    JLabel data53Label = new JLabel();
    JLabel data14Label = new JLabel();
    JLabel data24Label = new JLabel();
    JLabel data34Label = new JLabel();
    JLabel data44Label = new JLabel();
    JLabel data15Label = new JLabel();
    JLabel data25Label = new JLabel();
    JLabel data35Label = new JLabel();
    JLabel data45Label = new JLabel();
    JLabel data16Label = new JLabel();
    JLabel data26Label = new JLabel();
    JLabel data36Label = new JLabel();
    JLabel data46Label = new JLabel();

    JButton selectButton = new JButton("Select Atom");
    boolean atomSelected = false;
    String selectedAtom; 

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
		 GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);

    dataPane.setLayout(new GridBagLayout());
    dataPane.setForeground(Color.black);
    dataPane.setBackground(Color.lightGray);

    selectButton.addActionListener(new ActionListener(){
	    public void actionPerformed(ActionEvent e){
		if(atomSelected) {selectAll(); } else { selectAtom();}; }});

    addComponent(header,header1,0,0,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    addComponent(header,header2,1,0,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    addComponent(header,header3,2,0,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    addComponent(header,header4,3,0,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    addComponent(header,header5,4,0,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    addComponent(header,header6,5,0,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    addComponent(header,header7,6,0,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    addComponent(header,dataPane,0,2,5,1,1,1,
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
    data00Label.setText("Entries in red are redefined by entries in blue");
    addComponent(dataPane,data00Label,0,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);

    JButton doneButton = new JButton("Done");
    addComponent(header,doneButton,6,2,1,1,1,1,
		 GridBagConstraints.NONE,GridBagConstraints.CENTER);
    doneButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){ 
	setVisible(false); }});

    try{
      BufferedReader br = new BufferedReader(new FileReader("./.nwchemrc"));
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
    //System.out.println(nAtoms+" "+nCross+" "+nBonds+" "+nAngles+" "+nTorsions+" "+nImpropers+" "+nRules);
    checkRedefinitions();
    displayEntries();

    setLocation(25,225);	
    setSize(1000,800);
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
			if(card.length()>0){
			    if(card.charAt(0)!=' '){
				card2=br.readLine(); card3=br.readLine();
				Rules[nRules] = new Rule(card,card2,card3,ftype); nRules++;
				if(nRules>=mRules) { System.out.println("Increase mRules"); setVisible(false); };
			    };
			};
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
		    }; 
		    if(Torsions[i].type2.equals(Torsions[j].type3) && Torsions[i].type3.equals(Torsions[j].type2)) {
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
	atomList.removeAllElements();
	crossList.removeAllElements();
	bondList.removeAllElements();
	angleList.removeAllElements();
	torsionList.removeAllElements();
	improperList.removeAllElements();
	ruleList.removeAllElements();
	header.remove(atomPane);
	header.remove(crossPane);
	header.remove(bondPane);
	header.remove(anglePane);
	header.remove(torsionPane);
	header.remove(improperPane);
	header.remove(rulePane);
	int j;
	if(nAtoms>1){
	    cellRender aRender = new cellRender();
	    aList.setCellRenderer(aRender);
	    j=0;
	    for(int i=0; i<nAtoms; i++){
		if(Atoms[i].selected){
		    if(Atoms[i].redefined){aRender.setInactive(atomList.size());};
		    if(Atoms[i].redefining){aRender.setReplacing(atomList.size());};
		    atomList.addElement(file_type[Atoms[i].source]+": "+Atoms[i].type);
		    indexA[j]=i; j++;
		};
	    };
	    if(j>0) addComponent(header,atomPane,0,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.CENTER);
	};
	if(nCross>1){
	    cellRender cRender = new cellRender();
	    cList.setCellRenderer(cRender);
	    j=0;
	    for(int i=0; i<nCross; i++){
		if(Cross[i].selected){
		    if(Cross[i].redefined){cRender.setInactive(crossList.size());};
		    if(Cross[i].redefining){cRender.setReplacing(crossList.size());};
		    crossList.addElement(file_type[Cross[i].source]+": "+Cross[i].type1+"-"+Cross[i].type2);
		    indexC[j]=i; j++;
		};
	    };
	    if(j>0) addComponent(header,crossPane,1,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.CENTER); 
	};
	if(nBonds>1){
	    cellRender bRender = new cellRender();
	    bList.setCellRenderer(bRender);
	    j=0;
	    for(int i=0; i<nBonds; i++){
		if(Bonds[i].selected){
		    if(Bonds[i].redefined){bRender.setInactive(bondList.size());};
		    if(Bonds[i].redefining){bRender.setReplacing(bondList.size());};
		    bondList.addElement(file_type[Bonds[i].source]+": "+Bonds[i].type1+"-"+Bonds[i].type2);
		    indexB[j]=i; j++;
		};
	    };
	    if(j>0) addComponent(header,bondPane,2,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.CENTER);
	};
	if(nAngles>1){
	    cellRender hRender = new cellRender();
	    hList.setCellRenderer(hRender);
	    j=0;
	    for(int i=0; i<nAngles; i++){
		if(Angles[i].selected){
		    if(Angles[i].redefined){hRender.setInactive(angleList.size());};
		    if(Angles[i].redefining){hRender.setReplacing(angleList.size());};
		    angleList.addElement(file_type[Angles[i].source]+": "+Angles[i].type1+"-"+Angles[i].type2+"-"+Angles[i].type3);
		    indexH[j]=i; j++;
		};
	    };
	    if(j>0) addComponent(header,anglePane,3,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.CENTER); 
	};
	if(nTorsions>1){
	    cellRender tRender = new cellRender();
	    tList.setCellRenderer(tRender);
	    j=0;
	    for(int i=0; i<nTorsions; i++){
		if(Torsions[i].selected){
		    if(Torsions[i].redefined){tRender.setInactive(torsionList.size());};
		    if(Torsions[i].redefining){tRender.setReplacing(torsionList.size());};
		    torsionList.addElement(file_type[Torsions[i].source]+": "+Torsions[i].type1+"-"+Torsions[i].type2+"-"+Torsions[i].type3+"-"+Torsions[i].type4);
		    indexT[j]=i; j++;
		};
	    };
	    if(j>0) addComponent(header,torsionPane,4,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.CENTER);
	};
	if(nImpropers>1){
	    cellRender iRender = new cellRender();
	    iList.setCellRenderer(iRender);
	    j=0;
	    for(int i=0; i<nImpropers; i++){
		if(Impropers[i].selected){
		    if(Impropers[i].redefined){iRender.setInactive(improperList.size());};
		    if(Impropers[i].redefining){iRender.setReplacing(improperList.size());};
		    improperList.addElement(file_type[Impropers[i].source]+": "+Impropers[i].type1+"-"+Impropers[i].type2+"-"+Impropers[i].type3+"-"+Impropers[i].type4);
		    indexI[j]=i; j++;
		};
	    };
	    if(j>0) addComponent(header,improperPane,5,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.CENTER);
	};
	if(nRules>1){
	    cellRender rRender = new cellRender();
	    rList.setCellRenderer(rRender);
	    j=0;
	    for(int i=0; i<nRules; i++){
		if(Rules[i].selected){
		    if(Rules[i].redefined){rRender.setInactive(ruleList.size());};
		    if(Rules[i].redefining){rRender.setReplacing(ruleList.size());};
		    ruleList.addElement(file_type[Rules[i].source]+": "+Rules[i].type);
		    indexR[j]=i; j++;
		};
	    };
	    if(j>0) addComponent(header,rulePane,6,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.CENTER);
	};
	header.validate();
    }

    void selectAll(){
	if(nAtoms>1){ for(int i=0; i<nAtoms; i++) Atoms[i].selected=true; };
	if(nCross>1){ for(int i=0; i<nCross; i++) Cross[i].selected=true; };
	if(nBonds>1){ for(int i=0; i<nBonds; i++) Bonds[i].selected=true; };
	if(nAngles>1){ for(int i=0; i<nAngles; i++) Angles[i].selected=true; };
	if(nTorsions>1){ for(int i=0; i<nTorsions; i++) Torsions[i].selected=true; };
	if(nImpropers>1){ for(int i=0; i<nImpropers; i++) Impropers[i].selected=true; };
	if(nRules>1){ for(int i=0; i<nRules; i++) Rules[i].selected=true; };
	selectButton.setText("Select Atom"); atomSelected=false; displayEntries();
    }

    void selectAtom(){
	if(nCross>1){
	    for(int i=0; i<nCross; i++){
		Cross[i].selected=Cross[i].type1.equals(selectedAtom) || Cross[i].type2.equals(selectedAtom);
	    }
	};
	if(nBonds>1){
	    for(int i=0; i<nBonds; i++){
		Bonds[i].selected=Bonds[i].type1.equals(selectedAtom) || Bonds[i].type2.equals(selectedAtom);
	    }
	};
	if(nAngles>1){
	    for(int i=0; i<nAngles; i++){
		Angles[i].selected=Angles[i].type1.equals(selectedAtom) || Angles[i].type2.equals(selectedAtom) || Angles[i].type3.equals(selectedAtom);
	    }
	};
	if(nTorsions>1){
	    for(int i=0; i<nTorsions; i++){
		Torsions[i].selected=Torsions[i].type1.equals(selectedAtom) || Torsions[i].type2.equals(selectedAtom) || Torsions[i].type3.equals(selectedAtom) || Torsions[i].type4.equals(selectedAtom);
	    }
	};
	if(nImpropers>1){
	    for(int i=0; i<nImpropers; i++){
		Impropers[i].selected=Impropers[i].type1.equals(selectedAtom) || Impropers[i].type2.equals(selectedAtom) || Impropers[i].type3.equals(selectedAtom) || Impropers[i].type4.equals(selectedAtom);
	    }
	};
	if(nRules>1){
	    for(int i=0; i<nRules; i++){
		Rules[i].selected=Rules[i].type.equals(selectedAtom);
	    }
	}; 
	atomSelected=true; selectButton.setText("Select All Atoms"); displayEntries();
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
	int j;
	if(mouse.getModifiers()==MouseEvent.BUTTON1_MASK){
	    if(mouse.getSource()==aList){
		j=indexA[aList.getSelectedIndex()];
		data00Label.setText("Atom ");
		data11Label.setText("Type ");
		data12Label.setText(file_type[Atoms[j].source]+": "+Atoms[j].type);
		data21Label.setText("Atom number ");
		data22Label.setText(String.valueOf(Atoms[j].number));
		data31Label.setText("Mass");
		data32Label.setText(String.valueOf(Atoms[j].mass));
		data41Label.setText("Epsilon");
		data42Label.setText(String.valueOf(Atoms[j].epsilon));
		data43Label.setText(",  "+String.valueOf(Atoms[j].epsilon3));
		data51Label.setText("R*");
		data52Label.setText(String.valueOf(Atoms[j].rstar));
		data53Label.setText(",  "+String.valueOf(Atoms[j].rstar3));
		dataPane.removeAll();
		addComponent(dataPane,data00Label,0,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data11Label,1,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data12Label,2,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data21Label,1,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data22Label,2,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data31Label,1,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data32Label,2,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data41Label,1,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data42Label,2,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data43Label,3,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data51Label,1,4,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data52Label,2,4,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data53Label,3,4,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(header,selectButton,6,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.CENTER);
		selectedAtom=Atoms[j].type;
		dataPane.validate();
	    };
	    if(mouse.getSource()==bList){
		j=indexB[bList.getSelectedIndex()];
		data00Label.setText("Bond ");
		data11Label.setText("Type "+file_type[Bonds[j].source]+": "+Bonds[j].type1+"-"+Bonds[j].type2);
		data21Label.setText("Length");
		data22Label.setText(String.valueOf(Bonds[j].value));
		data31Label.setText("Force constant");
		data32Label.setText(String.valueOf(Bonds[j].force));
		data41Label.setText("Charge");
		data42Label.setText(String.valueOf(Bonds[j].charge));
		dataPane.removeAll();
		addComponent(dataPane,data00Label,0,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data11Label,1,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data21Label,1,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data22Label,2,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data31Label,1,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data32Label,2,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data41Label,1,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data42Label,2,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		if(!atomSelected) header.remove(selectButton);
		dataPane.validate();
	    };
	    if(mouse.getSource()==hList){
		j=indexH[hList.getSelectedIndex()];
		data00Label.setText("Angle ");
		data11Label.setText("Type "+file_type[Angles[j].source]+": "+Angles[j].type1+"-"+Angles[j].type2+"-"+Angles[j].type3);
		data21Label.setText("Angle");
		data22Label.setText(String.valueOf(Angles[j].value));
		data31Label.setText("Force constant");
		data32Label.setText(String.valueOf(Angles[j].force));
		dataPane.removeAll();
		addComponent(dataPane,data00Label,0,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data11Label,1,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data21Label,1,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data22Label,2,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data31Label,1,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data32Label,2,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		if(!atomSelected) header.remove(selectButton);
		dataPane.validate();
	    };
	    if(mouse.getSource()==tList){
		j=indexT[tList.getSelectedIndex()];
		data00Label.setText("Torsion ");
		data11Label.setText("Type "+file_type[Torsions[j].source]+": "+Torsions[j].type1+"-"+Torsions[j].type2+"-"+Torsions[j].type3+"-"+Torsions[j].type4);
		data21Label.setText("Angle");
		data22Label.setText(String.valueOf(Torsions[j].value[0]));
		data31Label.setText("Force constant");
		data32Label.setText(String.valueOf(Torsions[j].force[0]));
		data41Label.setText("Multiplicity");
		data42Label.setText(String.valueOf(Torsions[j].multiplicity[0]));
		dataPane.removeAll();
		addComponent(dataPane,data00Label,0,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data11Label,1,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data21Label,1,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data22Label,2,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data31Label,1,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data32Label,2,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data41Label,1,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data42Label,2,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		if(Torsions[j].number>0){
		    data23Label.setText(", "+String.valueOf(Torsions[j].value[1]));
		    data33Label.setText(", "+String.valueOf(Torsions[j].force[1]));
		    data43Label.setText(", "+String.valueOf(Torsions[j].multiplicity[1]));
		    addComponent(dataPane,data23Label,3,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		    addComponent(dataPane,data33Label,3,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		    addComponent(dataPane,data43Label,3,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		};
		if(Torsions[j].number>1){
		    data24Label.setText(", "+String.valueOf(Torsions[j].value[2]));
		    data34Label.setText(", "+String.valueOf(Torsions[j].force[2]));
		    data44Label.setText(", "+String.valueOf(Torsions[j].multiplicity[2]));
		    addComponent(dataPane,data24Label,4,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		    addComponent(dataPane,data34Label,4,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		    addComponent(dataPane,data44Label,4,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		};
		if(Torsions[j].number>2){
		    data25Label.setText(", "+String.valueOf(Torsions[j].value[3]));
		    data35Label.setText(", "+String.valueOf(Torsions[j].force[3]));
		    data45Label.setText(", "+String.valueOf(Torsions[j].multiplicity[3]));
		    addComponent(dataPane,data25Label,5,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		    addComponent(dataPane,data35Label,5,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		    addComponent(dataPane,data45Label,5,3,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		};
		if(!atomSelected) header.remove(selectButton);
		dataPane.validate();
	    };
	    if(mouse.getSource()==iList){
		j=indexI[iList.getSelectedIndex()];
		data00Label.setText("Improper ");
		data11Label.setText("Type "+file_type[Impropers[j].source]+": "+Impropers[j].type1+"-"+Impropers[j].type2+"-"+Impropers[j].type3+"-"+Impropers[j].type4);
		data21Label.setText("Angle");
		data22Label.setText(String.valueOf(Impropers[j].value));
		data31Label.setText("Force constant");
		data32Label.setText(String.valueOf(Impropers[j].force));
		dataPane.removeAll();
		addComponent(dataPane,data00Label,0,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data11Label,1,0,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data21Label,1,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data22Label,2,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data31Label,1,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		addComponent(dataPane,data32Label,2,2,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.WEST);
		if(!atomSelected) header.remove(selectButton);
		dataPane.validate();
	    };
	};
    }
    
    public void mouseEntered(MouseEvent mouse){
    }

    public void mouseExited(MouseEvent mouse){
    }
  
}
