import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

class nwchem_Fragment extends JFrame implements ActionListener, ChangeListener, WindowListener, MouseListener {
    
    Font defaultFont;

    String ffield = "unknown";
    String dir_s = " ";
    String dir_x = " ";
    String dir_u = " ";
    String dir_t = " ";

    char[] file_type = { 's', 'x', 'u', 't'};

    JPanel header = new JPanel();

    DefaultListModel fragmentList = new DefaultListModel();
    JList frgList = new JList(fragmentList);
    JScrollPane frgPane = new JScrollPane(frgList);

    DefaultListModel atomList = new DefaultListModel();
    JList atmList = new JList(atomList);
    JScrollPane atmPane = new JScrollPane(atmList);

    File[] files;
    String[] dirs = {" ", " ", " ", " "};
    FragmentDefinition[] frgDef;
    
    JFileChooser chooser;
    ExtensionFilter frgFilter;

    JLabel frLabel = new JLabel(" From ");
    boolean frBool = true;

    JLabel help1 = new JLabel("FRG list: LEFT selects fragment, RIGHT toggles between fragments");
    JLabel help2 = new JLabel("ATM list: in SELECT mode: LEFT splits atoms");
    JLabel help3 = new JLabel("ATM list: in SELECT mode: CTRL-SHIFT-LEFT selects atom to join");
    JLabel help4 = new JLabel("ATM list: in SELECT mode: LEFT joins atom to previous selected atom");
    JLabel help5 = new JLabel("ATM list: in SELECT mode: SHIFT-LEFT moves atom one position up");
    JLabel help6 = new JLabel("ATM list: in SELECT mode: CTRL-LEFT moves atom one position down");
    JLabel help7 = new JLabel("ATM list: in ORDER mode: LEFT selects atoms in order");
    JLabel help8 = new JLabel("WRITE writes new fragment file to NEW.frg");
    JLabel help9 = new JLabel("ORDER/SELECT toggles between ORDER and SELECT mode");

    JLabel xxx = new JLabel("   ");
    JLabel yyy = new JLabel("   ");

    boolean frRead=false, toRead=false;

    boolean ordering = false;
    int iorder=0;

    int frIndex, toIndex;
    int atmNumber=0;
    int[][] id;
    int[] idf,idt;

    int selected = -1;
    
    Fragment ToFrg = new Fragment();
    Fragment FrFrg = new Fragment();

    JButton writeButton = new JButton("Write");
    JButton orderButton = new JButton("Order");
    
    public nwchem_Fragment(){
	
	super("Fragment");
	
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
	
	frgList.addMouseListener(this);
	frgList.setVisibleRowCount(25);

	atmList.addMouseListener(this);
	atmList.setVisibleRowCount(25);
	
	JButton doneButton = new JButton("Done");
	
	addComponent(header,doneButton,7,0,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.CENTER);
	doneButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){ 
		    setVisible(false); }});
	
	addComponent(header,writeButton,6,0,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.CENTER);
	writeButton.addActionListener(this);

	addComponent(header,orderButton,5,0,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.CENTER);
	orderButton.addActionListener(this);

	addComponent(header,xxx,8,0,1,1,10,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);

	frLabel.setBackground(Color.yellow);
	frLabel.setForeground(Color.blue);

	addComponent(header,frLabel,1,0,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.CENTER);

	addComponent(header,help1,5,1,10,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	addComponent(header,help2,5,2,10,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	addComponent(header,help3,5,3,10,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	addComponent(header,help4,5,4,10,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	addComponent(header,help5,5,5,10,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	addComponent(header,help6,5,6,10,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	addComponent(header,help7,5,7,10,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	addComponent(header,help8,5,8,10,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	addComponent(header,help9,5,9,10,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	addComponent(header,yyy,6,10,1,1,10,25,
		     GridBagConstraints.NONE,GridBagConstraints.NORTHWEST);
	
	try{
	    BufferedReader br = new BufferedReader(new FileReader("./.nwchemrc"));
	    String card;
	    while((card=br.readLine()) != null){
		if(card.startsWith("ffield")) {ffield=card.substring(7,12).trim();};
		if(card.startsWith(ffield+"_s")) {dirs[0]=card.substring(card.indexOf(" ")+1,card.length());};
		if(card.startsWith(ffield+"_x")) {dirs[1]=card.substring(card.indexOf(" ")+1,card.length());};
		if(card.startsWith(ffield+"_u")) {dirs[2]=card.substring(card.indexOf(" ")+1,card.length());};
		if(card.startsWith(ffield+"_t")) {dirs[3]=card.substring(card.indexOf(" ")+1,card.length());};
	    };
	    br.close();
	} catch(Exception e) {e.printStackTrace();};
	
	int frgNumber=0;
	
	for(int idir=0; idir<4; idir++){
	    if(dirs[idir]!=" "){
		File dir = new File(dirs[idir]);
		File[] files = dir.listFiles();
		String name;
		for(int i=0; i<files.length; i++) {
		    name=files[i].getName(); 
		    if(name.toLowerCase().endsWith(".frg")) frgNumber++;
		};
	    };
	};

	frgDef = new FragmentDefinition[frgNumber];
	frgNumber=0;
	for(int idir=0; idir<4; idir++){
	    if(dirs[idir]!=" "){
		File dir = new File(dirs[idir]);
		File[] files = dir.listFiles();
		String name;
		for(int i=0; i<files.length; i++) {
		    name=files[i].getName(); 
		    if(name.toLowerCase().endsWith(".frg")) {
			frgDef[frgNumber] = new FragmentDefinition();
			frgDef[frgNumber].Name=name;
			frgDef[frgNumber].Dir=dirs[idir];
			fragmentList.addElement(file_type[idir]+" "+name.substring(0,name.indexOf(".frg")));
			frgNumber++;
		    };
		};
	    };
	};
	addComponent(header,frgPane,0,1,1,12,1,1,GridBagConstraints.NONE,GridBagConstraints.CENTER);
	
	setLocation(25,225);	
	setSize(800,800);
	setVisible(true);
	
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
    
    void atomListUpdate(){
	//	header.remove(atmPane);
	atomList.removeAllElements();
	for(int i=0; i<atmNumber; i++){
	    if(id[i][0]>=0 && id[i][1]>=0) {
		atomList.addElement(FrFrg.atom[id[i][0]].Name+">"+ToFrg.atom[id[i][1]].Name);
	    } else if(id[i][0]>=0 && id[i][1]<0) {
		atomList.addElement(FrFrg.atom[id[i][0]].Name+"     ");
	    } else {
		atomList.addElement("     "+ToFrg.atom[id[i][1]].Name);
	    }
	};
	atmList.clearSelection();
	addComponent(header,atmPane,1,1,2,12,2,2,GridBagConstraints.NONE,GridBagConstraints.CENTER);
	header.validate();
    };

    public void actionPerformed(ActionEvent e){
	if(e.getSource()==writeButton){
	    idf = new int[FrFrg.numAtoms+1];
	    idt = new int[ToFrg.numAtoms+1];
	    for(int k=0; k<atmNumber; k++){
		if(id[k][0]>=0){idf[id[k][0]+1]=k;};
		if(id[k][1]>=0){idt[id[k][1]+1]=k;};
	    };
	    try{
		PrintfWriter frgFile = new PrintfWriter(new FileWriter("NEW.frg"));
		frgFile.println("# Merged Fragment File");
		frgFile.printf("%5d",atmNumber);
		// count new number of bonds
		int number=0;
		boolean found=false;
		int[][] ida = new int[FrFrg.numBonds+ToFrg.numBonds][2];
		for(int i=0; i<FrFrg.numBonds; i++){
		    ida[number][0]=idf[FrFrg.bond[i].atomi];
		    ida[number][1]=idf[FrFrg.bond[i].atomj];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1]) found=true;
			if(ida[j][0]==ida[number][1] && ida[j][1]==ida[number][0]) found=true;
		    };
		    if(!found) number++;
		};
		for(int i=0; i<ToFrg.numBonds; i++){
		    ida[number][0]=idt[ToFrg.bond[i].atomi];
		    ida[number][1]=idt[ToFrg.bond[i].atomj];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1]) found=true;
			if(ida[j][0]==ida[number][1] && ida[j][1]==ida[number][0]) found=true;
		    };
		    if(!found) number++;
		};
		frgFile.printf("%5d",number);
		// count new number of zmatrixs
		number=0;
		ida = new int[FrFrg.numZmatrix+ToFrg.numZmatrix][4];
		for(int i=0; i<FrFrg.numZmatrix; i++){
		    ida[number][0]=idf[FrFrg.zmatrix[i].atomi];
		    ida[number][1]=idf[FrFrg.zmatrix[i].atomj];
		    ida[number][2]=idf[FrFrg.zmatrix[i].atomk];
		    ida[number][3]=idf[FrFrg.zmatrix[i].atoml];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][3]) found=true;
		    };
		    if(!found) number++;
		};
		for(int i=0; i<ToFrg.numZmatrix; i++){
		    ida[number][0]=idt[ToFrg.zmatrix[i].atomi];
		    ida[number][1]=idt[ToFrg.zmatrix[i].atomj];
		    ida[number][2]=idt[ToFrg.zmatrix[i].atomk];
		    ida[number][3]=idt[ToFrg.zmatrix[i].atoml];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][3]) found=true;
		    };
		    if(!found) number++;
		};
		frgFile.printf("%41d",number);
		frgFile.println();
		int polgrp= 1;
		int chgrp = 1;
		double charge1 = 0.0;
		double charge2 = 0.0;
		double charge3 = 0.0;
		for(int k=0; k<atmNumber; k++){
		    frgFile.printf("%5d",k+1);
		    if(id[k][0]>=0){
			frgFile.print(FrFrg.atom[id[k][0]].Name+" ");
			frgFile.print(FrFrg.atom[id[k][0]].Type1+" ");
			frgFile.print(FrFrg.atom[id[k][0]].Type1+" ");
			if(id[k][1]>=0){
			    frgFile.print(ToFrg.atom[id[k][1]].Type1+" ");
			} else {
			    frgFile.print(FrFrg.atom[id[k][0]].Type1+"D");
			};
			frgFile.printf("%4d",chgrp);
			frgFile.printf("%4d",polgrp);
			// frgFile.printf("%4d",FrFrg.atom[id[k][0]].cgroup);
			// frgFile.printf("%4d",FrFrg.atom[id[k][0]].pgroup);
			frgFile.printf("%4d",FrFrg.atom[id[k][0]].link);
			frgFile.printf("%4d",FrFrg.atom[id[k][0]].type);
			frgFile.println();
			frgFile.printf("%12.6f",FrFrg.atom[id[k][0]].q1);
			frgFile.printf("%12.5E",FrFrg.atom[id[k][0]].p1);
			frgFile.printf("%12.6f",FrFrg.atom[id[k][0]].q1);
			frgFile.printf("%12.5E",FrFrg.atom[id[k][0]].p1);
			charge1=charge1+FrFrg.atom[id[k][0]].q1;
			if(id[k][1]>=0){
			    frgFile.printf("%12.6f",ToFrg.atom[id[k][1]].q1);
			    frgFile.printf("%12.5E",ToFrg.atom[id[k][1]].p1);
			    charge3=charge3+ToFrg.atom[id[k][1]].q1;
			} else {
			    frgFile.printf("%12.6f",0.0);
			    frgFile.printf("%12.5E",0.0);
			};
			frgFile.println();
		    } else {
			frgFile.print(ToFrg.atom[id[k][1]].Name+" ");
			frgFile.print(ToFrg.atom[id[k][1]].Type1+"D");
			frgFile.print(ToFrg.atom[id[k][1]].Type1+"D");
			frgFile.print(ToFrg.atom[id[k][1]].Type1+" ");
			frgFile.printf("%4d",chgrp);
			frgFile.printf("%4d",polgrp);
			// frgFile.printf("%4d",ToFrg.atom[id[k][1]].cgroup);
			// frgFile.printf("%4d",ToFrg.atom[id[k][1]].pgroup);
			frgFile.printf("%4d",ToFrg.atom[id[k][1]].link);
			frgFile.printf("%4d",ToFrg.atom[id[k][1]].type);
			frgFile.println();
			frgFile.printf("%12.6f",0.0);
			frgFile.printf("%12.5E",0.0);
			frgFile.printf("%12.6f",0.0);
			frgFile.printf("%12.5E",0.0);
			frgFile.printf("%12.6f",ToFrg.atom[id[k][1]].q1);
			frgFile.printf("%12.5E",ToFrg.atom[id[k][1]].p1);
			charge3=charge3+ToFrg.atom[id[k][1]].q1;
			frgFile.println();
		    };
		    if(charge1==0.0 && charge3==0.0) chgrp++;
		};
		number=0;
	        int jfound=-1;
		ida = new int[FrFrg.numBonds+ToFrg.numBonds][2];
		for(int i=0; i<FrFrg.numBonds; i++){
		    ida[number][0]=idf[FrFrg.bond[i].atomi];
		    ida[number][1]=idf[FrFrg.bond[i].atomj];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1]) found=true;
			if(ida[j][0]==ida[number][1] && ida[j][1]==ida[number][0]) found=true;
		    };
		    if(!found) {
			frgFile.printf("%5d",number+1);
			frgFile.printf("%5d",ida[number][0]+1);
			frgFile.printf("%5d",ida[number][1]+1);
			frgFile.printf("%5d",FrFrg.bond[i].type);
			jfound=-1;
			for(int j=0; j<ToFrg.numBonds; j++){
			    if(idt[ToFrg.bond[j].atomi]==ida[number][0] && idt[ToFrg.bond[j].atomj]==ida[number][1]) jfound=j;
			    if(idt[ToFrg.bond[j].atomj]==ida[number][0] && idt[ToFrg.bond[j].atomi]==ida[number][1]) jfound=j;
			};
			if(jfound>=0) {
			    frgFile.printf("%3d",FrFrg.bond[i].source);
			    frgFile.printf("%1d",FrFrg.bond[i].source);
			    frgFile.printf("%1d",ToFrg.bond[jfound].source); frgFile.println();
			    frgFile.printf("%12.6f",FrFrg.bond[i].bond1);
			    frgFile.printf("%12.5E",FrFrg.bond[i].force1);
			    frgFile.printf("%12.6f",FrFrg.bond[i].bond2);
			    frgFile.printf("%12.5E",FrFrg.bond[i].force2);
			    frgFile.printf("%12.6f",ToFrg.bond[jfound].bond1);
			    frgFile.printf("%12.5E",ToFrg.bond[jfound].force1);
			} else {
			    frgFile.printf("%3d",FrFrg.bond[i].source);
			    frgFile.printf("%1d",FrFrg.bond[i].source);
			    if(id[ida[number][0]][1]>=0 && id[ida[number][1]][1]>=0) {
				frgFile.printf("%1d",1);
			    } else { 
				frgFile.printf("%1d",FrFrg.bond[i].source); 
			    };
			    frgFile.println();
			    frgFile.printf("%12.6f",FrFrg.bond[i].bond1);
			    frgFile.printf("%12.5E",FrFrg.bond[i].force1);
			    frgFile.printf("%12.6f",FrFrg.bond[i].bond2);
			    frgFile.printf("%12.5E",FrFrg.bond[i].force2);
			    frgFile.printf("%12.6f",FrFrg.bond[i].bond3);
			    if(id[ida[number][0]][1]>=0 && id[ida[number][1]][1]>=0) {
				frgFile.printf("%12.5E",0.0);
			    } else { 
				frgFile.printf("%12.5E",FrFrg.bond[i].force3);
			    };
			};
			frgFile.println();
			number++;
		    };
		};
		for(int i=0; i<ToFrg.numBonds; i++){
		    ida[number][0]=idt[ToFrg.bond[i].atomi];
		    ida[number][1]=idt[ToFrg.bond[i].atomj];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1]) found=true;
			if(ida[j][0]==ida[number][1] && ida[j][1]==ida[number][0]) found=true;
		    };
		    if(!found) {
			frgFile.printf("%5d",number+1);
			frgFile.printf("%5d",ida[number][0]+1);
			frgFile.printf("%5d",ida[number][1]+1);
			frgFile.printf("%5d",ToFrg.bond[i].type);
			if(id[ida[number][0]][0]>=0 && id[ida[number][1]][0]>=0) {
			    frgFile.printf("%3d",1); frgFile.printf("%1d",1);
			    frgFile.printf("%1d",ToFrg.bond[i].source); frgFile.println();
			    frgFile.printf("%12.6f",ToFrg.bond[i].bond1);
			    frgFile.printf("%12.5E",0.0);
			    frgFile.printf("%12.6f",ToFrg.bond[i].bond2);
			    frgFile.printf("%12.5E",0.0);
			} else {
			    frgFile.printf("%3d",ToFrg.bond[i].source); frgFile.printf("%1d",ToFrg.bond[i].source);
			    frgFile.printf("%1d",ToFrg.bond[i].source); frgFile.println();
			    frgFile.printf("%12.6f",ToFrg.bond[i].bond1);
			    frgFile.printf("%12.5E",ToFrg.bond[i].force1);
			    frgFile.printf("%12.6f",ToFrg.bond[i].bond2);
			    frgFile.printf("%12.5E",ToFrg.bond[i].force2);
			};
			frgFile.printf("%12.6f",ToFrg.bond[i].bond3);
			frgFile.printf("%12.5E",ToFrg.bond[i].force3);
			frgFile.println();
			number++;
		    };
		};
		number=0;
	        jfound=-1;
		ida = new int[FrFrg.numZmatrix+ToFrg.numZmatrix][4];
		for(int i=0; i<FrFrg.numZmatrix; i++){
		    ida[number][0]=idf[FrFrg.zmatrix[i].atomi];
		    ida[number][1]=idf[FrFrg.zmatrix[i].atomj];
		    ida[number][2]=idf[FrFrg.zmatrix[i].atomk];
		    ida[number][3]=idf[FrFrg.zmatrix[i].atoml];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2]
			   && ida[j][3]==ida[number][3]) found=true;
		    };
		    if(!found) {
			frgFile.printf("%5d",number+1);
			frgFile.printf("%5d",ida[number][0]+1);
			frgFile.printf("%5d",ida[number][1]+1);
			frgFile.printf("%5d",ida[number][2]+1);
			frgFile.printf("%5d",ida[number][3]+1);
			frgFile.printf("%12.6f",FrFrg.zmatrix[i].bond);
			frgFile.printf("%12.6f",FrFrg.zmatrix[i].angle);
			frgFile.printf("%12.6f",FrFrg.zmatrix[i].torsion);
			frgFile.println();
			number++;
		    };
		};
		for(int i=0; i<ToFrg.numZmatrix; i++){
		    ida[number][0]=idt[ToFrg.zmatrix[i].atomi];
		    ida[number][1]=idt[ToFrg.zmatrix[i].atomj];
		    ida[number][2]=idt[ToFrg.zmatrix[i].atomk];
		    ida[number][3]=idt[ToFrg.zmatrix[i].atoml];
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2]
			   && ida[j][3]==ida[number][3]) found=true;
			if(ida[j][0]==ida[number][3] && ida[j][1]==ida[number][2] && ida[j][2]==ida[number][1]
			   && ida[j][3]==ida[number][0]) found=true;
		    };
		    if(!found) {
			frgFile.printf("%5d",number+1);
			frgFile.printf("%5d",ida[number][0]+1);
			frgFile.printf("%5d",ida[number][1]+1);
			frgFile.printf("%5d",ida[number][2]+1);
			frgFile.printf("%5d",ida[number][3]+1);
			frgFile.printf("%12.6f",ToFrg.zmatrix[i].bond);
			frgFile.printf("%12.6f",ToFrg.zmatrix[i].angle);
			frgFile.printf("%12.6f",ToFrg.zmatrix[i].torsion);
			frgFile.println();
			number++;
		    };
		};
		frgFile.close();
	    } catch (Exception ee) { ee.printStackTrace(); };
	};
	if(e.getSource()==orderButton){
	    if(ordering) {
		orderButton.setText("Order");
		ordering=false;
		iorder=0;
	    } else {
		orderButton.setText("Select");
		ordering=true;
		iorder=0;
	    };
	};
    };
    
    public void stateChanged(ChangeEvent e){}
    
    public void windowClosing(WindowEvent event){}
    
    public void windowClosed(WindowEvent event){}
    
    public void windowDeiconified(WindowEvent event){}
    
    public void windowIconified(WindowEvent event){}
    
    public void windowActivated(WindowEvent event){}
    
    public void windowDeactivated(WindowEvent e){}
    
    public void windowOpened(WindowEvent event){}
    
    public void mouseClicked(MouseEvent mouse){}
    
    public void mousePressed(MouseEvent mouse){}
    
    public void mouseReleased(MouseEvent mouse){
	int j;
	String fileName;
	if(mouse.getModifiers()==(MouseEvent.BUTTON1_MASK+MouseEvent.CTRL_MASK+MouseEvent.SHIFT_MASK)){
	    if(mouse.getSource()==atmList){ 
		selected=atmList.getSelectedIndex();
	    } else {
		selected = -1;
	    };
	};
	if(mouse.getModifiers()==(MouseEvent.BUTTON1_MASK+MouseEvent.SHIFT_MASK)){
	    if(mouse.getSource()==atmList){
		j=atmList.getSelectedIndex();
		if(j>0){
		    int tempo = id[j][0]; id[j][0]=id[j-1][0]; id[j-1][0]=tempo;
		    tempo = id[j][1]; id[j][1]=id[j-1][1]; id[j-1][1]=tempo;
		};
		atomListUpdate();
		atmList.clearSelection();
	    };
	    selected = -1;
	};
	if(mouse.getModifiers()==(MouseEvent.BUTTON1_MASK+MouseEvent.CTRL_MASK)){
	    if(mouse.getSource()==atmList){
		j=atmList.getSelectedIndex();
		if(j<(atmNumber-1)){
		    int tempo = id[j][0]; id[j][0]=id[j+1][0]; id[j+1][0]=tempo;
		    tempo = id[j][1]; id[j][1]=id[j+1][1]; id[j+1][1]=tempo;
		};
		atomListUpdate();
		atmList.clearSelection();
	    };
	    selected = -1;
	};
	if(mouse.getModifiers()==MouseEvent.BUTTON1_MASK){
	    if(mouse.getSource()==frgList){
		j=frgList.getSelectedIndex();
		fileName=frgDef[j].Dir+frgDef[j].Name;
		FrFrg.read(fileName);
		frLabel.setText(" "+frgDef[j].Name.substring(0,frgDef[j].Name.indexOf(".frg")));
		frRead=true;
		if(frRead && !toRead){
		    id = new int[FrFrg.numAtoms][2];
		    atmNumber=0;
		    for(int i=0; i<FrFrg.numAtoms; i++){
			id[i][0]=i; id[i][1]=-1; atmNumber++;
		    };
		};
		if(frRead && toRead){
		    boolean[] frFnd = new boolean[FrFrg.numAtoms];
		    for(int i=0; i<FrFrg.numAtoms; i++){frFnd[i]=false;};
		    boolean[] toFnd = new boolean[ToFrg.numAtoms];
		    for(int i=0; i<ToFrg.numAtoms; i++){toFnd[i]=false;};
		    int num = FrFrg.numAtoms+ToFrg.numAtoms;
		    id = new int[num][2];
		    atmNumber=0;
		    for(int i=0; i<num; i++){ id[i][0]=-1; id[i][1]=-1;};
		    for(int i=0; i<FrFrg.numAtoms; i++){
			if(!frFnd[i]){
			    for(int k=0; k<ToFrg.numAtoms; k++){
			    	if(!toFnd[k]){
				    if(FrFrg.atom[i].Name.equals(ToFrg.atom[k].Name)){
					id[atmNumber][0]=i; id[atmNumber][1]=k; atmNumber++; frFnd[i]=true; toFnd[k]=true;
				    };
			    	};
			    };
			};
			if(!frFnd[i]){
			    id[atmNumber][0]=i; atmNumber++; frFnd[i]=true;
			};
		    };
		    for(int k=0; k<ToFrg.numAtoms; k++){
			if(!toFnd[k]){
			    id[atmNumber][1]=k; atmNumber++; toFnd[k]=true;
			};
		    };
		};
		selected = -1;
		atomListUpdate();
	    };
	    if(mouse.getSource()==atmList){
		j=atmList.getSelectedIndex();
		if(ordering){
		    if(iorder<atmNumber && j!=iorder){
			    int tempo = id[iorder][0]; id[iorder][0]=id[j][0]; id[j][0]=tempo;
			    tempo = id[iorder][1]; id[iorder][1]=id[j][1]; id[j][1]=tempo;
			    atomListUpdate();
		    };
		    iorder++;
		    if(iorder>=atmNumber){
			ordering=false;
			iorder=0;
			orderButton.setText("Order");
		    };
		} else {
		    if(id[j][0]>=0 && id[j][1]>=0){
			for(int k=atmNumber; k>j; k--){id[k][0]=id[k-1][0]; id[k][1]=id[k-1][1];};
			if(ToFrg.atom[id[j+1][1]].Name.charAt(4)==' '){
			    ToFrg.atom[id[j+1][1]].Name=ToFrg.atom[id[j+1][1]].Name.substring(0,4)+"t";
			};
			id[j][1]=-1; id[j+1][0]=-1;  atmNumber++;
		    };
		    if(selected>=0 && selected!=j){
			if(id[j][0]==-1 && id[j][1]>=0 && id[selected][0]>=0 && id[selected][1]==-1){
			    id[selected][1]=id[j][1];
			    if(FrFrg.atom[id[selected][0]].Name.substring(0,4).equals(ToFrg.atom[id[selected][1]].Name.substring(0,4)) &&
			       ToFrg.atom[id[selected][1]].Name.charAt(4)=='t') ToFrg.atom[id[selected][1]].Name=ToFrg.atom[id[selected][1]].Name.substring(0,4)+" ";
			    atmNumber--;
			    for(int k=j; k<atmNumber; k++){
				id[k][0]=id[k+1][0];
				id[k][1]=id[k+1][1];
			    };
			};
			if(id[j][0]>=0 && id[j][1]==-1 && id[selected][0]==-1 && id[selected][1]>=0){
			    id[j][1]=id[selected][1];
			    if(FrFrg.atom[id[j][0]].Name.substring(0,4).equals(ToFrg.atom[id[j][1]].Name.substring(0,4)) &&
			       ToFrg.atom[id[j][1]].Name.charAt(4)=='t') ToFrg.atom[id[j][1]].Name=ToFrg.atom[id[j][1]].Name.substring(0,4)+" ";
			    atmNumber--;
			    for(int k=selected; k<atmNumber; k++){
				id[k][0]=id[k+1][0];
				id[k][1]=id[k+1][1];
			    };
			};
		    };
		};
		selected= -1;
		atomListUpdate();
	    };
	};
	if(mouse.getModifiers()==MouseEvent.BUTTON3_MASK){
	    if(mouse.getSource()==frgList){
		if(frBool) {
		    frLabel.setBackground(Color.lightGray);
		    frLabel.setForeground(Color.darkGray);
		} else {
		    frLabel.setBackground(Color.yellow);
		    frLabel.setForeground(Color.blue);
		};
	      	frBool=!frBool;
	    };
	};
    };
    
    public void mouseEntered(MouseEvent mouse){}

    public void mouseExited(MouseEvent mouse){}
  
}


