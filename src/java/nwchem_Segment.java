import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

class nwchem_Segment extends JFrame implements ActionListener, ChangeListener, WindowListener, MouseListener {
    
    Font defaultFont;

    String ffield = "unknown";
    String dir_s = " ";
    String dir_x = " ";
    String dir_q = " ";
    String dir_u = " ";
    String dir_t = " ";
    String dir_c = " ";

    char[] file_type = { 's', 'x', 'q', 'u', 't', 'c'};

    JPanel header = new JPanel();

    DefaultListModel segmentList = new DefaultListModel();
    JList sgmList = new JList(segmentList);
    JScrollPane sgmPane = new JScrollPane(sgmList);

    DefaultListModel atomList = new DefaultListModel();
    JList atmList = new JList(atomList);
    JScrollPane atmPane = new JScrollPane(atmList);

    File[] files;
    String[] dirs = {" ", " ", " ", " ", " ", "./"};
    SegmentDefinition[] sgmDef;
    
    JFileChooser chooser;
    ExtensionFilter sgmFilter;

    JLabel frLabel = new JLabel(" From ");
    JLabel toLabel = new JLabel("  To  ");
    boolean frBool = true;

    JLabel help1 = new JLabel("SGM list: LEFT selects segment, RIGHT toggles between segments");
    JLabel help2 = new JLabel("ATM list: in SELECT mode: LEFT splits atom pairs");
    JLabel help3 = new JLabel("ATM list: in SELECT mode: LEFT selects first atom for join");
    JLabel help4 = new JLabel("ATM list: in SELECT mode: LEFT joins second atom to previous selected first atom");
    JLabel help5 = new JLabel("ATM list: in SELECT mode: CTRL-LEFT moves atom one position up");
    JLabel help6 = new JLabel("ATM list: in SELECT mode: CTRL-RIGHT moves atom one position down");
    JLabel help7 = new JLabel("ATM list: in ORDER mode: LEFT selects atoms in order");
    JLabel help8 = new JLabel("WRITE writes new segment file to NEW.sgm");
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
    
    Segment ToSgm = new Segment();
    Segment FrSgm = new Segment();

    JButton writeButton = new JButton("Write");
    JButton orderButton = new JButton("Order");

    String frName = new String(" ");
    String toName = new String(" ");
    
    public nwchem_Segment(){
	
	super("Segment");
	
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
	
	sgmList.addMouseListener(this);
	sgmList.setVisibleRowCount(25);

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
	toLabel.setBackground(Color.lightGray);
	frLabel.setForeground(Color.blue);
	toLabel.setForeground(Color.darkGray);

	addComponent(header,frLabel,1,0,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.CENTER);
	addComponent(header,toLabel,2,0,1,1,1,1,
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
		if(card.startsWith(ffield+"_q")) {dirs[2]=card.substring(card.indexOf(" ")+1,card.length());};
		if(card.startsWith(ffield+"_u")) {dirs[3]=card.substring(card.indexOf(" ")+1,card.length());};
		if(card.startsWith(ffield+"_t")) {dirs[4]=card.substring(card.indexOf(" ")+1,card.length());};
		if(card.startsWith(ffield+"_c")) {dirs[5]=card.substring(card.indexOf(" ")+1,card.length());};
	    };
	    br.close();
	} catch(Exception e) {e.printStackTrace();};
	
	int sgmNumber=0;
	
	for(int idir=0; idir<6; idir++){
	    if(dirs[idir]!=" "){
		File dir = new File(dirs[idir]);
		File[] files = dir.listFiles();
		String name;
		for(int i=0; i<files.length; i++) {
		    name=files[i].getName(); 
		    if(name.toLowerCase().endsWith(".sgm")) sgmNumber++;
		};
	    };
	};

	sgmDef = new SegmentDefinition[sgmNumber];
	sgmNumber=0;
	for(int idir=0; idir<6; idir++){
	    if(dirs[idir]!=" "){
		File dir = new File(dirs[idir]);
		File[] files = dir.listFiles();
		String name;
		for(int i=0; i<files.length; i++) {
		    name=files[i].getName(); 
		    if(name.toLowerCase().endsWith(".sgm")) {
			sgmDef[sgmNumber] = new SegmentDefinition();
			sgmDef[sgmNumber].Name=name;
			sgmDef[sgmNumber].Dir=dirs[idir];
			segmentList.addElement(file_type[idir]+" "+name.substring(0,name.indexOf(".sgm")));
			sgmNumber++;
		    };
		};
	    };
	};
	addComponent(header,sgmPane,0,1,1,12,1,1,GridBagConstraints.NONE,GridBagConstraints.CENTER);
	
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
		atomList.addElement(FrSgm.atom[id[i][0]].Name+">"+ToSgm.atom[id[i][1]].Name);
	    } else if(id[i][0]>=0 && id[i][1]<0) {
		atomList.addElement(FrSgm.atom[id[i][0]].Name+"     ");
	    } else {
		atomList.addElement("     "+ToSgm.atom[id[i][1]].Name);
	    }
	};
	atmList.clearSelection();
	addComponent(header,atmPane,1,1,2,12,2,2,GridBagConstraints.NONE,GridBagConstraints.CENTER);
	header.validate();
    };

    public void actionPerformed(ActionEvent e){
	if(e.getSource()==writeButton){
	    idf = new int[FrSgm.numAtoms+1];
	    idt = new int[ToSgm.numAtoms+1];
	    for(int k=0; k<atmNumber; k++){
		if(id[k][0]>=0){idf[id[k][0]+1]=k;};
		if(id[k][1]>=0){idt[id[k][1]+1]=k;};
	    };
	    try{
		PrintfWriter sgmFile = new PrintfWriter(new FileWriter("TEMP.sgm"));
		sgmFile.println("# Merged Segment File for "+frName+" to "+toName);
		sgmFile.println("    4.600000");
		sgmFile.printf("%5d",atmNumber);
		// count new number of bonds
		int number=0;
		boolean found=false;
		int[][] ida = new int[FrSgm.numBonds+ToSgm.numBonds][2];
		for(int i=0; i<FrSgm.numBonds; i++){
		    ida[number][0]=idf[FrSgm.bond[i].atomi];
		    ida[number][1]=idf[FrSgm.bond[i].atomj];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1]) found=true;
			if(ida[j][0]==ida[number][1] && ida[j][1]==ida[number][0]) found=true;
		    };
		    if(!found) number++;
		};
		for(int i=0; i<ToSgm.numBonds; i++){
		    ida[number][0]=idt[ToSgm.bond[i].atomi];
		    ida[number][1]=idt[ToSgm.bond[i].atomj];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1]) found=true;
			if(ida[j][0]==ida[number][1] && ida[j][1]==ida[number][0]) found=true;
		    };
		    if(!found) number++;
		};
		sgmFile.printf("%5d",number);
		// count new number of angles
		number=0;
		ida = new int[FrSgm.numAngles+ToSgm.numAngles][3];
		for(int i=0; i<FrSgm.numAngles; i++){
		    ida[number][0]=idf[FrSgm.angle[i].atomi];
		    ida[number][1]=idf[FrSgm.angle[i].atomj];
		    ida[number][2]=idf[FrSgm.angle[i].atomk];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2]) found=true;
			if(ida[j][0]==ida[number][2] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][0]) found=true;
		    };
		    if(!found) number++;
		};
		for(int i=0; i<ToSgm.numAngles; i++){
		    ida[number][0]=idt[ToSgm.angle[i].atomi];
		    ida[number][1]=idt[ToSgm.angle[i].atomj];
		    ida[number][2]=idt[ToSgm.angle[i].atomk];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2]) found=true;
			if(ida[j][0]==ida[number][2] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][0]) found=true;
		    };
		    if(!found) number++;
		};
		sgmFile.printf("%5d",number);
		// count new number of torsions
		number=0;
		ida = new int[FrSgm.numTorsions+ToSgm.numTorsions][4];
		for(int i=0; i<FrSgm.numTorsions; i++){
		    ida[number][0]=idf[FrSgm.torsion[i].atomi];
		    ida[number][1]=idf[FrSgm.torsion[i].atomj];
		    ida[number][2]=idf[FrSgm.torsion[i].atomk];
		    ida[number][3]=idf[FrSgm.torsion[i].atoml];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][3]) found=true;
			if(ida[j][0]==ida[number][3] && ida[j][1]==ida[number][2] && ida[j][2]==ida[number][1] && ida[j][3]==ida[number][0]) found=true;
		    };
		    if(!found) number++;
		};
		for(int i=0; i<ToSgm.numTorsions; i++){
		    ida[number][0]=idt[ToSgm.torsion[i].atomi];
		    ida[number][1]=idt[ToSgm.torsion[i].atomj];
		    ida[number][2]=idt[ToSgm.torsion[i].atomk];
		    ida[number][3]=idt[ToSgm.torsion[i].atoml];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][3]) found=true;
			if(ida[j][0]==ida[number][3] && ida[j][1]==ida[number][2] && ida[j][2]==ida[number][1] && ida[j][3]==ida[number][0]) found=true;
		    };
		    if(!found) number++;
		};
		sgmFile.printf("%5d",number);
		// count new number of impropers
		number=0;
		ida = new int[FrSgm.numImpropers+ToSgm.numImpropers][4];
		for(int i=0; i<FrSgm.numImpropers; i++){
		    ida[number][0]=idf[FrSgm.improper[i].atomi];
		    ida[number][1]=idf[FrSgm.improper[i].atomj];
		    ida[number][2]=idf[FrSgm.improper[i].atomk];
		    ida[number][3]=idf[FrSgm.improper[i].atoml];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][3]) found=true;
			if(ida[j][0]==ida[number][1] && ida[j][1]==ida[number][0] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][3]) found=true;
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][3] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][1]) found=true;
			if(ida[j][0]==ida[number][3] && ida[j][1]==ida[number][0] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][1]) found=true;
			if(ida[j][0]==ida[number][1] && ida[j][1]==ida[number][3] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][0]) found=true;
			if(ida[j][0]==ida[number][3] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][0]) found=true;
		    };
		    if(!found) number++;
		};
		for(int i=0; i<ToSgm.numImpropers; i++){
		    ida[number][0]=idt[ToSgm.improper[i].atomi];
		    ida[number][1]=idt[ToSgm.improper[i].atomj];
		    ida[number][2]=idt[ToSgm.improper[i].atomk];
		    ida[number][3]=idt[ToSgm.improper[i].atoml];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][3]) found=true;
			if(ida[j][0]==ida[number][1] && ida[j][1]==ida[number][0] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][3]) found=true;
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][3] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][1]) found=true;
			if(ida[j][0]==ida[number][3] && ida[j][1]==ida[number][0] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][1]) found=true;
			if(ida[j][0]==ida[number][1] && ida[j][1]==ida[number][3] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][0]) found=true;
			if(ida[j][0]==ida[number][3] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][0]) found=true;
		    };
		    if(!found) number++;
		};
		sgmFile.printf("%5d",number);
		// count new number of zmatrixs
		number=0;
		ida = new int[FrSgm.numZmatrix+ToSgm.numZmatrix][4];
		for(int i=0; i<FrSgm.numZmatrix; i++){
		    ida[number][0]=idf[FrSgm.zmatrix[i].atomi];
		    ida[number][1]=idf[FrSgm.zmatrix[i].atomj];
		    ida[number][2]=idf[FrSgm.zmatrix[i].atomk];
		    ida[number][3]=idf[FrSgm.zmatrix[i].atoml];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][3]) found=true;
		    };
		    if(!found) number++;
		};
		for(int i=0; i<ToSgm.numZmatrix; i++){
		    ida[number][0]=idt[ToSgm.zmatrix[i].atomi];
		    ida[number][1]=idt[ToSgm.zmatrix[i].atomj];
		    ida[number][2]=idt[ToSgm.zmatrix[i].atomk];
		    ida[number][3]=idt[ToSgm.zmatrix[i].atoml];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][3]) found=true;
		    };
		    if(!found) number++;
		};
		sgmFile.printf("%5d",number);
		sgmFile.printf("%5d",3);
		sgmFile.printf("%5d",1);
		sgmFile.println();
		sgmFile.println("    0.000000");
		int polgrp= 1;
		int chgrp = 1;
		double charge1 = 0.0;
		double charge2 = 0.0;
		double charge3 = 0.0;
		for(int k=0; k<atmNumber; k++){
		    sgmFile.printf("%5d",k+1);
		    if(id[k][0]>=0){
			sgmFile.print(FrSgm.atom[id[k][0]].Name+" ");
			sgmFile.printf("%5d",FrSgm.atom[id[k][0]].link);
			sgmFile.printf("%5d",FrSgm.atom[id[k][0]].type);
			sgmFile.printf("%5d",0);
			sgmFile.printf("%5d",chgrp);
			sgmFile.printf("%5d",polgrp);
			sgmFile.println();
			sgmFile.print("     "+FrSgm.atom[id[k][0]].Type1);
			sgmFile.printf("%12.6f",FrSgm.atom[id[k][0]].q1);
			sgmFile.printf("%12.5E",FrSgm.atom[id[k][0]].p1);
			sgmFile.println();
			sgmFile.print("     "+FrSgm.atom[id[k][0]].Type1);
			sgmFile.printf("%12.6f",FrSgm.atom[id[k][0]].q1);
			sgmFile.printf("%12.5E",FrSgm.atom[id[k][0]].p1);
			sgmFile.println();
			if(id[k][1]>=0){
			    sgmFile.print("     "+ToSgm.atom[id[k][1]].Type1);
			} else {
			    sgmFile.print("     "+FrSgm.atom[id[k][0]].Type1.substring(0,5)+"D");
			};
			charge1=charge1+FrSgm.atom[id[k][0]].q1;
			if(id[k][1]>=0){
			    sgmFile.printf("%12.6f",ToSgm.atom[id[k][1]].q1);
			    sgmFile.printf("%12.5E",ToSgm.atom[id[k][1]].p1);
			    charge3=charge3+ToSgm.atom[id[k][1]].q1;
			} else {
			    sgmFile.printf("%12.6f",0.0);
			    sgmFile.printf("%12.5E",0.0);
			};
			sgmFile.println();
		    } else {
			sgmFile.print(ToSgm.atom[id[k][1]].Name+" ");
			sgmFile.printf("%5d",ToSgm.atom[id[k][1]].link);
			sgmFile.printf("%5d",ToSgm.atom[id[k][1]].type);
			sgmFile.printf("%5d",0);
			sgmFile.printf("%5d",chgrp);
			sgmFile.printf("%5d",polgrp);
			sgmFile.println();
			sgmFile.print("     "+ToSgm.atom[id[k][1]].Type1.substring(0,5)+"D");
			sgmFile.printf("%12.6f",0.0);
			sgmFile.printf("%12.5E",0.0);
			sgmFile.println();
			sgmFile.print("     "+ToSgm.atom[id[k][1]].Type1.substring(0,5)+"D");
			sgmFile.printf("%12.6f",0.0);
			sgmFile.printf("%12.5E",0.0);
			sgmFile.println();
			sgmFile.print("     "+ToSgm.atom[id[k][1]].Type1);
			sgmFile.printf("%12.6f",ToSgm.atom[id[k][1]].q1);
			sgmFile.printf("%12.5E",ToSgm.atom[id[k][1]].p1);
			charge3=charge3+ToSgm.atom[id[k][1]].q1;
			sgmFile.println();
		    };
		    if(charge1==0.0 && charge3==0.0) chgrp++;
		};
		number=0;
	        int jfound=-1;
		ida = new int[FrSgm.numBonds+ToSgm.numBonds][2];
		for(int i=0; i<FrSgm.numBonds; i++){
		    ida[number][0]=idf[FrSgm.bond[i].atomi];
		    ida[number][1]=idf[FrSgm.bond[i].atomj];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1]) found=true;
			if(ida[j][0]==ida[number][1] && ida[j][1]==ida[number][0]) found=true;
		    };
		    if(!found) {
			sgmFile.printf("%5d",number+1);
			sgmFile.printf("%5d",ida[number][0]+1);
			sgmFile.printf("%5d",ida[number][1]+1);
			sgmFile.printf("%5d",FrSgm.bond[i].type);
			jfound=-1;
			for(int j=0; j<ToSgm.numBonds; j++){
			    if(idt[ToSgm.bond[j].atomi]==ida[number][0] && idt[ToSgm.bond[j].atomj]==ida[number][1]) jfound=j;
			    if(idt[ToSgm.bond[j].atomj]==ida[number][0] && idt[ToSgm.bond[j].atomi]==ida[number][1]) jfound=j;
			};
			if(jfound>=0) {
			    sgmFile.printf("%3d",FrSgm.bond[i].source);
			    sgmFile.printf("%1d",FrSgm.bond[i].source);
			    sgmFile.printf("%1d",ToSgm.bond[jfound].source); sgmFile.println();
			    sgmFile.printf("%12.6f",FrSgm.bond[i].bond1);
			    sgmFile.printf("%12.5E",FrSgm.bond[i].force1);
			    sgmFile.println();
			    sgmFile.printf("%12.6f",FrSgm.bond[i].bond2);
			    sgmFile.printf("%12.5E",FrSgm.bond[i].force2);
			    sgmFile.println();
			    sgmFile.printf("%12.6f",ToSgm.bond[jfound].bond1);
			    sgmFile.printf("%12.5E",ToSgm.bond[jfound].force1);
			} else {
			    sgmFile.printf("%3d",FrSgm.bond[i].source);
			    sgmFile.printf("%1d",FrSgm.bond[i].source);
			    if(id[ida[number][0]][1]>=0 && id[ida[number][1]][1]>=0) {
				sgmFile.printf("%1d",1);
			    } else { 
				sgmFile.printf("%1d",FrSgm.bond[i].source); 
			    };
			    sgmFile.println();
			    sgmFile.printf("%12.6f",FrSgm.bond[i].bond1);
			    sgmFile.printf("%12.5E",FrSgm.bond[i].force1);
			    sgmFile.println();
			    sgmFile.printf("%12.6f",FrSgm.bond[i].bond2);
			    sgmFile.printf("%12.5E",FrSgm.bond[i].force2);
			    sgmFile.println();
			    sgmFile.printf("%12.6f",FrSgm.bond[i].bond3);
			    if(id[ida[number][0]][1]>=0 && id[ida[number][1]][1]>=0) {
				sgmFile.printf("%12.5E",0.0);
			    } else { 
				sgmFile.printf("%12.5E",FrSgm.bond[i].force3);
			    };
			};
			sgmFile.println();
			number++;
		    };
		};
		for(int i=0; i<ToSgm.numBonds; i++){
		    ida[number][0]=idt[ToSgm.bond[i].atomi];
		    ida[number][1]=idt[ToSgm.bond[i].atomj];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1]) found=true;
			if(ida[j][0]==ida[number][1] && ida[j][1]==ida[number][0]) found=true;
		    };
		    if(!found) {
			sgmFile.printf("%5d",number+1);
			sgmFile.printf("%5d",ida[number][0]+1);
			sgmFile.printf("%5d",ida[number][1]+1);
			sgmFile.printf("%5d",ToSgm.bond[i].type);
			if(id[ida[number][0]][0]>=0 && id[ida[number][1]][0]>=0) {
			    sgmFile.printf("%3d",1); sgmFile.printf("%1d",1);
			    sgmFile.printf("%1d",ToSgm.bond[i].source); sgmFile.println();
			    sgmFile.printf("%12.6f",ToSgm.bond[i].bond1);
			    sgmFile.printf("%12.5E",0.0);
			    sgmFile.println();
			    sgmFile.printf("%12.6f",ToSgm.bond[i].bond2);
			    sgmFile.printf("%12.5E",0.0);
			} else {
			    sgmFile.printf("%3d",ToSgm.bond[i].source); sgmFile.printf("%1d",ToSgm.bond[i].source);
			    sgmFile.printf("%1d",ToSgm.bond[i].source); sgmFile.println();
			    sgmFile.printf("%12.6f",ToSgm.bond[i].bond1);
			    sgmFile.printf("%12.5E",ToSgm.bond[i].force1);
			    sgmFile.println();
			    sgmFile.printf("%12.6f",ToSgm.bond[i].bond2);
			    sgmFile.printf("%12.5E",ToSgm.bond[i].force2);
			};
			sgmFile.println();
			sgmFile.printf("%12.6f",ToSgm.bond[i].bond3);
			sgmFile.printf("%12.5E",ToSgm.bond[i].force3);
			sgmFile.println();
			number++;
		    };
		};
		number=0;
	        jfound=-1;
		ida = new int[FrSgm.numAngles+ToSgm.numAngles][3];
		for(int i=0; i<FrSgm.numAngles; i++){
		    ida[number][0]=idf[FrSgm.angle[i].atomi];
		    ida[number][1]=idf[FrSgm.angle[i].atomj];
		    ida[number][2]=idf[FrSgm.angle[i].atomk];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2]) found=true;
			if(ida[j][0]==ida[number][2] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][0]) found=true;
		    };
		    if(!found) {
			sgmFile.printf("%5d",number+1);
			sgmFile.printf("%5d",ida[number][0]+1);
			sgmFile.printf("%5d",ida[number][1]+1);
			sgmFile.printf("%5d",ida[number][2]+1);
			sgmFile.printf("%5d",FrSgm.angle[i].type);
			jfound=-1;
			for(int j=0; j<ToSgm.numAngles; j++){
			    if(idt[ToSgm.angle[j].atomi]==ida[number][0] && idt[ToSgm.angle[j].atomj]==ida[number][1]
			       && idt[ToSgm.angle[j].atomk]==ida[number][2]) jfound=j;
			    if(idt[ToSgm.angle[j].atomk]==ida[number][0] && idt[ToSgm.angle[j].atomj]==ida[number][1]
			       && idt[ToSgm.angle[j].atomi]==ida[number][2]) jfound=j;
			};
			if(jfound>=0) {
			    sgmFile.printf("%3d",FrSgm.angle[i].source);
			    sgmFile.printf("%1d",FrSgm.angle[i].source);
			    sgmFile.printf("%1d",ToSgm.angle[jfound].source); sgmFile.println();
			    sgmFile.printf("%10.6f",FrSgm.angle[i].angle1);
			    sgmFile.printf("%12.5E",FrSgm.angle[i].force1);
			    sgmFile.println();
			    sgmFile.printf("%10.6f",FrSgm.angle[i].angle2);
			    sgmFile.printf("%12.5E",FrSgm.angle[i].force2);
			    sgmFile.println();
			    sgmFile.printf("%10.6f",ToSgm.angle[jfound].angle1);
			    sgmFile.printf("%12.5E",ToSgm.angle[jfound].force1);
			} else {
			    sgmFile.printf("%3d",FrSgm.angle[i].source);
			    sgmFile.printf("%1d",FrSgm.angle[i].source);
			    if(id[ida[number][0]][1]>=0 && id[ida[number][1]][1]>=0 && id[ida[number][2]][1]>=0) {
				sgmFile.printf("%1d",1);
			    } else { 
				sgmFile.printf("%1d",FrSgm.angle[i].source); 
			    };
			    sgmFile.println();
			    sgmFile.printf("%10.6f",FrSgm.angle[i].angle1);
			    sgmFile.printf("%12.5E",FrSgm.angle[i].force1);
			    sgmFile.println();
			    sgmFile.printf("%10.6f",FrSgm.angle[i].angle2);
			    sgmFile.printf("%12.5E",FrSgm.angle[i].force2);
			    sgmFile.println();
			    sgmFile.printf("%10.6f",FrSgm.angle[i].angle3);
			    if(id[ida[number][0]][1]>=0 && id[ida[number][1]][1]>=0 && id[ida[number][2]][1]>=0) {
				sgmFile.printf("%12.5E",0.0);
			    } else {
				sgmFile.printf("%12.5E",FrSgm.angle[i].force3);
			    };
			};
			sgmFile.println();
			number++;
		    };
		};
		for(int i=0; i<ToSgm.numAngles; i++){
		    ida[number][0]=idt[ToSgm.angle[i].atomi];
		    ida[number][1]=idt[ToSgm.angle[i].atomj];
		    ida[number][2]=idt[ToSgm.angle[i].atomk];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2]) found=true;
			if(ida[j][0]==ida[number][2] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][0]) found=true;
		    };
		    if(!found) {
			sgmFile.printf("%5d",number+1);
			sgmFile.printf("%5d",ida[number][0]+1);
			sgmFile.printf("%5d",ida[number][1]+1);
			sgmFile.printf("%5d",ida[number][2]+1);
			sgmFile.printf("%5d",ToSgm.angle[i].type);
			if(id[ida[number][0]][0]>=0 && id[ida[number][1]][0]>=0 && id[ida[number][2]][0]>=0) {
			    sgmFile.printf("%3d",1); sgmFile.printf("%1d",1);
			    sgmFile.printf("%1d",ToSgm.angle[i].source); sgmFile.println();
			    sgmFile.printf("%10.6f",ToSgm.angle[i].angle1);
			    sgmFile.printf("%12.5E",0.0);
			    sgmFile.println();
			    sgmFile.printf("%10.6f",ToSgm.angle[i].angle2);
			    sgmFile.printf("%12.5E",0.0);
			} else {
			    sgmFile.printf("%3d",ToSgm.angle[i].source); sgmFile.printf("%1d",ToSgm.angle[i].source);
			    sgmFile.printf("%1d",ToSgm.angle[i].source); sgmFile.println();
			    sgmFile.printf("%10.6f",ToSgm.angle[i].angle1);
			    sgmFile.printf("%12.5E",ToSgm.angle[i].force1);
			    sgmFile.println();
			    sgmFile.printf("%10.6f",ToSgm.angle[i].angle2);
			    sgmFile.printf("%12.5E",ToSgm.angle[i].force2);
			};
			sgmFile.println();
			sgmFile.printf("%10.6f",ToSgm.angle[i].angle3);
			sgmFile.printf("%12.5E",ToSgm.angle[i].force3);
			sgmFile.println();
			number++;
		    };
		};
		number=0;
	        jfound=-1;
		ida = new int[FrSgm.numTorsions+ToSgm.numTorsions][4];
		for(int i=0; i<FrSgm.numTorsions; i++){
		    ida[number][0]=idf[FrSgm.torsion[i].atomi];
		    ida[number][1]=idf[FrSgm.torsion[i].atomj];
		    ida[number][2]=idf[FrSgm.torsion[i].atomk];
		    ida[number][3]=idf[FrSgm.torsion[i].atoml];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][3]) found=true;
			if(ida[j][0]==ida[number][3] && ida[j][1]==ida[number][2] && ida[j][2]==ida[number][1] && ida[j][3]==ida[number][0]) found=true;
		    };
		    if(!found) {
			sgmFile.printf("%5d",number+1);
			sgmFile.printf("%5d",ida[number][0]+1);
			sgmFile.printf("%5d",ida[number][1]+1);
			sgmFile.printf("%5d",ida[number][2]+1);
			sgmFile.printf("%5d",ida[number][3]+1);
			sgmFile.printf("%5d",FrSgm.torsion[i].type);
			jfound=-1;
			for(int j=0; j<ToSgm.numTorsions; j++){
			    if(idt[ToSgm.torsion[j].atomi]==ida[number][0] && idt[ToSgm.torsion[j].atomj]==ida[number][1]
			       && idt[ToSgm.torsion[j].atomk]==ida[number][2] && idt[ToSgm.torsion[j].atoml]==ida[number][3]) jfound=j;
			    if(idt[ToSgm.torsion[j].atoml]==ida[number][0] && idt[ToSgm.torsion[j].atomk]==ida[number][1]
			       && idt[ToSgm.torsion[j].atomj]==ida[number][2] && idt[ToSgm.torsion[j].atomi]==ida[number][3]) jfound=j;
			};
			if(jfound>=0) {
			    sgmFile.printf("%3d",FrSgm.torsion[i].source);
			    sgmFile.printf("%1d",FrSgm.torsion[i].source);
			    sgmFile.printf("%1d",ToSgm.torsion[jfound].source); sgmFile.println();
			    sgmFile.printf("%3d",FrSgm.torsion[i].multi1);
			    sgmFile.printf("%10.6f",FrSgm.torsion[i].torsion1);
			    sgmFile.printf("%12.5E",FrSgm.torsion[i].force1);
			    sgmFile.println();
			    sgmFile.printf("%3d",FrSgm.torsion[i].multi2);
			    sgmFile.printf("%10.6f",FrSgm.torsion[i].torsion2);
			    sgmFile.printf("%12.5E",FrSgm.torsion[i].force2);
			    sgmFile.println();
			    sgmFile.printf("%3d",ToSgm.torsion[jfound].multi1);
			    sgmFile.printf("%10.6f",ToSgm.torsion[jfound].torsion1);
			    sgmFile.printf("%12.5E",ToSgm.torsion[jfound].force1);
			} else {
			    sgmFile.printf("%3d",FrSgm.torsion[i].source);
			    sgmFile.printf("%1d",FrSgm.torsion[i].source);
			    if(id[ida[number][0]][1]>=0 && id[ida[number][1]][1]>=0 && id[ida[number][2]][1]>=0 && id[ida[number][3]][1]>=0) {
				sgmFile.printf("%1d",1);
			    } else { 
				sgmFile.printf("%1d",FrSgm.torsion[i].source); 
			    }; 
			    sgmFile.println();
			    sgmFile.printf("%3d",FrSgm.torsion[i].multi1);
			    sgmFile.printf("%10.6f",FrSgm.torsion[i].torsion1);
			    sgmFile.printf("%12.5E",FrSgm.torsion[i].force1);
			    sgmFile.println();
			    sgmFile.printf("%3d",FrSgm.torsion[i].multi2);
			    sgmFile.printf("%10.6f",FrSgm.torsion[i].torsion2);
			    sgmFile.printf("%12.5E",FrSgm.torsion[i].force2);
			    sgmFile.println();
			    sgmFile.printf("%3d",FrSgm.torsion[i].multi3);
			    sgmFile.printf("%10.6f",FrSgm.torsion[i].torsion3);
			    if(id[ida[number][0]][1]>=0 && id[ida[number][1]][1]>=0 && id[ida[number][2]][1]>=0 && id[ida[number][3]][1]>=0) {
				sgmFile.printf("%12.5E",0.0);
			    } else {
				sgmFile.printf("%12.5E",FrSgm.torsion[i].force3);
			    };
			};
			sgmFile.println();
			number++;
		    };
		};
		for(int i=0; i<ToSgm.numTorsions; i++){
		    ida[number][0]=idt[ToSgm.torsion[i].atomi];
		    ida[number][1]=idt[ToSgm.torsion[i].atomj];
		    ida[number][2]=idt[ToSgm.torsion[i].atomk];
		    ida[number][3]=idt[ToSgm.torsion[i].atoml];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][3]) found=true;
			if(ida[j][0]==ida[number][3] && ida[j][1]==ida[number][2] && ida[j][2]==ida[number][1] && ida[j][3]==ida[number][0]) found=true;
		    };
		    if(!found) {
			sgmFile.printf("%5d",number+1);
			sgmFile.printf("%5d",ida[number][0]+1);
			sgmFile.printf("%5d",ida[number][1]+1);
			sgmFile.printf("%5d",ida[number][2]+1);
			sgmFile.printf("%5d",ida[number][3]+1);
			sgmFile.printf("%5d",ToSgm.torsion[i].type);
			if(id[ida[number][0]][0]>=0 && id[ida[number][1]][0]>=0 && id[ida[number][2]][0]>=0 && id[ida[number][3]][0]>=0) {
			    sgmFile.printf("%1d",1); sgmFile.printf("%1d",1);
			    sgmFile.printf("%3d",ToSgm.torsion[i].source); sgmFile.println();
			    sgmFile.printf("%3d",ToSgm.torsion[i].multi1);
			    sgmFile.printf("%10.6f",ToSgm.torsion[i].torsion1);
			    sgmFile.printf("%12.5E",0.0);
			    sgmFile.println();
			    sgmFile.printf("%3d",ToSgm.torsion[i].multi2);
			    sgmFile.printf("%10.6f",ToSgm.torsion[i].torsion2);
			    sgmFile.printf("%12.5E",0.0);
			} else {
			    sgmFile.printf("%3d",ToSgm.torsion[i].source); sgmFile.printf("%1d",ToSgm.torsion[i].source);
			    sgmFile.printf("%1d",ToSgm.torsion[i].source); sgmFile.println();
			    sgmFile.printf("%3d",ToSgm.torsion[i].multi1);
			    sgmFile.printf("%10.6f",ToSgm.torsion[i].torsion1);
			    sgmFile.printf("%12.5E",ToSgm.torsion[i].force1);
			    sgmFile.println();
			    sgmFile.printf("%3d",ToSgm.torsion[i].multi2);
			    sgmFile.printf("%10.6f",ToSgm.torsion[i].torsion2);
			    sgmFile.printf("%12.5E",ToSgm.torsion[i].force2);
			};
			sgmFile.println();
			sgmFile.printf("%3d",ToSgm.torsion[i].multi3);
			sgmFile.printf("%10.6f",ToSgm.torsion[i].torsion3);
			sgmFile.printf("%12.5E",ToSgm.torsion[i].force3);
			sgmFile.println();
			number++;
		    };
		};
		number=0;
	        jfound=-1;
		ida = new int[FrSgm.numImpropers+ToSgm.numImpropers][4];
		for(int i=0; i<FrSgm.numImpropers; i++){
		    ida[number][0]=idf[FrSgm.improper[i].atomi];
		    ida[number][1]=idf[FrSgm.improper[i].atomj];
		    ida[number][2]=idf[FrSgm.improper[i].atomk];
		    ida[number][3]=idf[FrSgm.improper[i].atoml];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][3]) found=true;
			if(ida[j][0]==ida[number][1] && ida[j][1]==ida[number][0] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][3]) found=true;
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][3] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][1]) found=true;
			if(ida[j][0]==ida[number][3] && ida[j][1]==ida[number][0] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][1]) found=true;
			if(ida[j][0]==ida[number][1] && ida[j][1]==ida[number][3] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][0]) found=true;
			if(ida[j][0]==ida[number][3] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][0]) found=true;
		    };
		    if(!found) {
			sgmFile.printf("%5d",number+1);
			sgmFile.printf("%5d",ida[number][0]+1);
			sgmFile.printf("%5d",ida[number][1]+1);
			sgmFile.printf("%5d",ida[number][2]+1);
			sgmFile.printf("%5d",ida[number][3]+1);
			sgmFile.printf("%5d",FrSgm.improper[i].type);
			jfound=-1;
			for(int j=0; j<ToSgm.numImpropers; j++){
			    if(idt[ToSgm.improper[j].atomi]==ida[number][0] && idt[ToSgm.improper[j].atomj]==ida[number][1]
			       && idt[ToSgm.improper[j].atomk]==ida[number][2] && idt[ToSgm.improper[j].atoml]==ida[number][3]) jfound=j;
			    if(idt[ToSgm.improper[j].atomi]==ida[number][1] && idt[ToSgm.improper[j].atomj]==ida[number][0]
			       && idt[ToSgm.improper[j].atomk]==ida[number][2] && idt[ToSgm.improper[j].atoml]==ida[number][3]) jfound=j;
			    if(idt[ToSgm.improper[j].atomi]==ida[number][0] && idt[ToSgm.improper[j].atomj]==ida[number][3]
			       && idt[ToSgm.improper[j].atomk]==ida[number][2] && idt[ToSgm.improper[j].atoml]==ida[number][1]) jfound=j;
			    if(idt[ToSgm.improper[j].atomi]==ida[number][3] && idt[ToSgm.improper[j].atomj]==ida[number][0]
			       && idt[ToSgm.improper[j].atomk]==ida[number][2] && idt[ToSgm.improper[j].atoml]==ida[number][1]) jfound=j;
			    if(idt[ToSgm.improper[j].atomi]==ida[number][1] && idt[ToSgm.improper[j].atomj]==ida[number][3]
			       && idt[ToSgm.improper[j].atomk]==ida[number][2] && idt[ToSgm.improper[j].atoml]==ida[number][0]) jfound=j;
			    if(idt[ToSgm.improper[j].atomi]==ida[number][3] && idt[ToSgm.improper[j].atomj]==ida[number][1]
			       && idt[ToSgm.improper[j].atomk]==ida[number][2] && idt[ToSgm.improper[j].atoml]==ida[number][0]) jfound=j;
			};
			if(jfound>=0) {
			    sgmFile.printf("%3d",FrSgm.improper[i].source);
			    sgmFile.printf("%1d",FrSgm.improper[i].source);
			    sgmFile.printf("%1d",ToSgm.improper[jfound].source); sgmFile.println();
			    sgmFile.printf("%3d",FrSgm.improper[i].multi1);
			    sgmFile.printf("%10.6f",FrSgm.improper[i].improper1);
			    sgmFile.printf("%12.5E",FrSgm.improper[i].force1);
			    sgmFile.println();
			    sgmFile.printf("%3d",FrSgm.improper[i].multi2);
			    sgmFile.printf("%10.6f",FrSgm.improper[i].improper2);
			    sgmFile.printf("%12.5E",FrSgm.improper[i].force2);
			    sgmFile.println();
			    sgmFile.printf("%3d",ToSgm.improper[jfound].multi1);
			    sgmFile.printf("%10.6f",ToSgm.improper[jfound].improper1);
			    sgmFile.printf("%12.5E",ToSgm.improper[jfound].force1);
			} else {
			    sgmFile.printf("%3d",FrSgm.improper[i].source);
			    sgmFile.printf("%1d",FrSgm.improper[i].source);
			    if(id[ida[number][0]][1]>=0 && id[ida[number][1]][1]>=0 && id[ida[number][2]][1]>=0 && id[ida[number][3]][1]>=0) {
				sgmFile.printf("%1d",1);
			    } else {
				sgmFile.printf("%1d",FrSgm.improper[i].source); 
			    }; 
			    sgmFile.println();
			    sgmFile.printf("%3d",FrSgm.improper[i].multi1);
			    sgmFile.printf("%10.6f",FrSgm.improper[i].improper1);
			    sgmFile.printf("%12.5E",FrSgm.improper[i].force1);
			    sgmFile.println();
			    sgmFile.printf("%3d",FrSgm.improper[i].multi2);
			    sgmFile.printf("%10.6f",FrSgm.improper[i].improper2);
			    sgmFile.printf("%12.5E",FrSgm.improper[i].force2);
			    sgmFile.println();
			    sgmFile.printf("%3d",FrSgm.improper[i].multi3);
			    sgmFile.printf("%10.6f",FrSgm.improper[i].improper3);
			    if(id[ida[number][0]][1]>=0 && id[ida[number][1]][1]>=0 && id[ida[number][2]][1]>=0 && id[ida[number][3]][1]>=0) {
				sgmFile.printf("%12.5E",0.0);
			    } else {
				sgmFile.printf("%12.5E",FrSgm.improper[i].force3);
			    };
			};
			sgmFile.println();
			number++;
		    };
		};
		for(int i=0; i<ToSgm.numImpropers; i++){
		    ida[number][0]=idt[ToSgm.improper[i].atomi];
		    ida[number][1]=idt[ToSgm.improper[i].atomj];
		    ida[number][2]=idt[ToSgm.improper[i].atomk];
		    ida[number][3]=idt[ToSgm.improper[i].atoml];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][3]) found=true;
			if(ida[j][0]==ida[number][1] && ida[j][1]==ida[number][0] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][3]) found=true;
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][3] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][1]) found=true;
			if(ida[j][0]==ida[number][3] && ida[j][1]==ida[number][0] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][1]) found=true;
			if(ida[j][0]==ida[number][1] && ida[j][1]==ida[number][3] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][0]) found=true;
			if(ida[j][0]==ida[number][3] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2] && ida[j][3]==ida[number][0]) found=true;
		    };
		    if(!found) {
			sgmFile.printf("%5d",number+1);
			sgmFile.printf("%5d",ida[number][0]+1);
			sgmFile.printf("%5d",ida[number][1]+1);
			sgmFile.printf("%5d",ida[number][2]+1);
			sgmFile.printf("%5d",ida[number][3]+1);
			sgmFile.printf("%5d",ToSgm.improper[i].type);
			if(id[ida[number][0]][0]>=0 && id[ida[number][1]][0]>=0 && id[ida[number][2]][0]>=0 && id[ida[number][3]][0]>=0) {
			    sgmFile.printf("%3d",1); sgmFile.printf("%1d",1);
			    sgmFile.printf("%1d",ToSgm.improper[i].source); sgmFile.println();
			    sgmFile.printf("%3d",ToSgm.improper[i].multi1);
			    sgmFile.printf("%10.6f",ToSgm.improper[i].improper1);
			    sgmFile.printf("%12.5E",0.0);
			    sgmFile.println();
			    sgmFile.printf("%3d",ToSgm.improper[i].multi2);
			    sgmFile.printf("%10.6f",ToSgm.improper[i].improper2);
			    sgmFile.printf("%12.5E",0.0);
			} else {
			    sgmFile.printf("%3d",ToSgm.improper[i].source); sgmFile.printf("%1d",ToSgm.improper[i].source);
			    sgmFile.printf("%1d",ToSgm.improper[i].source); sgmFile.println();
			    sgmFile.printf("%3d",ToSgm.improper[i].multi1);
			    sgmFile.printf("%10.6f",ToSgm.improper[i].improper1);
			    sgmFile.printf("%12.5E",ToSgm.improper[i].force1);
			    sgmFile.println();
			    sgmFile.printf("%3d",ToSgm.improper[i].multi2);
			    sgmFile.printf("%10.6f",ToSgm.improper[i].improper2);
			    sgmFile.printf("%12.5E",ToSgm.improper[i].force2);
			};
			sgmFile.println();
			sgmFile.printf("%3d",ToSgm.improper[i].multi3);
			sgmFile.printf("%10.6f",ToSgm.improper[i].improper3);
			sgmFile.printf("%12.5E",ToSgm.improper[i].force3);
			sgmFile.println();
			number++;
		    };
		};
		number=0;
	        jfound=-1;
		ida = new int[FrSgm.numZmatrix+ToSgm.numZmatrix][4];
		for(int i=0; i<FrSgm.numZmatrix; i++){
		    ida[number][0]=idf[FrSgm.zmatrix[i].atomi];
		    ida[number][1]=idf[FrSgm.zmatrix[i].atomj];
		    ida[number][2]=idf[FrSgm.zmatrix[i].atomk];
		    ida[number][3]=idf[FrSgm.zmatrix[i].atoml];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2]
			   && ida[j][3]==ida[number][3]) found=true;
		    };
		    if(!found) {
			sgmFile.printf("%5d",number+1);
			sgmFile.printf("%5d",ida[number][0]+1);
			sgmFile.printf("%5d",ida[number][1]+1);
			sgmFile.printf("%5d",ida[number][2]+1);
			sgmFile.printf("%5d",ida[number][3]+1);
			sgmFile.printf("%12.6f",FrSgm.zmatrix[i].bond);
			sgmFile.printf("%12.6f",FrSgm.zmatrix[i].angle);
			sgmFile.printf("%12.6f",FrSgm.zmatrix[i].torsion);
			sgmFile.println();
			number++;
		    };
		};
		for(int i=0; i<ToSgm.numZmatrix; i++){
		    ida[number][0]=idt[ToSgm.zmatrix[i].atomi];
		    ida[number][1]=idt[ToSgm.zmatrix[i].atomj];
		    ida[number][2]=idt[ToSgm.zmatrix[i].atomk];
		    ida[number][3]=idt[ToSgm.zmatrix[i].atoml];
		    found=false;
		    for(int j=0; j<number; j++){
			if(ida[j][0]==ida[number][0] && ida[j][1]==ida[number][1] && ida[j][2]==ida[number][2]
			   && ida[j][3]==ida[number][3]) found=true;
			if(ida[j][0]==ida[number][3] && ida[j][1]==ida[number][2] && ida[j][2]==ida[number][1]
			   && ida[j][3]==ida[number][0]) found=true;
		    };
		    if(!found) {
			sgmFile.printf("%5d",number+1);
			sgmFile.printf("%5d",ida[number][0]+1);
			sgmFile.printf("%5d",ida[number][1]+1);
			sgmFile.printf("%5d",ida[number][2]+1);
			sgmFile.printf("%5d",ida[number][3]+1);
			sgmFile.printf("%12.6f",ToSgm.zmatrix[i].bond);
			sgmFile.printf("%12.6f",ToSgm.zmatrix[i].angle);
			sgmFile.printf("%12.6f",ToSgm.zmatrix[i].torsion);
			sgmFile.println();
			number++;
		    };
		};
		sgmFile.close();
	    } catch (Exception ee) { ee.printStackTrace(); };
	    Segment tempSgm = new Segment();
	    tempSgm.read("TEMP.sgm");
	    tempSgm.order();
	    tempSgm.write("NEW.sgm");
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
		System.out.println("SELECT TO JOIN atmList"); 
		selected=atmList.getSelectedIndex(); 
		j=atmList.getSelectedIndex();
		System.out.println("SELECT TO JOIN "+j);
	    } else {
		System.out.println("SELECT TO JOIN NOT atmList"); 
		selected = -1;
	    };
	};
	if(mouse.getModifiers()==(MouseEvent.BUTTON1_MASK+MouseEvent.CTRL_MASK)){
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
	if(mouse.getModifiers()==(MouseEvent.BUTTON3_MASK+MouseEvent.CTRL_MASK)){
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
	    if(mouse.getSource()==sgmList){
		j=sgmList.getSelectedIndex();
		fileName=sgmDef[j].Dir+sgmDef[j].Name;
		if(frBool){
		    FrSgm.read(fileName);
		    frLabel.setText(" "+sgmDef[j].Name.substring(0,sgmDef[j].Name.indexOf(".sgm")));
		    frRead=true;
		    frName=sgmDef[j].Name.substring(0,sgmDef[j].Name.indexOf(".sgm"));
		} else {
		    ToSgm.read(fileName);
		    toLabel.setText(" "+sgmDef[j].Name.substring(0,sgmDef[j].Name.indexOf(".sgm")));
		    toRead=true;
		    toName=sgmDef[j].Name.substring(0,sgmDef[j].Name.indexOf(".sgm"));
		};
		if(frRead && !toRead){
		    id = new int[FrSgm.numAtoms][2];
		    atmNumber=0;
		    for(int i=0; i<FrSgm.numAtoms; i++){
			id[i][0]=i; id[i][1]=-1; atmNumber++;
		    };
		};
		if(frRead && toRead){
		    boolean[] frFnd = new boolean[FrSgm.numAtoms];
		    for(int i=0; i<FrSgm.numAtoms; i++){frFnd[i]=false;};
		    boolean[] toFnd = new boolean[ToSgm.numAtoms];
		    for(int i=0; i<ToSgm.numAtoms; i++){toFnd[i]=false;};
		    int num = FrSgm.numAtoms+ToSgm.numAtoms;
		    id = new int[num][2];
		    atmNumber=0;
		    for(int i=0; i<num; i++){ id[i][0]=-1; id[i][1]=-1;};
		    for(int i=0; i<FrSgm.numAtoms; i++){
			if(!frFnd[i]){
			    for(int k=0; k<ToSgm.numAtoms; k++){
			    	if(!toFnd[k]){
				    if(FrSgm.atom[i].Name.equals(ToSgm.atom[k].Name)){
					id[atmNumber][0]=i; id[atmNumber][1]=k; atmNumber++; frFnd[i]=true; toFnd[k]=true;
				    };
			    	};
			    };
			};
			if(!frFnd[i]){
			    id[atmNumber][0]=i; atmNumber++; frFnd[i]=true;
			};
		    };
		    for(int k=0; k<ToSgm.numAtoms; k++){
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
			if(ToSgm.atom[id[j+1][1]].Name.charAt(4)==' '){
			    ToSgm.atom[id[j+1][1]].Name=ToSgm.atom[id[j+1][1]].Name.substring(0,4)+"t";
			};
			selected=-1;
			id[j][1]=-1; id[j+1][0]=-1;  atmNumber++;
		    }; 
		    if(id[j][0]>=0 && id[j][1]<0){
			selected=j;
		    };
		    if(id[j][0]<0 && id[j][1]>=0 && selected>=0){
			if(id[j][0]==-1 && id[j][1]>=0 && id[selected][0]>=0 && id[selected][1]==-1){
			    id[selected][1]=id[j][1];
			    if(FrSgm.atom[id[selected][0]].Name.substring(0,4).equals(ToSgm.atom[id[selected][1]].Name.substring(0,4)) &&
			       ToSgm.atom[id[selected][1]].Name.charAt(4)=='t') ToSgm.atom[id[selected][1]].Name=ToSgm.atom[id[selected][1]].Name.substring(0,4)+" ";
			    atmNumber--;
			    for(int k=j; k<atmNumber; k++){
				id[k][0]=id[k+1][0];
				id[k][1]=id[k+1][1];
			    };
			};
			if(id[j][0]>=0 && id[j][1]==-1 && id[selected][0]==-1 && id[selected][1]>=0){
			    System.out.println("Join "+j+" and "+selected);
			    id[j][1]=id[selected][1];
			    if(FrSgm.atom[id[j][0]].Name.substring(0,4).equals(ToSgm.atom[id[j][1]].Name.substring(0,4)) &&
			       ToSgm.atom[id[j][1]].Name.charAt(4)=='t') ToSgm.atom[id[j][1]].Name=ToSgm.atom[id[j][1]].Name.substring(0,4)+" ";
			    atmNumber--;
			    for(int k=selected; k<atmNumber; k++){
				id[k][0]=id[k+1][0];
				id[k][1]=id[k+1][1];
			    };
			};
			atomListUpdate();
		    };
		};
		//		selected= -1;
		atomListUpdate();
	    };
	};
	if(mouse.getModifiers()==MouseEvent.BUTTON3_MASK){
	    if(mouse.getSource()==sgmList){
		if(frBool) {
		    toLabel.setBackground(Color.yellow);
		    frLabel.setBackground(Color.lightGray);
		    toLabel.setForeground(Color.blue);
		    frLabel.setForeground(Color.darkGray);
		} else {
		    frLabel.setBackground(Color.yellow);
		    toLabel.setBackground(Color.lightGray);
		    frLabel.setForeground(Color.blue);
		    toLabel.setForeground(Color.darkGray);
		};
	      	frBool=!frBool;
	    };
	};
    };
    
    public void mouseEntered(MouseEvent mouse){}

    public void mouseExited(MouseEvent mouse){}
  
}


