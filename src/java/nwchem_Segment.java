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
    String dir_u = " ";
    String dir_t = " ";

    char[] file_type = { 's', 'x', 'u', 't'};

    JPanel header = new JPanel();

    DefaultListModel segmentList = new DefaultListModel();
    JList sgmList = new JList(segmentList);
    JScrollPane sgmPane = new JScrollPane(sgmList);

    DefaultListModel atomList = new DefaultListModel();
    JList atmList = new JList(atomList);
    JScrollPane atmPane = new JScrollPane(atmList);

    File[] files;
    String[] dirs = {" ", " ", " ", " "};
    SegmentDefinition[] sgmDef;
    
    JFileChooser chooser;
    ExtensionFilter sgmFilter;

    JLabel frLabel = new JLabel(" From ");
    JLabel toLabel = new JLabel("  To  ");
    boolean frBool = true;

    boolean frRead=false, toRead=false;

    int frIndex, toIndex;
    int atmNumber=0;
    int[][] id;
    
    Segment ToSgm = new Segment();
    Segment FrSgm = new Segment();
    
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
	sgmList.setVisibleRowCount(15);

	atmList.addMouseListener(this);
	atmList.setVisibleRowCount(15);
	
	JButton doneButton = new JButton("Done");
	
	addComponent(header,doneButton,5,0,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.CENTER);
	doneButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){ 
		    setVisible(false); }});

	frLabel.setBackground(Color.yellow);
	toLabel.setBackground(Color.lightGray);
	frLabel.setForeground(Color.blue);
	toLabel.setForeground(Color.darkGray);

	addComponent(header,frLabel,1,0,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.CENTER);
	addComponent(header,toLabel,2,0,1,1,1,1,
		     GridBagConstraints.NONE,GridBagConstraints.CENTER);
	
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
	
	int sgmNumber=0;
	
	for(int idir=0; idir<4; idir++){
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
	for(int idir=0; idir<4; idir++){
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
	addComponent(header,sgmPane,0,1,1,1,1,1,GridBagConstraints.NONE,GridBagConstraints.CENTER);
	
	setLocation(25,225);	
	setSize(1000,800);
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
	header.remove(atmPane);
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
	addComponent(header,atmPane,1,1,2,1,2,2,GridBagConstraints.NONE,GridBagConstraints.CENTER);
	header.validate();
    };

    public void actionPerformed(ActionEvent e){}
    
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
	if(mouse.getModifiers()==MouseEvent.BUTTON1_MASK){
	    if(mouse.getSource()==sgmList){
		j=sgmList.getSelectedIndex();
		fileName=sgmDef[j].Dir+sgmDef[j].Name;
		if(frBool){
		    FrSgm.SegmentRead(fileName);
		    frLabel.setText(" "+sgmDef[j].Name.substring(0,sgmDef[j].Name.indexOf(".sgm")));
		    frRead=true;
		} else {
		    ToSgm.SegmentRead(fileName);
		    toLabel.setText(" "+sgmDef[j].Name.substring(0,sgmDef[j].Name.indexOf(".sgm")));
		    toRead=true;
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
		atomListUpdate();
	    };
	    if(mouse.getSource()==atmList){
		System.out.println(" Mouse event on atmList");
		j=atmList.getSelectedIndex();
		System.out.println(" Event index is "+j+" "+id[j][0]+" "+id[j][1]);
	        if(id[j][0]>=0 && id[j][1]>=0){
		    for(int k=atmNumber; k>j; k--){id[k][0]=id[k-1][0]; id[k][1]=id[k-1][1];};
		    id[j][1]=-1; id[j+1][0]=-1;  atmNumber++;
		};
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


