import java.io.*;
public class Segment{

    public int numAtoms, numBonds, numAngles, numTorsions, numImpropers, numZmatrix;

    public AtomDefinition[] atom;
    public BondDefinition[] bond;
    public AngleDefinition[] angle;
    public TorsionDefinition[] torsion;
    public ImproperDefinition[] improper;
    public ZmatrixDefinition[] zmatrix;

    public Segment(){
    };

    void SegmentRead(String fileName){
	//	System.out.println("Reading file "+fileName);
	try{
	    BufferedReader br = new BufferedReader(new FileReader(fileName));
	    String card1;
	    String card2;
	    card1=br.readLine();
	    while(card1.startsWith("$") || card1.startsWith("#")){card1=br.readLine();};
	    numAtoms=Integer.parseInt(card1.substring(0,5).trim());
	    numBonds=Integer.parseInt(card1.substring(5,10).trim());
	    numAngles=Integer.parseInt(card1.substring(10,15).trim());
	    numTorsions=Integer.parseInt(card1.substring(15,20).trim());
	    numImpropers=Integer.parseInt(card1.substring(20,25).trim());
	    try{
		numZmatrix=Integer.parseInt(card1.substring(61,66).trim());
	    } catch(Exception ez) {numZmatrix=0;};
	    atom = new AtomDefinition[numAtoms];
	    bond = new BondDefinition[numBonds];
	    angle = new AngleDefinition[numAngles];
	    torsion = new TorsionDefinition[numTorsions];
	    improper = new ImproperDefinition[numImpropers];
	    zmatrix = new ZmatrixDefinition[numZmatrix];
	    for(int i=0; i<numAtoms; i++){
		card1=br.readLine();
		card2=br.readLine();
		atom[i] = new AtomDefinition(card1,card2);
	    };
	    for(int i=0; i<numBonds; i++){
		card1=br.readLine();
		card2=br.readLine();
		bond[i] = new BondDefinition(card1,card2);
	    };
	    for(int i=0; i<numAngles; i++){
		card1=br.readLine();
		card2=br.readLine();
		angle[i] = new AngleDefinition(card1,card2);
	    };
	    for(int i=0; i<numTorsions; i++){
		card1=br.readLine();
		card2=br.readLine();
		torsion[i] = new TorsionDefinition(card1,card2);
	    };
	    for(int i=0; i<numImpropers; i++){
		card1=br.readLine();
		card2=br.readLine();
		improper[i] = new ImproperDefinition(card1,card2);
	    };
	    for(int i=0; i<numZmatrix; i++){
		card1=br.readLine();
		zmatrix[i] = new ZmatrixDefinition(card1);
	    };
	    br.close();
	} catch(Exception e) {e.printStackTrace();};
    };

}
