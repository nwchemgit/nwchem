import java.io.*;
public class Segment{

    public int numAtoms, numBonds, numAngles, numTorsions, numImpropers;

    public AtomDefinition[] atom;
    public BondDefinition[] bond;

    public Segment(){
    };

    void SegmentRead(String fileName){
	System.out.println("Reading file "+fileName);
	try{
	    BufferedReader br = new BufferedReader(new FileReader(fileName));
	    String card1;
	    String card2;
	    card1=br.readLine();
	    while(card1.startsWith("$") || card1.startsWith("#")){card1=br.readLine();};
	    System.out.println(card1);
	    numAtoms=Integer.parseInt(card1.substring(1,5).trim());
	    numBonds=Integer.parseInt(card1.substring(6,10).trim());
	    numAngles=Integer.parseInt(card1.substring(11,15).trim());
	    numTorsions=Integer.parseInt(card1.substring(16,20).trim());
	    numImpropers=Integer.parseInt(card1.substring(21,25).trim());
	    System.out.println("Number of atoms is "+numAtoms);
	    atom = new AtomDefinition[numAtoms];
	    bond = new BondDefinition[numBonds];
	    for(int i=0; i<numAtoms; i++){
		card1=br.readLine();
		card2=br.readLine();
		atom[i] = new AtomDefinition(card1.substring(5,9),card1.substring(10,14),card1.substring(15,19),card1.substring(20,24));
	    };
	    for(int i=0; i<numBonds; i++){
		card1=br.readLine();
		card2=br.readLine();
	    };
	    for(int i=0; i<numAngles; i++){
		card1=br.readLine();
		card2=br.readLine();
	    };
	    for(int i=0; i<numTorsions; i++){
		card1=br.readLine();
		card2=br.readLine();
	    };
	    for(int i=0; i<numImpropers; i++){
		card1=br.readLine();
		card2=br.readLine();
	    };
	    br.close();
	} catch(Exception e) {e.printStackTrace();};
    };

}
