import java.io.*;
public class Segment{

    public int numAtoms, numBonds, numAngles, numTorsions, numImpropers;

    public AtomDefinition[] atom;

    public Segment(){
    };

    void SegmentRead(String fileName){
	System.out.println("Reading file "+fileName);
	try{
	    BufferedReader br = new BufferedReader(new FileReader(fileName));
	    String card;
	    card=br.readLine();
	    while(card.startsWith("$") || card.startsWith("#")){card=br.readLine();};
	    System.out.println(card);
	    numAtoms=Integer.parseInt(card.substring(1,5).trim());
	    numBonds=Integer.parseInt(card.substring(6,10).trim());
	    numAngles=Integer.parseInt(card.substring(11,15).trim());
	    numTorsions=Integer.parseInt(card.substring(16,20).trim());
	    numImpropers=Integer.parseInt(card.substring(21,25).trim());
	    System.out.println("Number of atoms is "+numAtoms);
	    atom = new AtomDefinition[numAtoms];
	    for(int i=0; i<numAtoms; i++){
		card=br.readLine();
		atom[i] = new AtomDefinition(card.substring(5,9));
		card=br.readLine();
	    };
	    br.close();
	} catch(Exception e) {e.printStackTrace();};
    };

}
