import java.io.*;
public class Segment{

    public int numAtoms, numBonds, numTorsions, numImpropers;

    int maxAtoms=100, maxBonds=100, maxTorsions=100, maxImpropers=100;

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
	    System.out.println("Number of atoms is "+numAtoms);
	    for(int i=0; i<numAtoms; i++){
		card=br.readLine();
		System.out.println(card.substring(5,9));
		card=br.readLine();
	    };
	    br.close();
	} catch(Exception e) {e.printStackTrace();};
    };

}
