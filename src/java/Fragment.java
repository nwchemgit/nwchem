import java.io.*;
public class Fragment{

    public int numAtoms, numBonds, numZmatrix;

    public AtomDefinition[] atom;
    public BondDefinition[] bond;
    public ZmatrixDefinition[] zmatrix;

    public Fragment(){
    };

    void read(String fileName){
	//	System.out.println("Reading file "+fileName);
	try{
	    BufferedReader br = new BufferedReader(new FileReader(fileName));
	    String card;
	    card=br.readLine();
	    while(card.startsWith("$") || card.startsWith("#")){card=br.readLine();};
	    numAtoms=Integer.parseInt(card.substring(0,5).trim());
	    int maxBonds=10*numAtoms;
	    int maxZmatrix=10*numAtoms;
	    try{
		numZmatrix=Integer.parseInt(card.substring(61,66).trim());
	    } catch(Exception ez) {numZmatrix=0;};
	    atom = new AtomDefinition[numAtoms];
	    bond = new BondDefinition[maxBonds];
	    zmatrix = new ZmatrixDefinition[maxZmatrix];
	    for(int i=0; i<numAtoms; i++){
		card=br.readLine();
		atom[i] = new AtomDefinition(card);
	    };
	    card=br.readLine();
	    System.out.println("Size is "+card.length());
	    //	    for(int i=0; i<numBonds; i++){
	    //	card=br.readLine();
	    //	bond[i] = new BondDefinition(card);
	    //};
	    //for(int i=0; i<numZmatrix; i++){
	    //	card=br.readLine();
	    //	zmatrix[i] = new ZmatrixDefinition(card);
	    //};
	    br.close();
	} catch(Exception e) {e.printStackTrace();};
    };

}
