public class Torsion{

    public String type1, type2, type3, type4;
    public int source;
    public boolean redefining, redefined, selected;
    public int number;
    public double[] value = new double[5];
    public double[] force = new double[5];
    public int[] multiplicity = new int[5];
    public boolean repeat = false;

    public Torsion(String card, int isource){
	source=isource;
	type1=card.substring(0,3);
	type2=card.substring(6,9);
	type3=card.substring(12,15);
	type4=card.substring(18,21);
	redefining=false; redefined=false; selected=true;
	number=0;
	multiplicity[number]=Integer.parseInt(card.substring(45,50).trim());
	value[number]=Double.valueOf(card.substring(24,33)).doubleValue();
	force[number]=Double.valueOf(card.substring(33,45)).doubleValue();
	redefining=false; redefined=false;
        repeat=(multiplicity[number]<0);
    }

    public void add(String card){
	number++;
	multiplicity[number]=Integer.parseInt(card.substring(45,50).trim());
	value[number]=Double.valueOf(card.substring(24,33)).doubleValue();
	force[number]=Double.valueOf(card.substring(33,45)).doubleValue();
        repeat=(multiplicity[number]<0);
    }
}
