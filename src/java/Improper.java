public class Improper{

    public double value, force;
    public String type1, type2, type3, type4;
    public int source;
    boolean redefined;

    public Improper(String card, int isource){
	source=isource;
	type1=card.substring(0,2);
	type2=card.substring(6,8);
	type3=card.substring(12,14);
	type4=card.substring(18,20);
	value=Double.valueOf(card.substring(24,33)).doubleValue();
	force=Double.valueOf(card.substring(33,45)).doubleValue();
	redefined=false;
    }

}
