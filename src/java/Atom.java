public class Atom{

    public double epsilon, epsilon3, rstar, rstar3, mass;
    public String type;
    public int number, source;
    boolean redefining, redefined, selected;

    public Atom(String card, String card2, int isource){
	source=isource;
	type=card.substring(0,3);
	mass=Double.valueOf(card.substring(3,15)).doubleValue();
	epsilon=Double.valueOf(card.substring(15,27)).doubleValue();
       	rstar=Double.valueOf(card.substring(27,39)).doubleValue();
	number=Integer.parseInt(card2.substring(10,15).trim());
       	epsilon3=Double.valueOf(card2.substring(15,27)).doubleValue();
       	rstar3=Double.valueOf(card2.substring(27,39)).doubleValue();
	redefining=false; redefined=false; selected=true;
    }

}
