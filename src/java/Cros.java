public class Cros{

    public double epsilon, epsilon3, rstar, rstar3;
    public String type1, type2;
    public int source;
    boolean redefining, redefined, selected;

    public Cros(String card, String card2, int isource){
	source=isource;
	type1=card.substring(0,3);
	type2=card.substring(6,9);
	redefining=false; redefined=false; selected=true;
    }

}
