public class Rule{

    public String type;
    public int source;
    boolean redefining, redefined, selected;

    public Rule(String card, String card2, String card3, int isource){
	source=isource;
	type=card.substring(0,3);
	redefining=false; redefined=false; selected=true;
    }

}
