public class Rule{

    public String type;
    public int source;
    boolean redefined;

    public Rule(String card, String card2, String card3, int isource){
	source=isource;
	type=card.substring(0,2);
	redefined=false;
    }

}
