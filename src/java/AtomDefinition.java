public class AtomDefinition{

    public String Name;
    public String Type1, Type2, Type3;
    public int cgroup, pgroup, link, type;
    public double q1,q2,q3,p1,p2,p3;

    public AtomDefinition(String card1, String card2){
	Name=card1.substring(5,10);
	Type1=card1.substring(11,17);
	Type2=card1.substring(17,23);
	Type3=card1.substring(23,29);
	cgroup=Integer.parseInt(card1.substring(29,33).trim());
	pgroup=Integer.parseInt(card1.substring(33,37).trim());
	link=Integer.parseInt(card1.substring(37,41).trim());
	type=Integer.parseInt(card1.substring(41,45).trim());
	q1=Double.valueOf(card2.substring(0,12)).doubleValue();
	q2=Double.valueOf(card2.substring(24,36)).doubleValue();
	q3=Double.valueOf(card2.substring(48,60)).doubleValue();
	p1=Double.valueOf(card2.substring(12,24)).doubleValue();
	p2=Double.valueOf(card2.substring(36,48)).doubleValue();
	p3=Double.valueOf(card2.substring(60,72)).doubleValue();
    }
    public AtomDefinition(String card){
	Name=card.substring(5,10);
	Type1=card.substring(11,16);
	link=Integer.parseInt(card.substring(17,22).trim());
	type=Integer.parseInt(card.substring(22,27).trim());
	cgroup=Integer.parseInt(card.substring(27,32).trim());
	pgroup=Integer.parseInt(card.substring(32,37).trim());
	q1=Double.valueOf(card.substring(37,49)).doubleValue();
	p1=Double.valueOf(card.substring(49,61)).doubleValue();
    };

}
