public class AtomDefinition{

    public String Name;
    public String Type1, Type2, Type3;
    public int cgroup, pgroup, link, type;
    public double q1,q2,q3,p1,p2,p3;

    public AtomDefinition(String card1, String card2, String card3, String card4){
	Name=card1.substring(5,10);
	link=Integer.parseInt(card1.substring(11,16).trim());
	type=Integer.parseInt(card1.substring(16,21).trim());
	cgroup=Integer.parseInt(card1.substring(26,31).trim());
	pgroup=Integer.parseInt(card1.substring(31,36).trim());
	Type1=card2.substring(5,11);
	Type2=card3.substring(5,11);
	Type3=card4.substring(5,11);
	q1=Double.valueOf(card2.substring(11,23)).doubleValue();
	q2=Double.valueOf(card3.substring(11,23)).doubleValue();
	q3=Double.valueOf(card4.substring(11,23)).doubleValue();
	p1=Double.valueOf(card2.substring(23,35)).doubleValue();
	p2=Double.valueOf(card3.substring(23,35)).doubleValue();
	p3=Double.valueOf(card4.substring(23,35)).doubleValue();
    }
    public AtomDefinition(String card1, String card2){
	Name=card1.substring(5,10);
	link=Integer.parseInt(card1.substring(11,16).trim());
	type=Integer.parseInt(card1.substring(16,21).trim());
	cgroup=Integer.parseInt(card1.substring(26,31).trim());
	pgroup=Integer.parseInt(card1.substring(31,36).trim());
	Type1=card2.substring(5,11);
	Type2=Type1;
	Type3=Type1;
	q1=Double.valueOf(card2.substring(11,23)).doubleValue();
	q2=q1;
	q3=q1;
	p1=Double.valueOf(card2.substring(23,35)).doubleValue();
	p2=p2;
	p3=p2;
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
