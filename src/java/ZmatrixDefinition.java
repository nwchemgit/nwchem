import java.io.*;

public class ZmatrixDefinition{

    public int atomi,atomj,atomk,atoml;
    public double bond,angle,torsion;

    public ZmatrixDefinition(String card){
	atomi=Integer.parseInt(card.substring(5,10).trim());
	atomj=Integer.parseInt(card.substring(10,15).trim());
	atomk=Integer.parseInt(card.substring(15,20).trim());
	atoml=Integer.parseInt(card.substring(20,25).trim());
	bond=Double.valueOf(card.substring(25,37)).doubleValue();
	angle=Double.valueOf(card.substring(37,49)).doubleValue();
	torsion=Double.valueOf(card.substring(49,61)).doubleValue();
    }

}
