import java.io.*;

public class BondDefinition{

    public int atomi,atomj;
    public int type, source;
    public double bond1,bond2,bond3,force1,force2,force3;

    public BondDefinition(String card1, String card2){
	atomi=Integer.parseInt(card1.substring(5,10).trim());
	atomj=Integer.parseInt(card1.substring(10,15).trim());
	type=Integer.parseInt(card1.substring(15,20).trim());
	source=Integer.parseInt(card1.substring(20,25).trim());
	bond1=Double.valueOf(card2.substring(0,12)).doubleValue();
	bond2=Double.valueOf(card2.substring(24,36)).doubleValue();
	bond3=Double.valueOf(card2.substring(48,60)).doubleValue();
	force1=Double.valueOf(card2.substring(12,24)).doubleValue();
	force2=Double.valueOf(card2.substring(36,48)).doubleValue();
	force3=Double.valueOf(card2.substring(60,72)).doubleValue();
    };

    public BondDefinition(String card){
	atomi=Integer.parseInt(card.substring(0,5).trim());
	atomj=Integer.parseInt(card.substring(5,10).trim());
    };

}
