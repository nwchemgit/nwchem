import java.io.*;

public class BondDefinition{

    public int atomi,atomj;
    public int type, source;
    public int source1,source2,source3;
    public double bond1,bond2,bond3,force1,force2,force3;

    public BondDefinition(String card1, String card2, String card3, String card4){
	atomi=Integer.parseInt(card1.substring(5,10).trim());
	atomj=Integer.parseInt(card1.substring(10,15).trim());
	type=Integer.parseInt(card1.substring(15,20).trim());
	source3=Integer.parseInt(card1.substring(24,25).trim());
	source1=source3;
	source2=source3;
	if(card1.charAt(22)!=' ') source1=Integer.parseInt(card1.substring(20,23).trim());
	if(card1.charAt(23)!=' ') source2=Integer.parseInt(card1.substring(23,24).trim());
	source=source1;
	bond1=Double.valueOf(card2.substring(0,12)).doubleValue();
	bond2=Double.valueOf(card3.substring(0,12)).doubleValue();
	bond3=Double.valueOf(card4.substring(0,12)).doubleValue();
	force1=Double.valueOf(card2.substring(12,24)).doubleValue();
	force2=Double.valueOf(card3.substring(12,24)).doubleValue();
	force3=Double.valueOf(card4.substring(12,24)).doubleValue();
    };

    public BondDefinition(String card1, String card2){
	atomi=Integer.parseInt(card1.substring(5,10).trim());
	atomj=Integer.parseInt(card1.substring(10,15).trim());
	type=Integer.parseInt(card1.substring(15,20).trim());
	source3=Integer.parseInt(card1.substring(24,25).trim());
	source1=source3;
	source2=source3;
	if(card1.charAt(22)!=' ') source1=Integer.parseInt(card1.substring(20,23).trim());
	if(card1.charAt(23)!=' ') source2=Integer.parseInt(card1.substring(23,24).trim());
	source=source1;
	bond1=Double.valueOf(card2.substring(0,12)).doubleValue();
	bond2=bond1;
	bond3=bond1;
	force1=Double.valueOf(card2.substring(12,24)).doubleValue();
	force2=force1;
	force3=force1;
    };

    public BondDefinition(String card){
	atomi=Integer.parseInt(card.substring(0,5).trim());
	atomj=Integer.parseInt(card.substring(5,10).trim());
    };

}
