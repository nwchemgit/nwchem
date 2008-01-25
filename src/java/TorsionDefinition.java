import java.io.*;

public class TorsionDefinition{

    public int atomi,atomj,atomk,atoml;
    public int type, source;
    public int multi1,multi2,multi3;
    public int source1,source2,source3;
    public double torsion1,torsion2,torsion3,force1,force2,force3;

    public TorsionDefinition(String card1, String card2, String card3, String card4){
	atomi=Integer.parseInt(card1.substring(5,10).trim());
	atomj=Integer.parseInt(card1.substring(10,15).trim());
	atomk=Integer.parseInt(card1.substring(15,20).trim());
	atoml=Integer.parseInt(card1.substring(20,25).trim());
	type=Integer.parseInt(card1.substring(25,30).trim());
	source3=Integer.parseInt(card1.substring(34,35).trim());
	source1=source3;
	source2=source3;
	if(card1.charAt(32)!=' ') source1=Integer.parseInt(card1.substring(30,33).trim());
	if(card1.charAt(33)!=' ') source2=Integer.parseInt(card1.substring(33,34).trim());
	source=source1;
	multi1=Integer.parseInt(card2.substring(0,3).trim());
	multi2=Integer.parseInt(card3.substring(0,3).trim());
	multi3=Integer.parseInt(card4.substring(0,3).trim());
	torsion1=Double.valueOf(card2.substring(3,13)).doubleValue();
	torsion2=Double.valueOf(card3.substring(3,13)).doubleValue();
	torsion3=Double.valueOf(card4.substring(3,13)).doubleValue();
	force1=Double.valueOf(card2.substring(13,25)).doubleValue();
	force2=Double.valueOf(card3.substring(13,25)).doubleValue();
	force3=Double.valueOf(card4.substring(13,25)).doubleValue();
    }

    public TorsionDefinition(String card1, String card2){
	atomi=Integer.parseInt(card1.substring(5,10).trim());
	atomj=Integer.parseInt(card1.substring(10,15).trim());
	atomk=Integer.parseInt(card1.substring(15,20).trim());
	atoml=Integer.parseInt(card1.substring(20,25).trim());
	type=Integer.parseInt(card1.substring(25,30).trim());
	source3=Integer.parseInt(card1.substring(34,35).trim());
	source1=source3;
	source2=source3;
	if(card1.charAt(32)!=' ') source1=Integer.parseInt(card1.substring(30,33).trim());
	if(card1.charAt(33)!=' ') source2=Integer.parseInt(card1.substring(33,34).trim());
	source=source1;
	multi1=Integer.parseInt(card2.substring(0,3).trim());
	multi2=multi1;
	multi3=multi1;
	torsion1=Double.valueOf(card2.substring(3,13)).doubleValue();
	torsion2=torsion1;
	torsion3=torsion1;
	force1=Double.valueOf(card2.substring(13,25)).doubleValue();
	force2=force1;
	force3=force1;
    }

}
