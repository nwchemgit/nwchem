import java.io.*;

public class AngleDefinition{

    public int atomi,atomj,atomk;
    public int type, source;
    public double angle1,angle2,angle3,force1,force2,force3;

    public AngleDefinition(String card1, String card2){
	atomi=Integer.parseInt(card1.substring(5,10).trim());
	atomj=Integer.parseInt(card1.substring(10,15).trim());
	atomk=Integer.parseInt(card1.substring(15,20).trim());
	type=Integer.parseInt(card1.substring(20,25).trim());
	source=Integer.parseInt(card1.substring(25,30).trim());
	angle1=Double.valueOf(card2.substring(0,10)).doubleValue();
	angle2=Double.valueOf(card2.substring(22,32)).doubleValue();
	angle3=Double.valueOf(card2.substring(44,54)).doubleValue();
	force1=Double.valueOf(card2.substring(10,22)).doubleValue();
	force2=Double.valueOf(card2.substring(32,44)).doubleValue();
	force3=Double.valueOf(card2.substring(54,66)).doubleValue();
    }

}
