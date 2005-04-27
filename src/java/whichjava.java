// $Id: whichjava.java,v 1.4 2005-04-27 19:02:26 marat Exp $

public class whichjava implements Runnable {
  public static void main(String[] args){
    System.out.println("The Java Version in your path is " + System.getProperty("java.version"));
    if (System.getProperty("java.version").indexOf("1.2")>=0){System.exit((int)0);};
    if (System.getProperty("java.version").indexOf("1.3")>=0){System.exit((int)0);};
    if (System.getProperty("java.version").indexOf("1.4")>=0){System.exit((int)0);};
    System.out.println (" This code requires Java 1.2 or greater"); System.exit((int)1);
  }

  public void init(){};

  public void start(){};

  public void run(){};

}

