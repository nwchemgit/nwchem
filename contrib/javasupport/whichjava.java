// $Id$

public class whichjava implements Runnable {
  public static void main(String[] args){
    System.out.println("The Java Version in your path is " + System.getProperty("java.version"));
    if (System.getProperty("java.version").indexOf("1.1")>=0){
      //      System.out.println(" should be zero status" +  System.getProperty("java.version").substring(0,3));
      System.exit((int)0);
    } else {
      //      System.out.println(" should be nonzero status " +  System.getProperty("java.version").substring(0,3));
      System.out.println (" This code requires Java 1.1 or greater");
      System.exit((int)1);
    };
  }

  public void init(){};

  public void start(){};

  public void run(){};

}

