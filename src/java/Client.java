import java.io.*;
import java.net.*;

public class Client{
  private static final int PORTNUM = 3333;

  public static void main(String[] arguments){
    Socket socket = null;
    try{
      socket = new Socket("127.0.0.1",PORTNUM);
      socket.close();
    } catch (IOException e) {e.printStackTrace();};
  }
}
