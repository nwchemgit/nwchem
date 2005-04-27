import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

public class IntegerField extends JTextField implements KeyListener {

  int var;
  int def;

  public IntegerField(int value, int deflt, int width){
    var=value;
    def=deflt;
    setText(Integer.toString(var));
    setColumns(width);
    addKeyListener(this);
  }

  public int getValue(){
    if(isValid()){ return var; } else { return def; }
  }

  public boolean isValid(){
    try{
      Integer.parseInt(getText());
    } catch (NumberFormatException e) {
      var=def;
      setText(Integer.toString(var));
      return false;
    };
    return true;
  }

  public void keyPressed(KeyEvent e){}

  public void keyReleased(KeyEvent e){
    if(!isValid()){
      selectAll();
      requestFocus();
    };
    var=Integer.parseInt(getText());
  }

  public void keyTyped(KeyEvent e) {};

}
