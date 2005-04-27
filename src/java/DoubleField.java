import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

public class DoubleField extends JTextField implements KeyListener{

  double var;
  double def;

  public DoubleField(double value, double deflt, int width){
    var=value;
    def=deflt;
    setText(Double.toString(var));
    setColumns(width);
    addKeyListener(this);
  }

  public double getValue(){
    if(isValid()){ return var; } else { return def; }
  }

  public boolean isValid(){
    try{
      Double.valueOf(getText()).doubleValue();
    } catch (NumberFormatException e) {
      var=def;
      setText(Double.toString(var));
      return false;
    };
    return true;
  }
  public void keyPressed(KeyEvent e) {}

  public void keyReleased(KeyEvent e) {
    if(!isValid()){
      selectAll();
      requestFocus();
    };
    var=Double.valueOf(getText()).doubleValue();
  }

  public void keyTyped(KeyEvent e) {}

}

