import java.io.*;
import javax.swing.*;

public class nwchem_Filter extends javax.swing.filechooser.FileFilter {
  public boolean accept(File name){
    boolean correct = false;
    if(name.toString().endsWith(".nw")) correct=true;
    if(name.toString().endsWith(".top")) correct=true;
    if(name.toString().endsWith(".frg")) correct=true;
    if(name.toString().endsWith(".sgm")) correct=true;
    if(name.toString().endsWith(".seq")) correct=true;
    if(name.toString().endsWith(".rst")) correct=true;
    if(name.toString().endsWith(".pdb")) correct=true;
    if(name.toString().endsWith(".rdf")) correct=true;
    if(name.toString().endsWith(".rdi")) correct=true;
    if(name.toString().endsWith(".coo")) correct=true;
    if(name.toString().endsWith(".sco")) correct=true;
    if(name.toString().endsWith(".vel")) correct=true;
    if(name.toString().endsWith(".svl")) correct=true;
    if(name.toString().endsWith(".prp")) correct=true;
    if(name.toString().endsWith(".qrs")) correct=true;
    if(name.toString().endsWith(".gib")) correct=true;
    if(name.toString().endsWith(".cnv")) correct=true;
    if(name.toString().endsWith(".fet")) correct=true;
    if(name.toString().endsWith(".acf")) correct=true;
    if(!name.isFile()) correct=true;
    return correct;
  }

  public String getDescription() {
    String type = new String("*.nw, *.top, *.frg, *.sgm, *.seq, *.rst, " +
			     "\n*.pdb, *.rdf, *.rdi, *.coo, *.sco, *.vel, " +
			     "\n*.svl, *.prp, *.qrs, *.gib, *.cnv, *.fet, " +
			     "\n*acf");
    return type; }
  
}
