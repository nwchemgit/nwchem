import javax.swing.filechooser.FileFilter;
import java.io.File;

public class ExtensionFilter extends FileFilter
{

  private String extension;

  public ExtensionFilter(String ext)
  {
    extension=ext.toLowerCase();
    System.out.println("Extension Filter for "+extension);
  }
  public boolean accept(File file)
  {
    return (file.isDirectory() || file.getName().toLowerCase().endsWith(extension));
  }
  public String getDescription()
  {
    return " ";
  }
}
