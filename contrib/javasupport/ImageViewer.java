import java.awt.*;
import java.net.*;

public class ImageViewer extends Canvas {

Image image;

public ImageViewer(String filename) {


  image = ImageLoader.getImage(this,"",filename);

}
public void paint(Graphics g) {
	g.drawImage(image, 0, 0, this);
}


}
