/**
 *
 * ImageLoader Class
 * <p>
 * Offers a way to load an image (JPEG or GIF) from a disk location relative
 * the a Java Class. 
 * These methods work when using the code as a StandAlone application, class files 
 * inside a Browser, and or and JAR file loaded into a Browser.
 *
 * To use, simply call the static method with something like this...
 * Image myImage = ImageLoader.getImage(myComponent, "images", "filename.jpg");
 * where myComponent is the component to display the image on.
 *
 * @version 1.0 May 1998
 * @author Chris Parkinson
 */

import java.awt.*;
import java.awt.image.*;
import java.util.*;
import java.io.*;
import java.net.URL;



public class ImageLoader 
{

 
/////////////////////////////////////////////////////////////////
// 
// Constructors
//
/////////////////////////////////////////////////////////////////
/**
 * Construct nothing !
 */
 public ImageLoader() {
 }


////////////////////////////////////////////////////////
//
// Get a Resource Image
//
////////////////////////////////////////////////////////
/**
 * Get a Resource Image.
 * Read in an image file from a resource ( a file located in the class path)
 * This method accesses images using the getResourceStream methods of the 
 * Class class.	
 * <p>
 * @param owner the object who wants the image icon.
 * @param pathname the pathname for the image .
 * @param filename the filename for the image (GIF or JPG format).
 * @return an Image
 */
 public static Image getImage(Component owner, String pathname, String filename)
{
 
   	byte imageBytes[] = getImageBytes(owner.getClass(), pathname, filename);
  	if (imageBytes!=null) return getImage(owner, imageBytes);
	return null;
}
	

/**
 * Get a the byte array of a Resource Image.
 * Read in an image file from a resource ( a file located in the class path)
 * and return the image byte array.
 * This method accesses images using the getResourceStream methods of the 
 * Class class.	The file to access must be in the same directory as the class owner.
 * <p>
 * @param classOwner the class who wants the image icon.
 * @param pathname the pathname for the image .
 * @param filename the filename for the image (GIF or JPG format).
 * @return an array of bytes.
 */
public static byte[] getImageBytes(Class classOwner,  String pathname, String
filename) {
	String fullFile = filename;
    if (pathname!=null)fullFile  = pathname + "/" + filename;

	//Try and load our images the correct way 	
	InputStream inputStream = classOwner.getResourceAsStream(fullFile);

	//Now we will load whatever we have
    try {
	   	byte array[] = new byte[inputStream.available()];
	   	inputStream.read(array);
	   	inputStream.close();
	  	return array;
	}catch(Exception e) {
		System.out.println("Error occurred trying to read a resource image ");
		System.out.println("Filename = " + fullFile);
		System.out.println("Class = " + classOwner);
		System.out.println(e.toString());
		return null;
	}
}

/**
 * Get an image from an array of bytes.
 * <p>
 * @param owner a component who owns the image
 * @param imageData an array of bytes containing an image
 * @return an Image
 */
public static Image getImage( Component owner, byte[] imageData) {
  	try {
  		Toolkit toolkit = Toolkit.getDefaultToolkit();
		Image image = toolkit.createImage(imageData);
		MediaTracker tracker = new MediaTracker(owner);
		tracker.addImage(image,0);
		tracker.waitForAll();
		return image;
	}catch(Exception e)  {}
	return null;
}

									  
		 
	 
}

