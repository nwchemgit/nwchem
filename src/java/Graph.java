// $Id: Graph.java,v 1.1 1999-08-10 14:34:40 d3j191 Exp $
import ptplot.*;
import java.awt.*;
import java.util.*;
import java.awt.event.*;

public class Graph extends Plot{

  public boolean first = true, refill = true;
  public double xmin, xmax, ymin, ymax, dx, dy;

  public void addData(int set, double xval, double yval, boolean connect){
    if(first){
      xmin=xval;
      xmax=xval;
      ymin=yval;
      ymax=yval;
    } else {
      dx=xmax-xmin;
      dy=ymax-ymin;
      if(xval<xmin) { xmin=xval-0.2*dx; refill=true; };
      if(xval>xmax) { xmax=xval+0.2*dx; refill=true; };
      if(yval<ymin) { ymin=yval-0.2*dy; refill=true; };
      if(yval>ymax) { ymax=yval+0.2*dy; refill=true; };
    };
    addPoint(set,xval,yval,connect);
    if(refill) {
      setXRange(xmin,xmax);
      setYRange(ymin,ymax);
      drawPlot(getGraphics(),true);
    };
    first=false;
    refill=false;
  }

  public void removeSet(int dataset){
    Vector pts = _points[dataset];
    pts.removeAllElements();
    _xTop = - Double.MAX_VALUE;
    _yTop = - Double.MAX_VALUE;
    _xBottom = Double.MAX_VALUE;
    _yBottom = Double.MAX_VALUE;
  }

}
