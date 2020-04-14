#!/usr/local/bin/python -t
# $Id$

# Gnuplot.py -- A pipe-based interface to the gnuplot plotting program.

# Copyright (C) 1998 Michael Haggerty <mhagger@blizzard.harvard.edu>.

"""A pipe-based interface to the gnuplot plotting program.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.  This program is distributed in the
hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the GNU General Public License for more details; it is
available at <http://www.fsf.org/copyleft/gpl.html>, or by writing to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
Boston, MA 02111-1307, USA.

Written by Michael Haggerty <mhagger@blizzard.harvard.edu>.  Inspired
by and partly derived from an earlier version by Konrad Hinsen
<hinsen@ibs.ibs.fr>.  If you find a problem or have a suggestion,
please let me know at <mhagger@blizzard.harvard.edu>.  Other feedback
is also welcome.

For information about how to use this module, see the comments below,
the documentation string for class Gnuplot, and the test code at the
bottom of the file.  You can run the test code by typing
'python Gnuplot.py'.

You should import this file with 'import Gnuplot', not with
'from Gnuplot import *'; otherwise you will have problems with
conflicting names (specifically, the Gnuplot module name conflicts
with the Gnuplot class name).  To obtain gnuplot itself, see
<http://www.cs.dartmouth.edu/gnuplot_info.html>.

Features:

 o  Allows the creation of two or three dimensional plots from
    python by piping commands to the 'gnuplot' program.
 o  A gnuplot session is an instance of class 'Gnuplot', so multiple
    sessions can be open at once:
        'g1 = Gnuplot.Gnuplot(); g2 = Gnuplot.Gnuplot()'
 o  The implicitly-generated gnuplot commands can be stored to a file
    instead of executed immediately:
        'g = Gnuplot.Gnuplot("commands.gnuplot")'
    The file can then be run later with gnuplot's 'load' command.
    Beware, however, if the plot commands depend on the existence of
    temporary files, because they might be deleted before you use
    the command file.
 o  Can pass arbitrary commands to the gnuplot command interpreter:
        'g("set pointsize 2")'
 o  A Gnuplot object knows how to plot objects of type 'PlotItem'.
    Any PlotItem can have optional `title' and/or 'with' suboptions.
    Builtin PlotItem types:
        
    * 'Data(array1)' -- data from a Python list or NumPy array
                      (permits additional option 'cols' )
    * 'File("filename")' -- data from an existing data file (permits
                      additional option 'using' )
    * 'Func("exp(4.0 * sin(x))")' -- functions (passed as a string
                      for gnuplot to evaluate)
    * 'GridData(m, x, y)' -- data tabulated on a grid of (x,y) values
                      (usually to be plotted in 3-D)

    See those classes for more details.

 o PlotItems are implemented as objects that can be assigned to
    variables (including their options) and plotted repeatedly ---
    this also saves much of the overhead of plotting the same data
    multiple times.
 o  Communication of data between python and gnuplot is via temporary
    files, which are deleted automatically when their associated
    'PlotItem' is deleted.  (Communication of commands is via a pipe.)
    The PlotItems currently in use by a Gnuplot object are stored in
    an internal list so that they won't be deleted prematurely.
 o  Can use 'replot' method to add datasets to an existing plot.
 o  Can make persistent gnuplot windows by using the constructor option
    `persist=1'.  Such windows stay around even after the gnuplot
    program is exited.  Note that only newer version of gnuplot support
    this option.
 o  Plotting to a postscript file is via new 'hardcopy' method, which
    outputs the currently-displayed plot to either a postscript
    printer or to a postscript file.
 o  There is a 'plot' command which is roughly compatible with the
    command from Konrad Hinsen's old 'Gnuplot.py'.

Restrictions:
    
 -  Relies on the numpy Python extension.  This can be obtained
    from LLNL (See ftp://ftp-icf.llnl.gov/pub/python/README.html).
    If you're interested in gnuplot, you would probably also want
    NumPy anyway.
 -  Probably depends on a unix-type environment.  Anyone who wants
    to remedy this situation should contact me.
 -  Only a small fraction of gnuplot functionality is implemented as
    explicit method functions.  However, you can give arbitrary
    commands to gnuplot manually; for example:
        'g = Gnuplot.Gnuplot()',
        'g('set style data linespoints')',
        'g('set pointsize 5')',
    etc.  I might add a more organized way of setting arbitrary
    options.
 -  There is no provision for missing data points in array data
    (which gnuplot would allow by specifying '?' as a data point).
    I can't think of a clean way to implement this; maybe one could
    use NaN for machines that support IEEE floating point.
 -  There is no supported way to change the plotting options of
    PlotItems after they have been created.

Bugs:

 -  No attempt is made to check for errors reported by gnuplot (but
    they will appear on stderr).
 -  All of these classes perform their resource deallocation when
    '__del__' is called.  If you delete things explicitly, there will
    be no problem.  If you don't, an attempt is made to delete
    remaining objects when the interpreter is exited, but this is
    not completely reliable, so sometimes temporary files will be
    left around.  If anybody knows how to fix this problem, please
    let me know.
"""

__version__ = "1.1a"
__cvs_version__ = "CVS version $Revision: 1.1 $"

import sys, os, string, tempfile, numpy


# Set after first call of test_persist().  This will be set from None
# to 0 or 1 upon the first call of test_persist(), then the stored
# value will be used thereafter.  To avoid the test, type 1 or 0 on
# the following line corresponding to whether your gnuplot is new
# enough to understand the -persist option.
_recognizes_persist = None

# After a hardcopy is produced, we have to set the terminal type back
# to `on screen'.  If you are using unix, then `x11' is probably
# correct.  If not, change the following line to the terminal type you
# use.
_default_term = 'x11'

# Gnuplot can plot to a printer by using `set output "| ..."' where
# ... is the name of a program that sends its stdin to a printer.  On
# my machine the appropriate program is `lpr', as set below.  On your
# computer it may be something different (like `lp'); you can set that
# by changing the variable below.  You can also use the following
# variable to add options to the print command.
_default_lpr = '| lpr'

def test_persist():
    """Test and report whether gnuplot recognizes the option '-persist'.

    Test if gnuplot is new enough to know the option '-persist'.  It it
    isn't, it will emit an error message with '-persist' in the first
    line.

    """

    global _recognizes_persist
    if _recognizes_persist is None:
        g = os.popen('echo | gnuplot -persist 2>&1', 'r')
        response = g.readlines()
        g.close()
        _recognizes_persist = ((not response)
                               or (string.find(response[0], '-persist') == -1))
    return _recognizes_persist


class OptionException(Exception):
    """raised for unrecognized option(s)"""
    pass

class DataException(Exception):
    """raised for data in the wrong format"""
    pass


class PlotItem:
    """Plotitem represents an item that can be plotted by gnuplot.

    For the finest control over the output, you can create the
    PlotItems yourself with additional keyword options, or derive new
    classes from PlotItem.

    Members:
    
    'basecommand' -- a string holding the elementary argument that
                     must be passed to gnuplot's `plot' command for
                     this item; e.g., 'sin(x)' or '"filename.dat"'.
    'options' -- a list of strings that need to be passed as options
                 to the plot command, in the order required; e.g.,
                 ['title "data"', 'with linespoints'].
    'title' -- the title requested (undefined if not requested).  Note
               that `title=None' implies the `notitle' option,
               whereas omitting the title option implies no option
               (the gnuplot default is then used).
    'with' -- the string requested as a `with' option (undefined if
              not requested)

    """

    def __init__(self, basecommand, **keyw):
        self.basecommand = basecommand
        self.options = []
        if keyw.has_key('title'):
            self.title = keyw['title']
            del keyw['title']
            if self.title is None:
                self.options.append('notitle')
            else:
                self.options.append('title "' + self.title + '"')
        if keyw.has_key('with_'):
            self.with_ = keyw['with_']
            del keyw['with_']
            self.options.append('with ' + self.with_)
        if keyw:
            raise OptionException(keyw)

    def command(self):
        """Build the 'plot' command to be sent to gnuplot.

        Build and return the 'plot' command, with options, necessary
        to display this item.

        """

        if self.options:
            return self.basecommand + ' ' + string.join(self.options)
        else:
            return self.basecommand

    # if the plot command requires data to be put on stdin (i.e.,
    # `plot "-"'), this method should put that data there.
    def pipein(self, file):
        pass


class Func(PlotItem):
    """Represents a mathematical expression to plot.

    Func represents a mathematical expression that is to be computed by
    gnuplot itself, as in the example

        gnuplot> plot sin(x)

    The argument to the contructor is a string which is a expression.
    Example:

        g.plot(Func("sin(x)", with_="line 3"))

    or the shorthand example:

        g.plot("sin(x)")

    """

    def __init__(self, funcstring, **keyw):
        apply(PlotItem.__init__, (self, funcstring), keyw)


class AnyFile:
    """An AnyFile represents any kind of file to be used by gnuplot.

    An AnyFile represents a file, but presumably one that holds data
    in a format readable by gnuplot.  This class simply remembers the
    filename; the existence and format of the file are not checked
    whatsoever.  Note that this is not a PlotItem, though it is used by
    the 'File' PlotItem.  Members:

    'self.filename' -- the filename of the file

    """

    def __init__(self, filename):
        self.filename = filename


class TempFile(AnyFile):
    """A TempFile is a file that is automatically deleted.

    A TempFile points to a file.  The file is deleted automatically
    when the TempFile object is deleted.

    WARNING: whatever filename you pass to this constructor **WILL BE
    DELETED** when the TempFile object is deleted, even if it was a
    pre-existing file! This is intended to be used as a parent class of
    TempArrayFile.

    """

    def __del__(self):
        os.unlink(self.filename)


def write_array(f, set,
                item_sep=' ',
                nest_prefix='', nest_suffix='\n', nest_sep=''):
    """Write an array of arbitrary dimension to a file.

    A general recursive array writer.  The last four parameters allow a
    great deal of freedom in choosing the output format of the
    array.  The defaults for those parameters give output that is
    gnuplot-readable.  But using, for example, ( ',', '{', '}', ',\\n'
    ) would output an array in a format that Mathematica could
    read.  item_sep should not contain '%' (or if it does, it should be
    escaped to '%%' ) since item_sep is put into a format string.

    """

    if len(set.shape) == 1:
        (columns,) = set.shape
        assert columns > 0
        fmt = string.join(['%s'] * columns, item_sep)
        f.write(nest_prefix)
        f.write(fmt % tuple(set.tolist()))
        f.write(nest_suffix)
    elif len(set.shape) == 2:
        # This case could be done with recursion, but `unroll' for
        # efficiency.
        (points, columns) = set.shape
        assert points > 0
        assert columns > 0
        fmt = string.join(['%s'] * columns, item_sep)
        f.write(nest_prefix + nest_prefix)
        f.write(fmt % tuple(set[0].tolist()))
        f.write(nest_suffix)
        for point in set[1:]:
            f.write(nest_sep + nest_prefix)
            f.write(fmt % tuple(point.tolist()))
            f.write(nest_suffix)
        f.write(nest_suffix)
    else:
        assert set.shape[0] > 0
        f.write(nest_prefix)
        write_array(f, set[0], item_sep, nest_prefix, nest_suffix, nest_sep)
        for subset in set[1:]:
            f.write(nest_sep)
            write_array(f, subset, item_sep, nest_prefix, nest_suffix, nest_sep)
        f.write(nest_suffix)


class ArrayFile(AnyFile):
    """A file to which, upon creation, an array is written.

    When an ArrayFile is constructed, it creates a file and fills it
    with the contents of a 2-d or 3-d numpy array in the format
    expected by gnuplot.  Specifically, for 2-d, the file organization
    is for example:

        set[0,0] set[0,1] ...
        set[1,0] set[1,1] ...

    etc.  For 3-d, it is for example:

        set[0,0,0] set[0,0,1] ...
        set[0,1,0] set[0,1,1] ...

        set[1,0,0] set[1,0,1] ...
        set[1,1,0] set[1,1,1] ...

    etc.

    The filename can be specified, otherwise a random filename is
    chosen.  The file is NOT deleted automatically.

    """

    def __init__(self, set, filename=None):
        if not filename:
            filename = tempfile.mktemp()
        f = open(filename, 'w')
        write_array(f, set)
        f.close()
        AnyFile.__init__(self, filename)


class TempArrayFile(ArrayFile, TempFile):
    """An ArrayFile that is deleted automatically."""

    def __init__(self, set, filename=None):
        ArrayFile.__init__(self, set, filename)


class File(PlotItem):
    """A PlotItem representing a file that contains gnuplot data.

    File is a PlotItem that represents a file that should be plotted
    by gnuplot.  The file can either be a string holding the filename
    of an existing file, or it can be anything derived from 'AnyFile'.

    """

    def __init__(self, file, using=None, **keyw):
        """Construct a File object.

        '<file>' can be either a string holding the filename of an
        existing file, or it can be an object of a class derived from
        'AnyFile' (such as a 'TempArrayFile').  Keyword arguments
        recognized (in addition to those recognized by 'PlotItem' ):

            'using=<n>' -- plot that column against line number
            'using=<tuple>' -- plot using a:b:c:d etc.
            'using=<string>' -- plot `using <string>' (allows gnuplot's
                              arbitrary column arithmetic)

        Note that the 'using' option is interpreted by gnuplot, so
        columns must be numbered starting with 1.  Other keyword
        arguments are passed along to PlotItem.  The default 'title'
        for an AnyFile PlotItem is 'notitle'.

        """

        if isinstance(file, AnyFile):
            self.file = file
            # If no title is specified, then use `notitle' for
            # TempFiles (to avoid using the temporary filename as the
            # title.)
            if isinstance(file, TempFile) and not keyw.has_key('title'):
                keyw['title'] = None
        elif type(file) == type(""):
            self.file = AnyFile(file)
        else:
            raise OptionException
        apply(PlotItem.__init__, (self, '"' + self.file.filename + '"'), keyw)
        self.using = using
        if self.using is None:
            pass
        elif type(self.using) == type(""):
            self.options.insert(0, "using " + self.using)
        elif type(self.using) == type(()):
            self.options.insert(0,
                                "using " +
                                string.join(map(repr, self.using), ':'))
        elif type(self.using) == type(1):
            self.options.insert(0, "using " + repr(self.using))
        else:
            raise OptionException('using=' + repr(self.using))


class Data(File):
    """Allows data from memory to be plotted with Gnuplot.

    Takes a numeric array from memory and outputs it to a temporary
    file that can be plotted by gnuplot.

    """

    def __init__(self, *set, **keyw):
        """Construct a Data object from a numeric array.

        Create a Data object (which is a type of PlotItem) out of one
        or more Float Python numpy arrays (or objects that can be
        converted to a Float numpy array).  If the routine is passed
        one array, the last index ranges over the values comprising a
        single data point (e.g., [x, y, and sigma]) and the rest of
        the indices select the data point.  If the routine is passed
        more than one array, they must have identical shapes, and then
        each data point is composed of one point from each array.
        I.e., 'Data(x,x**2)' is a PlotItem that represents x squared
        as a function of x.  For the output format, see the comments
        in ArrayFile.

        The array is first written to a temporary file, then that file
        is plotted.  Keyword arguments recognized (in addition to those
        recognized by PlotItem):

            cols=<tuple> -- write only the specified columns from each
                            data point to the file.  Since cols is
                            used by python, the columns should be
                            numbered in the python style (starting
                            from 0), not the gnuplot style (starting
                            from 1).

        The data are immediately written to the temp file; no copy is
        kept in memory.

        """

        if len(set) == 1:
            # set was passed as a single structure
            set = numpy.asarray(set, dtype=numpy.float32)
        else:
            # set was passed column by column (for example, Data(x,y))
            set = numpy.asarray(set, dtype=numpy.float32)
            dims = len(set.shape)
            # transpose so that the last index selects x vs. y:
            set = numpy.transpose(set, (dims-1,) + tuple(range(dims-1)))
        if keyw.has_key('cols') and keyw['cols'] is not None:
            set = numpy.take(set, keyw['cols'], -1)
            del keyw['cols']
        apply(File.__init__, (self, TempArrayFile(set)), keyw)


class GridData(File):
    """Holds data representing a function of two variables, for use in splot.

    GridData represents a function that has been tabulated on a
    rectangular grid.  It is a PlotItem, so GridData objects can be
    plotted by Gnuplot.  The data are written to a file but not stored
    in memory.

    """

    def __init__(self, data, xvals=None, yvals=None, **keyw):
        """GridData constructor.

        Arguments:

            'data' -- a 2-d array with dimensions (numx,numy)
            'xvals' -- a 1-d array with dimension (numx)
            'yvals' -- a 1-d array with dimension (numy)

        'data' is meant to hold the values of a function f(x,y) tabulated
        on a grid of points, such that 'data[i,j] == f(xvals[i],
        yvals[j])'.  These data are written to a datafile as 'x y f(x,y)'
        triplets that can be used by gnuplot's splot command.  Thus if you
        have three arrays in the above format and a Gnuplot instance
        called g, you can plot your data by typing for example:

            g.splot(Gnuplot.GridData(data,xvals,yvals))

        If 'xvals' and/or 'yvals' are omitted, integers (starting with
        0) are used for that coordinate.  The data are written to a
        temporary file; no copy of the data is kept in memory.

        """

        data = numpy.asarray(data, dtype=numpy.float32)
        assert len(data.shape) == 2
        (numx, numy) = data.shape

        if xvals is None:
            xvals = numpy.arange(numx)
        else:
            xvals = numpy.asarray(xvals, dtype=numpy.float32)
            assert len(xvals.shape) == 1
            assert xvals.shape[0] == numx

        if yvals is None:
            yvals = numpy.arange(numy)
        else:
            yvals = numpy.asarray(yvals, dtype=numpy.float32)
            assert len(yvals.shape) == 1
            assert yvals.shape[0] == numy

        set = numpy.transpose(
            numpy.array(
                (numpy.transpose(numpy.resize(xvals, (numy, numx))),
                 numpy.resize(yvals, (numx, numy)),
                 data)), (1,2,0))

        apply(File.__init__, (self, TempArrayFile(set)), keyw)


def grid_function(f, xvals, yvals):
    """Compute a function on a grid.

    'xvals' and 'yvals' should be 1-D arrays listing the values of x
    and y at which f should be tabulated.  f should be a function
    taking two floating point arguments.  The return value is a matrix
    M where M[i,j] = f(xvals[i],yvals[j]), which can for example be
    used in the 'GridData' constructor.

    Note that f is evaluated at each pair of points using a Python loop,
    which can be slow if the number of points is large.  If speed is an
    issue, you are better off computing functions matrix-wise using
    numpy's built-in ufuncs.

    """

    m = numpy.zeros((len(xvals), len(yvals)), dtype=numpy.float32)
    for xi in range(len(xvals)):
        x = xvals[xi]
        for yi in range(len(yvals)):
            y = yvals[yi]
            m[xi,yi] = f(x,y)
    return m


class Gnuplot:
    """gnuplot plotting object.

    A Gnuplot represents a running gnuplot program and a pipe to
    communicate with it.  It keeps a reference to each of the
    PlotItems used in the current plot, so that they (and their
    associated temporary files) are not deleted prematurely.  The
    communication is one-way; gnuplot's text output just goes to
    stdout with no attempt to check it for error messages.

    Members:

    'gnuplot' -- the pipe to gnuplot or a file gathering the commands
    'itemlist' -- a list of the PlotItems that are associated with the
                  current plot.  These are deleted whenever a new plot
                  command is issued via the `plot' method.
    'debug' -- if this flag is set, commands sent to gnuplot will also
               be echoed to stderr.
    'plotcmd' -- 'plot' or 'splot', depending on what was the last
                 plot command.

    Methods:

    '__init__' -- if a filename argument is specified, the commands
                  will be written to that file instead of being piped
                  to gnuplot immediately.
    'plot' -- clear the old plot and old PlotItems, then plot the
              arguments in a fresh plot command.  Arguments can be: a
              PlotItem, which is plotted along with its internal
              options; a string, which is plotted as a Func; or
              anything else, which is plotted as a Data.
    'hardcopy' -- replot the plot to a postscript file (if filename
                  argument is specified) or pipe it to lpr othewise.
                  If the option `color' is set to true, then output
                  color postscript.
    'replot' -- replot the old items, adding any arguments as
                additional items as in the plot method.
    'refresh' -- issue (or reissue) the plot command using the current
                 PlotItems.
    '__call__' -- pass an arbitrary string to the gnuplot process,
                  followed by a newline.
    'xlabel', 'ylabel', 'title' -- set attribute to be a string.
    'interact' -- read lines from stdin and send them, one by one, to
                  the gnuplot interpreter.  Basically you can type
                  commands directly to the gnuplot command processor
                  (though without command-line editing).
    'load' -- load a file (using the gnuplot `load' command).
    'save' -- save gnuplot commands to a file (using gnuplot `save'
              command) If any of the PlotItems is a temporary file, it
              will be deleted at the usual time and the save file might
              be pretty useless :-).
    'clear' -- clear the plot window (but not the itemlist).
    'reset' -- reset all gnuplot settings to their defaults and clear
               the current itemlist.
    'set_string' -- set or unset a gnuplot option whose value is a
                    string.
    '_clear_queue' -- clear the current PlotItem list.
    '_add_to_queue' -- add the specified items to the current
                       PlotItem list.

    """

    def __init__(self, filename=None, persist=0, debug=0):
        """Create a Gnuplot object.

        'Gnuplot(filename=None, persist=0, debug=0)':

        Create a 'Gnuplot' object.  By default, this starts a gnuplot
        process and prepares to write commands to it.  If 'filename'
        is specified, the commands are instead written to that file
        (i.e., for later use using 'load').  If 'persist' is set,
        gnuplot will be started with the '-persist' option (which
        creates a new X11 plot window for each plot command).  (This
        option is not available on older versions of gnuplot.)  If
        'debug' is set, the gnuplot commands are echoed to stderr as
        well as being send to gnuplot.

        """

        if filename:
            # put gnuplot commands into a file:
            self.gnuplot = open(filename, 'w')
        else:
            if persist:
                if not test_persist():
                    raise OptionException(
                        '-persist does not seem to be supported '
                        'by your version of gnuplot!')
                self.gnuplot = os.popen('gnuplot -persist', 'w')
            else:
                self.gnuplot = os.popen('gnuplot', 'w')
        self._clear_queue()
        self.debug = debug
        self.plotcmd = 'plot'

    def __del__(self):
        self('quit')
        self.gnuplot.close()

    def __call__(self, s):
        """Send a command string to gnuplot.

        '__call__(s)': send the string s as a command to gnuplot,
        followed by a newline and flush.  All interaction with the
        gnuplot process is through this method.

        """

        self.gnuplot.write(s + "\n")
        self.gnuplot.flush()
        if self.debug:
            # also echo to stderr for user to see:
            sys.stderr.write("gnuplot> %s\n" % (s,))

    def refresh(self):
        """Refresh the plot, using the current PlotItems.

        Refresh the current plot by reissuing the gnuplot plot command
        corresponding to the current itemlist.

        """

        plotcmds = []
        for item in self.itemlist:
            plotcmds.append(item.command())
        self(self.plotcmd + ' ' + string.join(plotcmds, ', '))
        for item in self.itemlist:
            item.pipein(self.gnuplot)

    def _clear_queue(self):
        """Clear the PlotItems from the queue."""

        self.itemlist = []

    def _add_to_queue(self, items):
        """Add a list of items to the itemlist, but don't plot them.

        'items' is a list or tuple of items, each of which should be a
        'PlotItem' of some kind, a string (interpreted as a function
        string for gnuplot to evaluate), or a numpy array (or
        something that can be converted to a numpy array).

        """

        for item in items:
            if isinstance(item, PlotItem):
                self.itemlist.append(item)
            elif type(item) is type(""):
                self.itemlist.append(Func(item))
            else:
                # assume data is an array:
                self.itemlist.append(Data(item))

    def plot(self, *items):
        """Draw a new plot.

        'plot(item, ...)': Clear the current plot and create a new 2-d
        plot containing the specified items.  Arguments can be of the
        following types:

        'PlotItem' (e.g., 'Data', 'File', 'Func', 'GridData') -- This
                   is the most flexible way to call plot because the
                   PlotItems can contain suboptions.  Moreover,
                   PlotItems can be saved to variables so that their
                   lifetime is longer than one plot command--thus they
                   can be replotted with minimal overhead.

        'string' (i.e., "sin(x)") -- The string is interpreted as
                 'Func(string)' (a function that is computed by
                 gnuplot).

        Anything else -- The object, which should be convertible to an
                         array, is converted to a Data() item, and
                         thus plotted as data.  If the conversion
                         fails, an exception is raised.

        """

        # remove old files:
        self.plotcmd = 'plot'
        self._clear_queue()
        self._add_to_queue(items)
        self.refresh()

    def splot(self, *items):
        """Draw a new three-dimensional plot.

        'splot(item, ...)' -- Clear the current plot and create a new
                3-d plot containing the specified items.  Arguments can
                be of the following types:
        'PlotItem' (e.g., 'Data', 'File', 'Func', 'GridData' ) -- This
                is the most flexible way to call plot because the
                PlotItems can contain suboptions.  Moreover, PlotItems
                can be saved to variables so that their lifetime is
                longer than one plot command--thus they can be
                replotted with minimal overhead.

        'string' (i.e., "sin(x*y)") -- The string is interpreted as a
                'Func()' (a function that is computed by gnuplot).

        Anything else -- The object is converted to a Data() item, and
                thus plotted as data.  Note that each data point
                should normally have at least three values associated
                with it (i.e., x, y, and z).  If the conversion fails,
                an exception is raised.

        """

        # remove old files:
        self.plotcmd = 'splot'
        self._clear_queue()
        self._add_to_queue(items)
        self.refresh()

    def replot(self, *items):
        """Replot the data, possibly adding new PlotItems.

        Replot the existing graph, using the items in the current
        itemlist.  If arguments are specified, they are interpreted as
        additional items to be plotted alongside the existing items on
        the same graph.  See 'plot' for details.

        """

        self._add_to_queue(items)
        self.refresh()

    def interact(self):
        """Allow user to type arbitrary commands to gnuplot.

        Read stdin, line by line, and send each line as a command to
        gnuplot.  End by typing C-d.

        """

        sys.stderr.write("Press C-d to end interactive input\n")
        while 1:
            sys.stderr.write("gnuplot>>> ")
            line = sys.stdin.readline()
            if not line:
                break
            if line[-1] == "\n": line = line[:-1]
            self(line)

    def clear(self):
        """Clear the plot window (without affecting the current itemlist)."""

        self('clear')

    def reset(self):
        """Reset all gnuplot settings to their defaults and clear itemlist."""

        self('reset')
        self.itemlist = []

    def load(self, filename):
        """Load a file using gnuplot's `load' command."""

        self('load "%s"' % (filename,))

    def save(self, filename):
        """Save the current plot commands using gnuplot's `save' command."""

        self('save "%s"' % (filename,))

    def set_string(self, option, s=None):
        """Set a string option, or if s is omitted, unset the option."""

        if s is None:
            self('set %s' % (option,))
        else:
            self('set %s "%s"' % (option, s))

    def xlabel(self, s=None):
        """Set the plot's xlabel."""

        self.set_string('xlabel', s)

    def ylabel(self, s=None):
        """Set the plot's ylabel."""

        self.set_string('ylabel', s)

    def title(self, s=None):
        """Set the plot's title."""

        self.set_string('title', s)

    def hardcopy(self, filename=None, eps=0, color=0, enhanced=1):
        """Create a hardcopy of the current plot.

        Create a postscript hardcopy of the current plot.  If a
        filename is specified, save the output in that file; otherwise
        print it immediately using lpr.  If eps is specified, generate
        encapsulated postscript.  If color is specified, create a
        color plot.  If enhanced is specified (the default), then
        generate enhanced postscript.  (Some old gnuplot versions do
        not support enhanced postscript; if this is the case set
        enhanced=0.)  Note that this command will return immediately
        even though it might take gnuplot a while to actually finish
        working.

        """

        if filename is None:
            filename = _default_lpr
        setterm = ['set', 'term', 'postscript']
        if eps: setterm.append('eps')
        else: setterm.append(' ')
        if enhanced: setterm.append('enhanced')
        if color: setterm.append('color')
        self(string.join(setterm))
        self.set_string('output', filename)
        self.refresh()
        self('set term %s' % _default_term)
        self.set_string('output')


# The following is a command defined for compatibility with Hinson's
# old Gnuplot.py module.  Its use is deprecated.

# When the plot command is called and persist is not available, the
# plotters will be stored here to prevent their being closed:
_gnuplot_processes = []

def plot(*items, **keyw):
    """plot data using gnuplot through Gnuplot.

    This command is roughly compatible with old Gnuplot plot command.
    It is provided for backwards compatibility with the old functional
    interface only.  It is recommended that you use the new
    object-oriented Gnuplot interface, which is much more flexible.

    It can only plot numpy array data.  In this routine an NxM array
    is plotted as M-1 separate datasets, using columns 1:2, 1:3, ...,
    1:M.

    Limitations:

        - If persist is not available, the temporary files are not
          deleted until final python cleanup.

    """

    newitems = []
    for item in items:
        # assume data is an array:
        item = numpy.asarray(item, dtype=numpy.float32)
        dim = len(item.shape)
        if dim == 1:
            newitems.append(Data(item[:, numpy.newaxis], with_='lines'))
        elif dim == 2:
            if item.shape[1] == 1:
                # one column; just store one item for tempfile:
                newitems.append(Data(item, with_='lines'))
            else:
                # more than one column; store item for each 1:2, 1:3, etc.
                tempf = TempArrayFile(item)
                for col in range(1, item.shape[1]):
                    newitems.append(File(tempf, using=(1,col+1), with_='lines'))
        else:
            raise DataException("Data array must be 1 or 2 dimensional")
    items = tuple(newitems)
    del newitems

    if keyw.has_key('file'):
        g = Gnuplot()
        # setup plot without actually plotting (so data don't appear
        # on the screen):
        g._add_to_queue(items)
        g.hardcopy(keyw['file'])
        # process will be closed automatically
    elif test_persist():
        g = Gnuplot(persist=1)
        apply(g.plot, items)
        # process will be closed automatically
    else:
        g = Gnuplot()
        apply(g.plot, items)
        # prevent process from being deleted:
        _gnuplot_processes.append(g)


# Demo code
if __name__ == '__main__':
    from numpy import *
    import sys

    # A straightforward use of gnuplot.  The `debug=1' switch is used
    # in these examples so that the commands that are sent to gnuplot
    # are also output on stderr.
    g1 = Gnuplot(debug=1)
    g1.title('A simple example') # (optional)
    g1('set style data linespoints') # give gnuplot an arbitrary command
    # Plot a list of (x, y) pairs (tuples or a numpy array would
    # also be OK):
    g1.plot([[0.,1.1], [1.,5.8], [2.,3.3], [3.,4.2]])

    # Plot one dataset from an array and one via a gnuplot function;
    # also demonstrate the use of item-specific options:
    g2 = Gnuplot(debug=1)
    x = arange(10, dtype=numpy.float32)
    y1 = x**2
    # Notice how this plotitem is created here but used later?  This
    # is convenient if the same dataset has to be plotted multiple
    # times, because the data need only be written to a temporary file
    # once.
    d = Data(x, y1,
             title="calculated by python",
             with_="points pointsize 3 pointtype 3")
    g2.title('Data can be computed by python or gnuplot')
    g2.xlabel('x')
    g2.ylabel('x squared')
    # Plot a function alongside the Data PlotItem defined above:
    g2.plot(Func("x**2", title="calculated by gnuplot"), d)

    # Save what we just plotted as a color postscript file:
    print("\n******** Generating postscript file 'gnuplot_test1.ps' ********\n")
    g2.hardcopy('gnuplot_test_plot.ps', color=1)

    # Demonstrate a 3-d plot:
    g3 = Gnuplot(debug=1)
    # set up x and y values at which the function will be tabulated:
    x = arange(35)/2.0
    y = arange(30)/10.0 - 1.5
    # Make a 2-d array containing a function of x and y.  First create
    # xm and ym which contain the x and y values in a matrix form that
    # can be `broadcast' into a matrix of the appropriate shape:
    xm = x[:,newaxis]
    ym = y[newaxis,:]
    m = (sin(xm) + 0.1*xm) - ym**2
    g3('set style data lines')
    g3('set hidden')
    g3('set contour base')
    g3.xlabel('x')
    g3.ylabel('y')
    g3.splot(GridData(m,x,y))

    # Delay so the user can see the plots:
    sys.stderr.write("Three plots should have appeared on your screen "
                     "(they may be overlapping).\n"
                     "Please press return to continue...\n")
    sys.stdin.readline()

    # ensure processes and temporary files are cleaned up:
    del g1, g2, g3, d

    # Enable the following code to test the old-style gnuplot interface
    if 0:
        # List of (x, y) pairs
        plot([(0.,1),(1.,5),(2.,3),(3.,4)])

        # List of y values, file output
        print("\n            Generating postscript file 'gnuplot_test2.ps'\n")
        plot([1, 5, 3, 4], file='gnuplot_test2.ps')

        # Two plots; each given by a 2d array
        x = arange(10, typecode=Float)
        y1 = x**2
        y2 = (10-x)**2
        plot(transpose(array([x, y1])), transpose(array([x, y2])))

