'''
Stuff for driving MS office applications from Python using COM

Currently just Excel but Word will come soon.
'''

from win32com.client import Dispatch
from types import *
from string import uppercase

class Excel:
    '''
    Wrapper for MS Excel derived from that in Python Programing on Win32
    '''
    def __init__(self,filename=None):
        '''
        Open a new Excel spreadsheet optionally associated with a file
        '''
        self.xlApp = Dispatch("Excel.Application")
        if filename:
            self.filename = filename
            self.xlBook = self.xlApp.Workbooks.Open(filename)
        else:
            self.xlBook = self.xlApp.Workbooks.Add()
            self.filename = ''
        self.open = 1

    def save(self, newfilename=None):
        '''
        Save the workbook either to the default file, another file,
        or let Excel query the user where to save it.
        '''
        if newfilename:
            self.filename = newfilename
            self.xlBook.SaveAs(newfilename)
        else:
            self.xlBook.Save()

    def close(self):
        self.xlBook.Close(SaveChanges=0)
        del self.xlApp
        self.open = 0

    def getCell(self, row, col, sheet=1):
        '''
        Returns the value in cell (row,col) or None if it is blank.
        '''
        xlSheet = self.xlBook.Worksheets(sheet)
        return xlSheet.Cells(row,col).Value

    def getCellFormula(self, row, col, sheet=1):
        '''
        Returns the formula in cell (row,col) or the value if
        there is no formula.  If there is no value nor formula,
        None is returned.
        '''
        xlSheet = self.xlBook.Worksheets(sheet)
        result = xlSheet.Cells(row,col).Formula
        if result == '':   # A blank field seems to return a blank string
            result = None
        return result

    def setCell(self, value, row, col, sheet=1):
        '''
        Sets the value in cell (row,col).
        '''
        xlSheet = self.xlBook.Worksheets(sheet)
        xlSheet.Cells(row,col).Value = value

    def getRange(self, row1, col1, row2=None, col2=None, sheet=1):
        '''
        Returns the data in the given range as a 2d array (i.e., as
        a tuple of tuples).  If the bottom corner is not specified
        or is incompletely specified, assume a dimension of 1.
        '''
        if not row2:
            row2 = row1
        if not col2:
            col2 = col1
        xlSheet = self.xlBook.Worksheets(sheet)
        cell1 = xlSheet.Cells(row1,col1)
        cell2 = xlSheet.Cells(row2,col2)
        return xlSheet.Range(cell1,cell2).Value

    def matrixDimensions(self, data):
        '''
        Determine the dimemension of the matrix data which can be a
        scalar, vector, or 2-D matrix. Allows for string data, or for
        matrices in which the first row or column are strings labels
        for series ... so look at the last row to determine the length
        of a row (= number of columns). If the data is a vector then
        it is taken as a row-vector in order to be consistent with how
        the default extension happens when assigning a simple list or
        vector into a rectangular range in Excel.
        '''
        last = None
        n = m = 1
        try:
            n = len(data)
            last = data[-1]
        except TypeError:
            n = m = 1  # We have a scalar

        if last:
            if type(last) == StringType:
                m = n   # Row-vector of strings 
                n = 1
            else:
                try:
                    m = len(last)
                except TypeError:
                    m = n # Row-vector of scalars
                    n = 1
        return (n,m)

    def setRange(self, data, row1=1, col1=1, row2=None, col2=None, sheet=1):
        '''
        Set the range of cells to the given data.
        
        If both corners of the range are specified, the corresponding
        piece of data is copied into the range.  If data is too small,
        then it is mystically extended to fill the range.  E.g., you
        can fill a range from a scalar or a vector.  Vectors are treated
        as row-vectors when filling a rectangular region.
        
        Optionally, you specify only the top-left corner of range in
        row1, cell1 and specify row2<=0 - the other coordinate is figured
        out from the dimension of the data.  This can always be overriden by
        specifying the full range coordinates.

        If no coordinates are given, the data is put into the top left
        of the spreadsheet.

        Returns the range that was set.
        '''
        (n,m) = self.matrixDimensions(data)
        if not row2:
            row2 = row1 + n - 1
        if not col2:
            col2 = col1 + m - 1
        xlSheet = self.xlBook.Worksheets(sheet)
        cell1 = xlSheet.Cells(row1,col1)
        cell2 = xlSheet.Cells(row2,col2)
        xlSheet.Range(cell1,cell2).Value = data
        return (row1, col1, row2, col2)

    def getContiguousRange(self, row1, col1, sheet=1):
        '''
        Returns data in the range which forms a continguous
        block with top-left corner in cell (row1,col1).
        
        Starting from the specified cell, scan down/across
        the first column/row and identify the range bordered
        by blank cells.  Blanks within the region will be
        set to None.
        '''
        xlSheet = self.xlBook.Worksheets(sheet)
        row2 = row1
        while xlSheet.Cells(row2+1,col1).Value not in [None,'']:
            row2 = row2 + 1

        col2 = col1
        while xlSheet.Cells(row1,col2+1).Value not in [None,'']:
            col2 = col2 + 1
        return self.getRange(row1, col1, row2, col2, sheet=sheet)

    def selectRange(self, row1, col1, row2=None, col2=None, sheet=1):
        '''
        Select the range of cells on the specified sheet.  It also
        has to select that sheet as the active worksheet.
        '''
        if not row2:
            row2 = row1
        if not col2:
            col2 = col1
        xlSheet = self.xlBook.Worksheets(sheet)
        xlSheet.Select()
        cell1 = xlSheet.Cells(row1,col1)
        cell2 = xlSheet.Cells(row2,col2)
        xlSheet.Range(cell1,cell2).Select()

    def chartRange(self,  row1, col1, row2, col2, sheet=1,
                   **keys):
        '''
        Chart the data in the specified range.  Additional options
        are processed by chartSelectedRange.
        '''
        self.selectRange(row1, col1, row2, col2, sheet=sheet)
        keys['sheet'] = sheet
        apply(self.chartSelectedRange, (), keys)

    def chartSelectedRange(self,  
                           title=None, xlabel=None, ylabel=None,
                           plotby='columns',
                           charttype='xy',
                           sheet=1,
                           xmin=None, xmax=None,
                           ymin=None, ymax=None,
                           xlog=0, ylog=0):
        '''
        The interface to Excel charts.  Just a few of the capabilities
        are exposed here.

        [The first of a set of options is the default]
        
        plotby = 'columns' ... data series run down columns
        .      = 'rows'    ... across rows

        charttype = 'xy'   ... XY scatter plot with lines and points.
        .           First series is X.  Others are y1, y2, etc.
        .         = 'surface' ... Surfce plot of a scalar function of
        .           two variables.  Data should be a grid of the function.
        .         = 'contour' or 'colorcontour' ... Contour plot of a scalar
        .           function of two variables.  Data should be a grid of
        .           values.
        xmin and xmax = min/max values of the x or category axis
        .         It defaults to autoscale by Excel.  This only applies to
        .         XY plots (since the surfce/contor plots do not use
        .         values for the category axes ... they use string labels)
        ymin and ymax = min/max values of the y or value axis
        .         It defaults to auto by Excel.  Applies to all charts.
        xlog = 0 ... use a linear  for x or category axis.
        .    = 1 ... use a log  (values must be positive)
        .         This only applies to XY plots.
        ylog = 0 ... use a linear  for the value or Y axis
        .    = 1 ... use a log .
        .         Applies to all charts
        
        If the first element of each data series is a string, it is
        used to label the series.  If this string is representable as
        a numerical value you must precede it with a single quote to
        force Excel to treat it as a string.  Note that you must use
        strings.  If you use numbers it will be interpreted as data
        and incorporated into the plot.  For the 2-D plots (xy,
        surface, contour) you can border the actual data on left and
        on the top with strings to label axes.
        '''
        charttypes = {'xy':74, 'surface':83, 'colorcontour':85, 'contour':56}
        try:
            charttype = charttypes[charttype]
        except KeyError:
            print('Excel.chartSelectedRange: Unkown charttype', charttype, ' defaulting to XY')
            charttype = charttypes['xy']

        # Make the chart and set how the data will be interpreted
        # Taking a reference to the active chart does not seemt to work???
        self.xlApp.Charts.Add()
        self.xlApp.ActiveChart.ChartType = charttype
        xlRows=1
        xlColumns=2
        if plotby == 'rows':
            self.xlApp.ActiveChart.PlotBy = xlRows
        elif plotby == 'columns':
            self.xlApp.ActiveChart.PlotBy = xlColumns
        else:
            print('Excel.chartSelectedRange: Unknown plotby', charttype, ' defaulting to columns')
            self.xlApp.ActiveChart.PlotBy = xlColumns

        # Set the title and axis labels
        if title:
            self.xlApp.ActiveChart.HasTitle = 1
            self.xlApp.ActiveChart.ChartTitle.Characters.Text = title
        xlCategory=1
        xlValue=2
        #xlSeries=3
        xlPrimary=1
        #xlSecondary=2
        if xlabel:
            self.xlApp.ActiveChart.Axes(xlCategory,xlPrimary).HasTitle = 1
            self.xlApp.ActiveChart.Axes(xlCategory,xlPrimary).AxisTitle.Characters.Text = xlabel
        if ylabel:
            self.xlApp.ActiveChart.Axes(xlValue,xlPrimary).HasTitle = 1
            self.xlApp.ActiveChart.Axes(xlValue,xlPrimary).AxisTitle.Characters.Text = ylabel

        # Set the axis scale and log options
        xlLinear = 0xffffefdc
        xlLogarithmic=0xffffefdb
        if ymin != None:
            self.xlApp.ActiveChart.Axes(xlValue).MinimumScale = ymin
        if ymax != None:
            self.xlApp.ActiveChart.Axes(xlValue).MaximumScale = ymax
        if ylog:
            self.xlApp.ActiveChart.Axes(xlValue).ScaleType = xlLogarithmic
        if charttype == charttypes['xy']:
            if xmin != None:
                self.xlApp.ActiveChart.Axes(xlCategory).MinimumScale = xmin
            if xmax != None:
                self.xlApp.ActiveChart.Axes(xlCategory).MaximumScale = xmax
            if xlog:
                self.xlApp.ActiveChart.Axes(xlCategory).ScaleType = xlLogarithmic
        # A legend is kinda useful
        self.xlApp.ActiveChart.HasLegend = 1

    def chartData(self, data, row1=1, col1=1, sheet=1, **keys):
        '''
        Simplest interface for creating a chart.  Data is a matrix
        of data.  Paste it into a sheet and plot it.  All arguments
        except the data can be defaulted.  Optional arguments are passed
        to the actual charting function.
        '''
        (n,m) = self.matrixDimensions(data)
        row2 = row1 + n - 1
        col2 = col1 + m - 1
        self.setRange(data, row1, col1, row2, col2, sheet=sheet)
        keys['sheet'] = sheet
        apply(self.chartRange, (row1, col1, row2, col2), keys)

    def a1(self, row, col, absrow=0, abscol=0):
        '''
        Return a string that may be used to adress the cell in
        a formula.  The row and/or column adress may be made absolute
        by setting absrow/col to true values.
        
        Internally we are adressing cells in the spreadsheet using
        integers (row,col), which is what Excel calls R1C1 style
        references.  But, unless the user has turned-on R1C1 style
        adressing (unlikely!) this will not work in formulae
        so we must translate to the usual adressing style, called A1,
        which uses letters for the columns and numbers for the rows,
        writing the column index first.

        E.g., A1 = R1C1 = (1,1), and B3 = R3C2 = (3,2).

        Absolute adresses are preceded with a $ symbol.
        '''
        ar = ac = ''
        if absrow: ar = '$'
        if abscol: ac = '$'
        if col < 1 or col > 256:
            raise RangeError('column index must be in [1,256]')
        (c1,c2) = divmod(col-1,26)
        if c1:
            c = uppercase[c1] + uppercase[c2]
        else:
            c = uppercase[c2]
        r = str(row)
        return ac + c + ar + r
        
    def visible(self):
        '''
        Make the spreadsheet visible.
        '''
        self.xlApp.Visible = 1

    def invisible(self):
        '''
        Make the spreadsheet invisible.
        '''
        self.xlApp.Visible = 0

    def isvisible(self):
        '''
        Returns true if the spreadsheet is visible.
        '''
        return self.xlApp.Visible

    def __del__(self):
        '''
        Destructor ... may be uncessary but it cannot hurt.
        '''
        if self.open:
            self.close()

if __name__ == "__main__":
    from math import *
    import time
    # Make a worksheet and test set/getCell
    xls = Excel()
    print(' Setting cell(2,2) to "Hi"')
    xls.setCell("Hi", 2, 2)
    print(xls.getCell(2,2))
    print(' Setting cell(1,2) to "(1,2)"')
    xls.setCell("(1,2)", 1, 2)
    print(' Setting cell(2,1) to "(1,2)"')
    xls.setCell("(2,1)", 2, 1)
    xls.visible()

    # Test setting a range to a scalar and getting contiguous range
    print(' Setting 9,1,12,2 to 0')
    xls.setRange(0,9,1,12,2)
    print(' Getting same contiguous range back ... expecting matrix(4,2)=0')
    value = xls.getContiguousRange(9,1)
    print(value)

    # Test setting/getting a range from/to a matrix
    n = 3
    m = 5
    x = [0]*n
    for i in range(n):
        x[i] = [0]*m
        for j in range(m):
            x[i][j] = i + j
    print(' Setting range (3:,4:) to ')
    print(x)
    xls.setRange(x,3,4)  # Auto determination of the bottom corner
    print(' Got back from same range ',3,3,3+n-1,4+m-1)
    y = xls.getRange(3,4,3+n-1,4+m-1)
    print(y)

    # Add names for the series that will eventually become the chart
    names = []
    for i in range(m):
        names.append("y%d" % i)
    xls.setRange(names,2,4)

    # Test selecting a range
    print(' Selecting range ', 3,3,3+n-1,4+m-1)
    xls.selectRange(3,4,3+n-1,4+m-1)

    # Test general matrix
    xls.setRange([[1,2],[5,6],["hi","bye"]],1,10,3,11)
    
    # Test making an x-y plot (changes the range selection)
    xls.chartRange(2,4,3+n-1,4+m-1,
                   title='THIS IS THE TITLE',
                   xlabel='XXXXX',
                   ylabel='YYYYY')
    
    # Test making an x-y plot just from the data ... use a
    # second sheet and the simple chart interface
    print(' Creating chart of sin(x) and cos(x) using second sheet')
    n = 20
    m = 3
    h = 2*pi/(n-1)
    data = range(n+1)
    data[0] = ['x', 'sin', 'cos']
    for i in range(n):
        x = i*h
        data[i+1] = (x,sin(x),cos(x))
    xls.chartData(data,sheet=2)

    # Try using a formula to add up the absolute values of the data
    # Use absolute values for the rows but not the columns to test
    # reuse of the formula.
    formula = '=sum('+xls.a1(2,2,absrow=1)+':'+xls.a1(21,2,absrow=1)+')'
    print(' The formula is ', formula)
    xls.setCell('Total',23,1,sheet=2)
    xls.setCell(formula,23,2,sheet=2)
    xls.setCell(formula,23,3,sheet=2)
    # Getting the cell contents back will get the value not the formula
    print(' The formula from the sheet is ', xls.getCellFormula(23,2,sheet=2))
    print(' The value of the formula (sum of sin) is ', xls.getCell(23,2,sheet=2))
    print(' The formula from where there is only the value "Total" is', xls.getCellFormula(23,1,sheet=2))
    print(' The formula from where there is nothing ', xls.getCellFormula(23,4,sheet=2))
    print(' The value from where there is nothing ', xls.getCell(23,4,sheet=2))

    # Make a surface plot by creating a 2-D grid bordered on the
    # left and top with strings to indicate the values.  Note the
    # use of a single quote before the value in the labels in
    # order to force Excel to treat them as strings.
    print(' Create surface chart of exp(-0.1*r*r)*cos(1.3*r)')
    n = 10
    h = 2*pi/(n-1)
    data = range(n+1)
    data[0] = range(n+1)
    data[0][0] = ' '
    for i in range(n):
        x = i*h-pi
        data[i+1] = range(n+1)
        data[0][i+1] = data[i+1][0] = ("'%5.2f" % x)
        for j in range(n):
            y = j*h-pi
            r = sqrt(x*x+y*y)
            data[i+1][j+1] = exp(-0.1*r*r)*cos(1.3*r)
    # Specify (row1,col1) to avoid overwriting the previous data
    # Also, specify axis ranges to make the animation smoother.
    xls.chartData(data,1,5,sheet=2,charttype='surface',
                  ymin=-1,ymax=1) 

    # Animate the chart by periodically updating the data range.
    nloop = 60
    for loop in range(1,nloop+1):
        phase = loop*2*pi/(nloop-1)
        for i in range(n):
            x = i*h-pi
            for j in range(n):
                y = j*h-pi
                r = sqrt(x*x+y*y)
                data[i+1][j+1] = exp(-0.1*r*r)*cos(1.3*r+phase)
        time.sleep(0.5)
        xls.setRange(data,1,5,sheet=2)

    # Finally make a chart with all options set
    print(' Creating chart of sin(x) and cos(x) using second sheet')
    n = 81
    data = range(n+1)
    data[0] = ['Age', 'Wisdom']
    for i in range(n):
        data[i+1] = [i+1, 1.0 + 100.0*exp(-((i-40)**2)/400.0)]
    xls.chartData(data,1,1,sheet=3,plotby='columns',charttype='xy',
                  xmin=1,xmax=80,ymin=1,ymax=100,ylog=1,
                  title='Wisdom vs. Age', xlabel='Age/years',
                  ylabel='Wisdom')

