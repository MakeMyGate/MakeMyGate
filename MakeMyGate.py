# -*- coding: utf-8 -*-
"""
@authors: 
tomasz marchlewski 
    (marchlewski.tomasz@gmail.com, marchlewski@slcj.uw.edu.pl) 
roman szenborn 
    (romanszenborn@gmail.com)

@version: beta 2.0

--Quick guide--
1. General info
2. MakeMyGate requirements
3. Installing python and modules
    a. Anaconda distribution
    b. Miniconda 
    c. Package manager

1. MakeMyGate is a program dedicated to visualise and
slice 2D coincidence matrices. It strongly relies on
pyqtgraph library (info here: http://www.pyqtgraph.org/).
MMG is designed to easily slice coincidence matrix with
background subtraction.

2. MakeMyGate requires:
- python2.7
- python modules: numpy, scipy, pyqtgraph*
    (they are not included by default,
    altought numpy and scipy are included in Anaconda)

3.a. Easiest way to obtain python2.7 with numpy and scipy
is to download and install Anaconda distribution,
which is available here https://www.continuum.io/downloads
Be sure to download version with python2.7.
Anaconda distribution has numpy and scipy by default,
still you need to install pyqtgraph by typing in teminal:
"[sudo] conda install pyqtgraph"
Conda can be used to install other required modules.

3.b. If you don't want to download and install whole Anaconda 
package (for example if you don't have much disk space) you 
can go for Miniconda distribution, which is 
available here: https://conda.io/miniconda.html
After installing it type in terminal:
"[sudo] conda install numpy scipy pyqtgraph"

3.c. Another option is to use apt, yum or another 
package manager. First make sure that you have python 2.7
and pip. If you don't, install it with availible 
package manager. NumPy and SciPy also can be installed 
via manager. To install pyqtgraph type in terminal:
"[sudo] pip install pyqtgraph"
Pip can be used to install other required 
python modules.
   
* - pyqtgraph needs at least one of these modules: 
PyQt4, PyQt5 or PySide. This requirement should be met
automatically when you install pyqtgraph. It happens that 
the user himself must take care to install at least 
one of these libraries using pip, conda or 
other package manager. 


Special thanks to Wouter for introducing us to pyqtgraph
"""



from __future__ import division
import pyqtgraph as pg
from pyqtgraph.widgets.GraphicsLayoutWidget import GraphicsLayoutWidget
from pyqtgraph.Qt import QtCore, QtGui
from pyqtgraph.parametertree import Parameter, ParameterTree
import numpy as np 
import sys
import struct as sct
import platform
from scipy.signal import find_peaks_cwt
from scipy.optimize import leastsq


### roi information ##
class roi(object):
    def __init__(self, roiCenter, roiWidth, color, fill, roiType):
        self.roiCenter = roiCenter
        self.roiWidth = roiWidth
        self.color = color
        self.fill = fill
        self.roiType = roiType
        self.addRoiToPlot()
        self.addLabelToPlot()
       
    def addRoiToPlot(self):
        self.roiRegion = pg.LinearRegionItem([
            self.roiCenter-self.fill*self.roiWidth//2,
            self.roiCenter+self.fill*self.roiWidth//2],
            brush=self.color)
        if self.roiType == 'group': #group roi must be below roi+ and roi-
            self.roiRegion.setZValue(-5)
        window.vbUpper.addItem(self.roiRegion)
        self.roiRegion.sigRegionChanged.connect(self.removeThisRoiOnShake)
        
    def askWidth(self):
        a = int(self.roiRegion.getRegion()[0])
        b = int(self.roiRegion.getRegion()[1])
        return np.absolute(b - a) + 1
        
    def addLabelToPlot(self):
        self.roiLabel = pg.TextItem(
            text = str(self.askWidth()), color=(200, 200, 200), angle=0)
        self.roiLabel.setZValue(20)
        top = window.vbUpper.getViewBox().viewRange()[1][1]
        position = int(self.roiRegion.getRegion()[0])
        self.roiLabel.setPos(position, top)
        window.vbUpper.addItem(self.roiLabel)
        
    def updateRoiLabel(self):
        self.roiLabel.setText(str(self.askWidth()))
        top = window.vbUpper.getViewBox().viewRange()[1][1]
        position = int(self.roiRegion.getRegion()[0])
        self.roiLabel.setPos(position, top)

    def removeThisRoi(self):
        window.vbUpper.removeItem(self.roiRegion)
        window.vbUpper.removeItem(self.roiLabel)
        
    def removeThisRoiOnShake(self):
        if window.isShakeRemoveActive == True:
            window.vbUpper.removeItem(self.roiRegion)
            window.vbUpper.removeItem(self.roiLabel)
            if self.roiType == 'plus':
                window.plusRoiList.remove(self)
            if self.roiType == 'minus':
                window.minusRoiList.remove(self)
            if self.roiType == 'group':
                window.groupRoiList.remove(self)
        
    def sliceMatrix(self):
        a = int(self.roiRegion.getRegion()[0])
        b = int(self.roiRegion.getRegion()[1] + 1)
        return np.sum(window.matrix[:,a:b], axis = 1)
        
    def isInRegion(self, region):
        if (self.roiRegion.getRegion()[0] > region[0] and \
            self.roiRegion.getRegion()[1] < region[1]):
            return True
        else:
            return False

### Additional Spectrum Object ###
class SubWindow(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self, parent)
        self.ifLowerDisplay = False
        self.ifUpperDisplay = False
        self.ifSave = False
        self.scalingFactorUpper = 1.0
        self.scalingFactorLower = 1.0
        self.nameToDisplay = 'name'
        self.newSpectrumUpper = pg.PlotCurveItem(
            np.arange(0,4097), np.zeros(4096), stepMode = True)
        self.newSpectrumLower = pg.PlotCurveItem(
            np.arange(0,4097), np.zeros(4096), stepMode = True)
        self.createWindow()
        self.p.sigTreeStateChanged.connect(self.change)

    def closeEvent(self, event):
        if self.ifSave:
            event.accept()
        else:
            self.cancelButtonFunct()
        
    def clearSpectra(self):
        window.vbUpper.removeItem(self.newSpectrumUpper)
        window.vbLower.removeItem(self.newSpectrumLower)
        window.vbUpper.legend.removeItem(self.newSpectrumUpper.name())        
        window.vbLower.legend.removeItem(self.newSpectrumLower.name())
        self.removeFromMenu()
        
    def removeFromMenu(self):
        window.spectrumMenu.removeAction(self.spectrumButton)

    def okButtonFunct(self):
        buttonText = 'remove ' + str(self.nameToDisplay)
        self.spectrumButton = QtGui.QAction(
            'test button', window, triggered=self.clearSpectra)
        self.spectrumButton.setText(buttonText)
        window.spectrumMenu.addAction(self.spectrumButton)
        self.ifSave = True
        self.newSpectrumUpper.setData(
            np.arange(0,len(self.loadedSpe) + 1), 
            self.loadedSpe*self.scalingFactorUpper,
            name = self.nameToDisplay)
        self.newSpectrumLower.setData(
            np.arange(0,len(self.loadedSpe) + 1), 
            self.loadedSpe*self.scalingFactorLower,
            name = self.nameToDisplay)
        self.close()
    
    def cancelButtonFunct(self):
        window.vbUpper.removeItem(self.newSpectrumUpper)
        window.vbLower.removeItem(self.newSpectrumLower)
        window.displaySpectra.pop()
        try:
            window.vbUpper.legend.removeItem(self.lastUpperSpe)
        except:
            0
        try:
            window.vbLower.legend.removeItem(self.lastLowerSpe)
        except:
            0
        self.close()

    def createWindow(self):
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)

        ### param tree definition
        params = [{
            'name': 'Display options', 'type': 'group',
            'children': [
                {'name': 'Load SPE', 'type': 'action'},
                {'name': 'Display name', 'type': 'str', 
                 'value': self.nameToDisplay, 
                 'tip': "set spectrum's name (optional)"},
                {'name': 'Spectrum scaling Upper', 'type': 'str', 
                 'value': self.scalingFactorUpper, 
                 'tip': "multiply upper spectrum. 'auto' is an option."},
                {'name': 'Spectrum scaling Lower', 'type': 'str', 
                 'value': self.scalingFactorLower, 
                 'tip': "multiply lower spectrum. 'auto' is an option."},
                {'name': 'Color', 'type': 'color', 'value': "FFFFFF", 
                 'tip': "Pick spectrum's color"},
                {'name': 'Upper', 'type': 'bool', 'value': False, 
                 'tip': "Display spectrum on matrix projection"},
                {'name': 'Lower', 'type': 'bool', 'value': False, 
                 'tip': "Display spectrum on gated spectrum"}]}]
        self.p = Parameter.create(name='params', type='group', children=params)      
        self.t = ParameterTree()
        self.t.setParameters(self.p, showTop=False)
        self.layout.addWidget(self.t)
            
        self.resize(400,340)

        okButton = QtGui.QPushButton('OK', clicked = self.okButtonFunct)
        cancelButton = QtGui.QPushButton(
            'Cancel', clicked = self.cancelButtonFunct)
        self.layout.addWidget(okButton)
        self.layout.addWidget(cancelButton)  
        
    def change(self, param, changes):
        for param, change, data in changes:
            path = self.p.childPath(param)
            
            ### Load SPE button clicked
            if path[1] == 'Load SPE':
                fileName = QtGui.QFileDialog.getOpenFileName(
                    self, "Open file","",
                    "Radware SPE (*.spe);;Text file(*.txt);;Any file(*)")
                if (fileName[-4:] == '.spe' or fileName[-4:] == '.err'):
                    with open(fileName,'rb') as f:
                        loadedSpeRaw = f.read()
                    self.loadedSpe = np.array(
                        sct.unpack('<'+4096*'f',loadedSpeRaw[36:-4]))
                elif(fileName[-4:] == '.txt'):
                    with open(fileName, 'r') as f:                            
                        self.loadedSpe = f.read().split()
                else:
                    print 'unknown format'
                self.newSpectrumUpper.setData(
                    np.arange(0,len(self.loadedSpe) + 1), 
                    self.loadedSpe, name = self.nameToDisplay)
                self.newSpectrumLower.setData(
                    np.arange(0,len(self.loadedSpe) + 1), 
                    self.loadedSpe, name = self.nameToDisplay)
                self.nameToDisplay = fileName
                     
            ### Display name
            elif path[1] == 'Display name':
                self.nameToDisplay = str(data)
                self.newSpectrumLower.setData(
                    np.arange(0,len(self.loadedSpe) + 1), 
                    self.loadedSpe*self.scalingFactorLower,
                    name = self.nameToDisplay)
                self.newSpectrumUpper.setData(
                    np.arange(0,len(self.loadedSpe) + 1), 
                    self.loadedSpe*self.scalingFactorUpper,
                    name = self.nameToDisplay)    
            ### Scaling factors
            elif path[1] == 'Spectrum scaling Upper':
                if data == 'auto':
                    originalSpeMax = float(np.max(window.matrixProjectionX))
                    newSpeMax = float(np.max(self.loadedSpe))
                    self.scalingFactorUpper = originalSpeMax/newSpeMax
                    print 'Auto scaling factor for projection spectrum: ' 
                    + str(self.scalingFactorUpper)
                else:
                    try:
                        self.scalingFactorUpper = float(data)
                    except:
                        print 'Input error: must be number(float) or "auto"'
                        return
                self.newSpectrumUpper.setData(
                    np.arange(0,len(self.loadedSpe) + 1), 
                    self.loadedSpe*self.scalingFactorUpper,
                    name = self.nameToDisplay)

            elif path[1] == 'Spectrum scaling Lower':
                if data == 'auto':
                    originalSpeMax = float(np.max(window.dataToPlot))
                    newSpeMax = float(np.max(self.loadedSpe))
                    self.scalingFactorLower = originalSpeMax/newSpeMax
                    print 'Auto scaling factor for gated spectrum: ' 
                    + str(self.scalingFactorLower)
                else:
                    try:
                        self.scalingFactorLower = float(data)
                    except:
                        print 'Input error: must be number(float) or "auto"'
                        return
                self.newSpectrumLower.setData(
                    np.arange(0,len(self.loadedSpe) + 1), 
                    self.loadedSpe*self.scalingFactorLower,
                    name = self.nameToDisplay)
                
            ### Color
            elif path[1] == 'Color':
                self.spectrumColor = data.getRgb()[:3]
                self.newSpectrumUpper.setPen(self.spectrumColor)                       
                self.newSpectrumLower.setPen(self.spectrumColor)
            
            ### Upper chekbox
            elif path[1] == 'Upper':
                self.ifUpperDisplay = data
                if data:
                    window.vbUpper.addItem(self.newSpectrumUpper)
                    self.lastUpperSpe =  self.newSpectrumUpper.name()
                else:
                    window.vbUpper.removeItem(self.newSpectrumUpper)
                    window.vbUpper.legend.removeItem(self.lastUpperSpe)
    
            ### Lower checkbox
            elif path[1] == 'Lower':
                self.ifLowerDisplay = data
                if data:
                    window.vbLower.addItem(self.newSpectrumLower)
                    self.lastLowerSpe =  self.newSpectrumLower.name()
                else:
                    window.vbLower.removeItem(self.newSpectrumLower)
                    window.vbLower.legend.removeItem(self.lastLowerSpe)
                    #window.vbLower.legend.items.pop()

## PeakFind parameters window 
class pfParamsWindow(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self, parent)
        self.ifSave = False
        self.fwhmLow = window.minPeakWidth
        self.fwhmHigh = window.maxPeakWidth
        self.noiseLevel = window.noisePeakWidth
        self.createWindow()
        self.p.sigTreeStateChanged.connect(self.change)
        
    def closeEvent(self, event):
        print 'pf params changed'
        window.peakFindFunct()

    def okButtonFunct(self):
        window.minPeakWidth = self.fwhmLow
        window.maxPeakWidth = self.fwhmHigh
        window.noisePeakWidth = self.noiseLevel
        self.close()
    
    def cancelButtonFunct(self):
        self.close()

    def createWindow(self):
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)

        ### param tree definition
        params = [{
            'name': 'Peak find parameters', 'type': 'group',
            'children': [
                {'name': 'min expected FWHM (channels)', 'type': 'float', 
                 'value': self.fwhmLow, 'step' : 0.5, 
                 'tip':"lowest possible FWHM given in channels"},
                {'name': 'max expected FWHM (channels)', 'type': 'float', 
                 'value': self.fwhmHigh, 'step' : 0.5, 
                 'tip': "highest possible FWHM given in channels"},
                {'name': 'Noise level', 'type': 'float', 
                 'value': self.noiseLevel, 'step' : 0.05, 
                 'tip': "Cuts off noises. Higher value = more peaks"}]}]
        self.p = Parameter.create(name='params', type='group', children=params)      
        self.t = ParameterTree()
        self.t.setParameters(self.p, showTop=False)
        self.layout.addWidget(self.t)
            
        self.resize(400,340)

        okButton = QtGui.QPushButton('OK', clicked = self.okButtonFunct)
        cancelButton = QtGui.QPushButton(
            'Cancel', clicked = self.cancelButtonFunct)
        self.layout.addWidget(okButton)
        self.layout.addWidget(cancelButton)  
        
    def change(self, param, changes):
        for param, change, data in changes:
            path = self.p.childPath(param)
                                
            ### Display name
            if path[1] == 'FWHM low (channels)':
                self.fwhmLow = data
                
            ### Upper chekbox
            elif path[1] == 'FWHM high (channels)':
                self.fwhmHigh = data
    
            ### Lower checkbox
            elif path[1] == 'Noise level':
                self.noiseLevel = data

### loading custom matrix
class loadCustomMatrix(QtGui.QWidget): #under development
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self, parent)
        self.matSizeX = 4096
        self.matSizeY = 4096
        self.dataOrder = 'C'
        self.dataType = 'H'
        self.dataEndian = '<'
        self.skipFirstBytes = 0
        self.skipLastBytes = 0
        self.createWindow()
        self.p.sigTreeStateChanged.connect(self.change)
        
    def closeEvent(self, event):
        print 'Load custom matix: closing'

    def okButtonFunct(self):
        self.close()
    
    def cancelButtonFunct(self):
        self.close()    
        
    def createWindow(self):
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)

        ### param tree definition
        params = [{
            'name': 'Load custom matrix', 'type': 'group',
            'children': [
                {'name': 'Load matrix', 'type': 'action'},
                {'name': 'Dimension X', 'type': 'int', 
                 'value': self.matSizeX, 'step' : 1, 
                 'tip': "matrix X length"},
                {'name': 'Dimension Y', 'type': 'int', 
                 'value': self.matSizeY, 'step' : 1, 
                 'tip': "matrix Y length"},
                {'name': 'Order', 'type': 'str', 
                 'value': self.dataOrder, 
                 'tip': "Data order: C, read rows or Fortran, read columns"},
                {'name': 'Data type', 'type': 'str', 
                 'value': self.dataType, 
                 'tip': "Type of single matrix cell (ex. H - \
                 2byte unsigned intigers, I - 4byte unsigned intigers)"},
                {'name': 'Endian type', 'type': 'str', 
                 'value': self.dataEndian, 'tip': "Data endian type"},
                {'name': 'Skip ... first bytes', 
                 'type': 'int', 'value': self.skipFirstBytes, 'step' : 1, 
                 'tip': "Don't read first ... bytes of loaded file"},
                {'name': 'Skip ... last bytes', 'type': 'int', 
                 'value': self.skipLastBytes, 'step' : 1, 
                 'tip': "Don't read last ... bytes of loaded file"}]}]
        self.p = Parameter.create(name='params', type='group', children=params)      
        self.t = ParameterTree()
        self.t.setParameters(self.p, showTop=False)
        self.layout.addWidget(self.t)
            
        self.resize(500,640)
        
        textbrowser = QtGui.QTextBrowser()
        textbrowser.append(
            '    Visit https://docs.python.org/2/library/struct.html for more \
            detailed informations about data types and endians.\n\n\
            Example: in case of .mat 4096x4096 matrix data type is H, \
            because matrix cells are 2 byte unsigned intigers. \
            Endian should be set to "<" (little endian).')
        textbrowser.setOpenExternalLinks(True)
        self.layout.addWidget(textbrowser)
        
        okButton = QtGui.QPushButton('OK', clicked = self.okButtonFunct)
        cancelButton = QtGui.QPushButton(
            'Cancel', clicked = self.cancelButtonFunct)
        self.layout.addWidget(okButton)
        self.layout.addWidget(cancelButton)  

    def readBinaryMatrix(self, f):
        f.seek(0)
        dataFormat = self.dataEndian + self.dataType
        self.matrix = np.fromfile(f, dtype=dataFormat)
        self.matrix = self.matrix[:]
        if self.skipFirstBytes:
            self.matrix = self.matrix[self.skipFirstBytes:]
        if self.skipLastBytes:
            self.matrix = self.matrix[:-self.skipLastBytes]          
        if self.dataOrder == 'F':
            self.matrix = self.matrix.reshape(
                (self.matSizeX,self.matSizeY), order="F")
        else:
            self.matrix.shape = (self.matSizeX,self.matSizeY)
        window.matrix = self.matrix
        window.showMatrix(1)
        
    def readNonBinaryMatrix(self, f):
        print 'under construction'
                        
    def change(self, param, changes):
        for param, change, data in changes:
            path = self.p.childPath(param)
            
            ### Load Matrix button clicked
            if path[1] == 'Load matrix':
                fileName = QtGui.QFileDialog.getOpenFileName(
                    self, "Open file","", "Any file(*)")
                ### checking if file is binary or text type
                with open(fileName,'rb') as f:
                    isBinary = False
                    for block in f:
                        if '\0' in block:
                            isBinary = True
                            break      
                    if isBinary:
                        self.readBinaryMatrix(f)
                    else:
                        self.readNonBinaryMatrix(f)
                print 'is ' + str(fileName) + ' binary: ' + str(isBinary)
                     
            ### Dimension X
            elif path[1] == 'Dimension X':
                self.matSizeX = int(data)
    
            ### Dimension Y
            elif path[1] == 'Dimension Y':
                self.matSizeY = int(data)

            ### Data order
            elif path[1] == 'Order':
                self.dataOrder = str(data)
                
            ### Data type
            elif path[1] == 'Data type':
                self.dataType = str(data)
            
            ### Endian type
            elif path[1] == 'Endian type':
                self.endianType = str(data)
    
            ### skip first bytes of file
            elif path[1] == 'Skip ... first bytes':
                self.skipFirstBytes = int(data)

            ### skip last bytes of file
            elif path[1] == 'Skip ... last bytes':
                self.skipLastBytes = int(data)
                                                      
### Main window and functions ###
class MainWindow(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self, parent)
        self.refreshTime = 200 #refresh interval in ms
        self.energyCalibAxis = 0.5 #default energy calibration: 0.5keV/channel
        self.plusRoiList  = [] #list of plus rois
        self.minusRoiList = [] #list of minus rois
        self.groupRoiList = [] #list of groups
        self.displaySpectra = [] #additional spectra to display
        self.programRunning = True #start/stop function needs this
        self.isShakeRemoveActive = False 
        self.legendVisible = False
        self.setupUserInterface() #creates GUI
        self.minPeakWidth = 5 #for peak find
        self.maxPeakWidth = 25 #for peak find
        self.noisePeakWidth = 0.1 #for peak find
        self.peakFindActive = False #for peak find auto refresh
        self.peaksLabelsUpper = [] #for peak find
        self.peaksLabelsLower = [] #for peak find
        self.ifTranspose = False #start with untransposed matrix
        self.additionalFunctionsMenu() #functions not usable for most users        
        
    def setupUserInterface(self):
        """ Initialise the User Interface """  
        # upper frame
        self.upperFrame = QtGui.QFrame()
        self.upperFrameLayout = QtGui.QHBoxLayout()
        self.upperFrame.setLayout(self.upperFrameLayout)
        self.upperFrame.setLineWidth(0)
        self.upperFrame.setFrameStyle(QtGui.QFrame.Panel)
        self.upperFrameLayout.setContentsMargins(0,0,0,0)
  
        # upper frame contents    
        self.viewUpper = GraphicsLayoutWidget()
        self.upperFrameLayout.addWidget(self.viewUpper)
        self.vbUpper = pg.PlotItem(title='Matrix projection')
        self.energyAxisUpper = self.vbUpper.axes["top"]["item"]#additional axis
        self.energyAxisUpper.setScale(self.energyCalibAxis)
        self.energyAxisUpper.show()
        self.viewUpper.addItem(self.vbUpper)
               
        # lower frame
        self.lowerFrame = QtGui.QFrame()
        self.lowerFrameLayout = QtGui.QHBoxLayout()
        self.lowerFrame.setLayout(self.lowerFrameLayout)
        self.lowerFrame.setLineWidth(0)
        self.lowerFrame.setFrameStyle(QtGui.QFrame.Panel)
        self.lowerFrameLayout.setContentsMargins(0,0,0,0)
        
        # lower frame content        
        self.viewLower = GraphicsLayoutWidget()
        self.lowerFrameLayout.addWidget(self.viewLower)
        self.vbLower = pg.PlotItem(title='Gated spectrum')
        self.energyAxisLower = self.vbLower.axes["top"]["item"]#additional axis
        self.energyAxisLower.setScale(self.energyCalibAxis)
        self.energyAxisLower.show()
        self.viewLower.addItem(self.vbLower)
        
        # UI window (containing left and right frames)
        UIwindow = QtGui.QWidget(self)
        UIwindowLayout = QtGui.QHBoxLayout()
        UIwindowSplitter = QtGui.QSplitter(QtCore.Qt.Vertical)
        UIwindowLayout.addWidget(UIwindowSplitter)
        UIwindow.setLayout(UIwindowLayout)
        self.setCentralWidget(UIwindow)
        UIwindowSplitter.addWidget(self.upperFrame)
        UIwindowSplitter.addWidget(self.lowerFrame)
        
        #Legend items
        self.vbUpper.addLegend()
        self.vbUpper.legend.anchor((1,0), (1,0))
        self.vbUpper.legend.hide()
        
        self.vbLower.addLegend()  
        self.vbLower.legend.anchor((1,0), (1,0))
        self.vbLower.legend.hide()

        # Status bar
        self.windowStatusBar = QtGui.QStatusBar()
        self.currentNameStatus = QtGui.QLabel(' ') #matrix name and path
        self.transposeStatus = QtGui.QLabel(' ') #is transposed?
        self.moveToRemoveStatus = QtGui.QLabel(' ') #is move to remove ON?
        self.windowStatusBar.insertPermanentWidget(
            0, self.currentNameStatus, 2)
        self.windowStatusBar.insertPermanentWidget(
            1, self.transposeStatus, 2)
        self.windowStatusBar.insertPermanentWidget(
            2, self.moveToRemoveStatus, 2)               
        self.setStatusBar(self.windowStatusBar)
          
        # Application window
        self.windowTitle = 'MakeMyGate v2.0'
        self.setWindowTitle(self.windowTitle)
        self.resize(1300,600)
   
        # Window menus      
        self.createMenus()    
        self.createActions()
    
    def createMenus(self):
        # Menus list
        menubar = self.menuBar()
        self.fileMenu = menubar.addMenu('File')
        self.roiMenu = menubar.addMenu('ROI')
        self.optionsMenu = menubar.addMenu('Options')
        self.spectrumMenu = menubar.addMenu('Spectrums')
        self.utilitiesMenu = menubar.addMenu('Fitting')
        self.aboutMenu = menubar.addMenu('About')

    def createActions(self):   
        # file Menu
        self.exitAct = QtGui.QAction(
            "Quit", self, shortcut="Ctrl+Q",
            statusTip="Exit the application")
        self.loadMatrix = QtGui.QAction("Load matrix", self, shortcut="Ctrl+L")
        self.saveSpe = QtGui.QAction("Save SPE", self, shortcut="Ctrl+S")
        self.saveRoiListToFile = QtGui.QAction("Save ROIs to file", self)
        self.loadRoiList = QtGui.QAction("Load ROIs from file", self)
        self.loadCustomMatrix = QtGui.QAction("Load custom matrix", self)
        fileMenuActions = [
            self.loadMatrix, self.exitAct, self.saveSpe, 
            self.saveRoiListToFile, self.loadRoiList,
            self.loadCustomMatrix]
        fileMenuActFuncs = [
            self.loadMatrixFunct, self.close, self.saveSpeFunct,
            self.saveRoiListToFileFunct, self.loadRoiListFunct,
            self.loadCustomMatrixFunct]
        for i in xrange(len(fileMenuActions)):
            action = fileMenuActions[i]
            function = fileMenuActFuncs[i]
            action.triggered[()].connect(function)
        
        self.fileMenu.addAction(self.loadMatrix)
        self.fileMenu.addAction(self.loadCustomMatrix)
        self.fileMenu.addAction(self.saveSpe)
        self.fileMenu.addAction(self.saveRoiListToFile)
        self.fileMenu.addAction(self.loadRoiList)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.exitAct)

        # ROI menu
        self.addRoiPlus = QtGui.QAction("Add ROI+", self, shortcut="Ctrl++")
        self.addRoiMinus = QtGui.QAction("Add ROI-", self, shortcut="Ctrl+-")
        self.removeRoiPlus = QtGui.QAction(
            "Remove last ROI+", self, shortcut="Ctrl+Shift++")
        self.removeRoiMinus = QtGui.QAction(
            "Remove last ROI-", self, shortcut="Ctrl+Shift+-")
        self.removeAllPlusRois = QtGui.QAction("Remove every ROI+", self)
        self.removeAllMinusRois = QtGui.QAction("Remove every ROI-", self)
        self.removeAllRois = QtGui.QAction("Remove every ROI", self)
        self.addGroupRoi = QtGui.QAction("Add group", self, shortcut="Ctrl+*")
        self.removeGroupRoi = QtGui.QAction(
            "Remove last group", self, shortcut="Ctrl+Shift+*")
        self.removeAllGroupRoi = QtGui.QAction("Remove every group", self)
        self.shakeRoiRemove = QtGui.QAction("Move ROI to remove it: OFF", self)
        roiMenuActions = [
            self.addRoiPlus, self.addRoiMinus, self.removeRoiPlus, 
            self.removeRoiMinus, self.removeAllPlusRois, 
            self.removeAllMinusRois, self.removeAllRois, self.addGroupRoi, 
            self.removeGroupRoi, self.removeAllGroupRoi, self.shakeRoiRemove]
        roiMenuActFuncs = [
            self.addRoiPlusFunct, self.addRoiMinusFunct, 
            self.removeRoiPlusFunct, 
            self.removeRoiMinusFunct, self.removeAllPlusRoisFunct, 
            self.removeAllMinusRoisFunct, self.removeAllRoisFunct, 
            self.addGroupRoiFunct, self.removeGroupRoiFunct, 
            self.removeAllGroupRoiFunct, self.shakeRoiRemoveFunct]
        for i in xrange(len(roiMenuActions)):
            action = roiMenuActions[i]
            function = roiMenuActFuncs[i]
            action.triggered[()].connect(function)    
        
        self.roiMenu.addAction(self.addRoiPlus)
        self.roiMenu.addAction(self.addRoiMinus)
        self.roiMenu.addAction(self.removeRoiPlus)
        self.roiMenu.addAction(self.removeRoiMinus)
        self.roiMenu.addAction(self.addGroupRoi)
        self.roiMenu.addAction(self.removeGroupRoi)
        self.roiMenu.addSeparator()
        self.roiMenu.addAction(self.removeAllPlusRois)
        self.roiMenu.addAction(self.removeAllMinusRois)
        self.roiMenu.addAction(self.removeAllGroupRoi)
        self.roiMenu.addAction(self.removeAllRois)   
        self.roiMenu.addSeparator()
        self.roiMenu.addAction(self.shakeRoiRemove)

        ## Options menu
        self.setRefreshInterval = QtGui.QAction(
            "Set refresh interval", self)
        self.startStopRefresh = QtGui.QAction(
            "Start/Stop refreshing", self, shortcut="Ctrl+X")
        self.setCalibration = QtGui.QAction("Set energy calibration", self)
        self.peakFind = QtGui.QAction("Peak find", self, shortcut="Ctrl+P")
        self.peakFindParams = QtGui.QAction(
            "Adjust peak find parameters", self)
        self.transposeMatrix = QtGui.QAction(
            "Transpose Matrix", self, shortcut="Ctrl+T")
        self.displayLegend = QtGui.QAction("Display legend", self)
        optionsMenuActions = [
            self.setRefreshInterval, self.startStopRefresh,
            self.setCalibration, self.peakFind, self.peakFindParams,
            self.transposeMatrix, self.displayLegend]
        optionsMenuFuncs = [
            self.setRefreshIntervalFunct, self.startStopRefreshFunct,
            self.setCalibrationFunct, self.peakFindFunct, 
            self.peakFindParamsFunct,
            self.transposeMatrixFunct, self.displayLegendFunct]
        for i in xrange(len(optionsMenuActions)):
            action = optionsMenuActions[i]
            function = optionsMenuFuncs[i]
            action.triggered[()].connect(function)  
        
        self.optionsMenu.addAction(self.setRefreshInterval)
        self.optionsMenu.addAction(self.startStopRefresh)      
        self.optionsMenu.addAction(self.displayLegend)
        self.optionsMenu.addAction(self.setCalibration)
        self.optionsMenu.addAction(self.transposeMatrix)
        self.optionsMenu.addSeparator()
        self.optionsMenu.addAction(self.peakFind)
        self.optionsMenu.addAction(self.peakFindParams)

        # Additional spectrums menu
        self.addSpectrum = QtGui.QAction("Display additional spectrum", self)
        self.removeSpectra = QtGui.QAction("Remove all spectra", self)
        spectrumMenuActions = [self.addSpectrum, self.removeSpectra]
        spectrumMenuFuncs = [self.addSpectrumFunct, self.removeSpectrumFunct]
        for i in xrange(len(spectrumMenuActions)):
            action = spectrumMenuActions[i]
            function = spectrumMenuFuncs[i]
            action.triggered[()].connect(function)
        
        self.spectrumMenu.addAction(self.addSpectrum)
        self.spectrumMenu.addSeparator()
        self.spectrumMenu.addAction(self.removeSpectra)
        self.spectrumMenu.addSeparator()
        
        # Utilities menu        
        self.utiRoiAdd = QtGui.QAction("Add Fit ROI", self, shortcut="F2")
        self.utiRoiRemove = QtGui.QAction("Remove Fit ROI", self)
        self.fitPeak = QtGui.QAction("Fit peak in ROI", self, shortcut="F4")
        self.fitNextPeak = QtGui.QAction("Fit 2nd peak in ROI", self)
        self.areaUnderPeak = QtGui.QAction("Calculate area", self)
        self.bgRoi = QtGui.QAction("Add background ROI", self, shortcut="F3")
        self.bgRoiRemove = QtGui.QAction("Remove background ROI", self)
        self.pasternakShap = QtGui.QAction("Paternak Shape", self)
        self.pasternakSingls = QtGui.QAction("Paternak Singlsh", self)
        utilitiesMenuActions = [
            self.utiRoiAdd, self.utiRoiRemove, 
            self.fitPeak, self.fitNextPeak,
            self.areaUnderPeak, self.bgRoi,
            self.bgRoiRemove, self.pasternakShap,
            self.pasternakSingls]
        utilitiesMenuFuncs = [
            self.utiRoiFunct, self.utiRoiRemoveFunct, 
            self.fitPeakFunct, self.fitNextPeakFunct,
            self.areaUnderPeakFunct, self.bgRoiFunct,
            self.bgRoiRemoveFunct, self.pasternakShape,
            self.pasternakSinglsh]
        for i in xrange(len(utilitiesMenuActions)):
            action = utilitiesMenuActions[i]
            function = utilitiesMenuFuncs[i]
            action.triggered[()].connect(function)
            
        self.utilitiesMenu.addAction(self.utiRoiAdd)         
        self.utilitiesMenu.addAction(self.utiRoiRemove)
        self.utilitiesMenu.addAction(self.bgRoi)
        self.utilitiesMenu.addAction(self.bgRoiRemove)
        self.utilitiesMenu.addSeparator()
        self.utilitiesMenu.addAction(self.fitPeak)
        self.utilitiesMenu.addAction(self.fitNextPeak)
        self.utilitiesMenu.addAction(self.areaUnderPeak)
#        self.utilitiesMenu.addAction(self.pasternakShap)
#        self.utilitiesMenu.addAction(self.pasternakSingls)
                
        # About page
        self.aboutAct = QtGui.QAction("&About", self, shortcut='F1',
                                      triggered=self.onAbout)
        self.aboutMenu.addAction(self.aboutAct)

    def additionalFunctionsMenu(self):
        print 'experimental functions visible in menu'
        self.utilitiesMenu.addAction(self.pasternakShap)
        self.utilitiesMenu.addAction(self.pasternakSingls)
        ## reading file with custom matrix data
        # MakeMyGate_mattype.inp
        basicMatType = [['2byte uint matrix','mat',4096,4096,'C','H','<',0,0],
                        ['4byte uint matrix','m4b',4096,4096,'C','I','<',0,0]]
        try:
            self.mattypeFile = basicMatType
            with open ('MakeMyGate_mattype.inp','r') as f:
                mattypeFile = np.loadtxt(f, dtype=str, delimiter='\n')
            mattypeFile = np.array(mattypeFile).reshape(-1,9)
            self.mattypeFile.extend(mattypeFile.tolist())
        except IOError:
            print (
            'MakeMyGate_mattype.inp not found.\n'\
            + 'Only basic type of matrices available.')
        except ValueError:
            print (
            'MakeMyGate_mattype.inp structure not matching required pattern.\n'
            +' Only basic types of matrices available')
        except:
            print (
            'unknown error occurred while trying\n'\
            +'to load MakeMyGate_mattype.inp')
            
    def onAbout(self):
        """ About message""" 
        QtGui.QMessageBox.about(
            self, 'Info page',
            """
            <b>MakeMyGate alpha2.0</b>
            <p>Tomasz Marchlewski (marchlewski@slcj.uw.edu.pl, 
            marchlewski.tomasz@gmail.com)</p>
            <p>Roman Szenborn (roman.szenborn@gmail.com)</p>
            <p>May 2015</p>
            """)

    def loadMatrixFunct(self):
        fileTypes = ''
        for fileType in self.mattypeFile:
            fileTypes += str(fileType[0]) + ' *.' + str(fileType[1]) \
                +'(*.' + str(fileType[1]) + ');;'
        fileTypes = fileTypes[:-2] #removes last(harmful)";;" from string
        fileName, filter = QtGui.QFileDialog.getOpenFileNameAndFilter(
            self, "Load matrix", "", fileTypes)
        print 'Loading matrix:'
        print fileName, filter        
        
        #Status bar and transpose flag update
        self.currentNameStatus.setText(str(fileName))
        self.transposeStatus.setText(' ')
        self.ifTranspose = False

        #matching filter with known matrix formats        
        for possibleFilter in self.mattypeFile:
            if (str(filter).startswith(possibleFilter[0])):
                with open(fileName, 'rb') as f:
                    dataType = possibleFilter[6]+possibleFilter[5]
                    self.matrix = np.fromfile(f,dtype=dataType)
                self.matrix.shape = (possibleFilter[2], possibleFilter[3])
                self.showMatrix(1)
                        
    def loadCustomMatrixFunct(self):
        self.customMatLoad = loadCustomMatrix()
        self.customMatLoad.show()
                
    ## show matrix projections after loading
    def showMatrix(self, *args):
        if args: #load new matrix
            self.matrixProjectionX = np.sum(self.matrix, axis = 0)
            self.matrixProjectionY = np.sum(self.matrix, axis = 1)
            self.removeAllRoisFunct()
            self.vbUpper.clear()
            self.vbLower.clear()
            self.upperSpe = pg.PlotCurveItem(
                np.arange(0, len(self.matrixProjectionX)+1), 
                self.matrixProjectionX,stepMode=True,
                name = 'mat proj')
            try:
                self.vbUpper.legend.removeItem('mat proj')
            except:
                0
            self.lowerSpe = pg.PlotCurveItem(
                np.arange(0, len(self.matrixProjectionY)+1), 
                self.matrixProjectionY,stepMode=True,
                name = 'gated spe')
            try:
                self.vbLower.legend.removeItem('gated spe')
            except:
                0
            self.vbUpper.addItem(self.upperSpe)
            self.vbLower.addItem(self.lowerSpe)
            self.dataToPlot = self.matrixProjectionY
        else: #just refresh the view after transpose
            self.matrixProjectionX = np.sum(self.matrix, axis = 0)
            self.vbUpper.removeItem(self.upperSpe)
            self.upperSpe = pg.PlotCurveItem(
                np.arange(0, len(self.matrixProjectionX)+1), 
                self.matrixProjectionX,stepMode=True)
            self.vbUpper.addItem(self.upperSpe)          
            self.matrixProjectionY = np.sum(self.matrix, axis = 1)
            self.vbLower.removeItem(self.lowerSpe)
            self.lowerSpe = pg.PlotCurveItem(
                np.arange(0, len(self.matrixProjectionY)+1), 
                self.matrixProjectionY,stepMode=True)
            self.vbLower.addItem(self.lowerSpe)
            self.dataToPlot = self.matrixProjectionY
                                            
    def saveSpeFunct(self): # saves gated spe, error spe and rois list to file
        fileTypes = ("Radware SPE (*.spe);;Text file (*.txt)")
        FileName, filter = QtGui.QFileDialog.getSaveFileNameAndFilter(
            self, 'Save file', '', fileTypes)
        print FileName, filter        
        if FileName == '':
            print 'no file'
            return
        if filter == 'Radware SPE (*.spe)':
            if platform.system() == 'Windows':
                FileName1 = FileName[:-4] + str('.spe')
                FileName2 = FileName[:-4] + str('.err')
                self.saveRoiListToFileFunct(FileName[:-4])
            elif platform.system() == 'Linux':
                FileName1 = FileName + str('.spe')
                FileName2 = FileName + str('.err') 
                self.saveRoiListToFileFunct(FileName)
            else:
                print 'saving spe (Linux naming system)'
                FileName1 = FileName + str('.spe')
                FileName2 = FileName + str('.err') 
                self.saveRoiListToFileFunct(FileName)
            inSpeName = str(8*' ') + str(FileName1)
            speHeader = sct.pack(
                'I8s6I',24,inSpeName[-8:],4096,1,1,1,24,16384)
            if len(self.groupRoiList) == 0:                
                speToPack = self.caclGatedSpe()
                errToPack = self.calcErrSpe()
            else:
                speToPack = self.calcGatedSpeGroups()
                errToPack = self.calcErrSpeGroups()
            packedSpe = sct.pack('<'+4096*'f',*speToPack)
            speEnding = sct.pack('I',16384)
            with open(FileName1, 'wb') as f:
                f.write(speHeader + packedSpe + speEnding)
            packedErr = sct.pack('<'+4096*'f',*errToPack)
            with open(FileName2, 'wb') as f:
                f.write(speHeader + packedErr + speEnding)
        if filter == 'Text file (*.txt)':
            FileName1 = FileName + str('.txt')
            toSaveText1 = self.caclGatedSpe()
            np.savetxt(FileName1, toSaveText1, fmt='%s', delimiter='/u')    
            self.saveRoiListToFileFunct(FileName)
            toSaveText2 = self.calcErrSpe()
            FileName2 = FileName + str('err.txt')
            np.savetxt(FileName2, toSaveText2, fmt='%s', delimiter='/u')

    def addRoiPlusFunct(self):
        viewBoxRange = np.array(self.vbUpper.getViewBox().viewRange()[0])
        roiCenter = (viewBoxRange[1] + viewBoxRange[0])//2
        roiWidth = viewBoxRange[1] - viewBoxRange[0]
        color = (255,0,0,90)
        newRoi = roi(roiCenter, roiWidth, color, 0.03, 'plus')
        self.plusRoiList.append(newRoi)
        
    def addRoiMinusFunct(self):
        viewBoxRange = np.array(self.vbUpper.getViewBox().viewRange()[0])
        roiCenter = (viewBoxRange[1] + viewBoxRange[0])//2
        roiWidth = viewBoxRange[1] - viewBoxRange[0]
        color = (0,0,255,90)
        newRoi = roi(roiCenter, roiWidth, color, 0.03, 'minus')
        self.minusRoiList.append(newRoi)
        
    def removeRoiPlusFunct(self):
        self.plusRoiList.pop().removeThisRoi()
        
    def removeRoiMinusFunct(self):
        self.minusRoiList.pop().removeThisRoi()     

    def addGroupRoiFunct(self):
        viewBoxRange = np.array(self.vbUpper.getViewBox().viewRange()[0])
        roiCenter = (viewBoxRange[1] + viewBoxRange[0])//2
        roiWidth = viewBoxRange[1] - viewBoxRange[0]
        color = (0,255,0,90)
        newRoi = roi(roiCenter, roiWidth, color, 0.03, 'group')
        self.groupRoiList.append(newRoi)
        print (self.groupRoiList[-1].askWidth())
    
    def removeGroupRoiFunct(self):
        self.groupRoiList.pop().removeThisRoi()
        
    def removeAllPlusRoisFunct(self):
        for i in xrange(len(self.plusRoiList)):
            self.plusRoiList.pop().removeThisRoi()
            
    def removeAllMinusRoisFunct(self):
        for i in xrange(len(self.minusRoiList)):
            self.minusRoiList.pop().removeThisRoi()

    def removeAllGroupRoiFunct(self):
        for i in xrange(len(self.groupRoiList)):
            self.groupRoiList.pop().removeThisRoi()

    def removeAllRoisFunct(self):
        self.removeAllMinusRoisFunct()
        self.removeAllPlusRoisFunct()
        self.removeAllGroupRoiFunct()

    def shakeRoiRemoveFunct(self):
        if self.isShakeRemoveActive == False:
            print 'remove ROI on move: ON'
            self.isShakeRemoveActive = True
            self.moveToRemoveStatus.setText('Move ROI to remove it: ON')
            self.shakeRoiRemove.setText("Move ROI to remove it: ON")
        else:
            print 'remove ROI on move: OFF'
            self.isShakeRemoveActive = False
            self.moveToRemoveStatus.setText(' ')            
            self.shakeRoiRemove.setText("Move ROI to remove it: OFF")
   
    def setRefreshIntervalFunct(self):
        DialogWindow = QtGui.QInputDialog(self)
        Text, ok = DialogWindow.getText(
            self, "set refresh interval in ms",
            "Time in ms", 
            QtGui.QLineEdit.Normal,
            "0 for instant refreshing")
        if ok and len(Text):
            self.resfreshTime = int(Text)
            timer.setInterval(self.resfreshTime)
        else:
            print 'canceled or input error'
        
    def startStopRefreshFunct(self):
        if self.programRunning: 
            self.programRunning = False
            print 'Stop'
            iniText = self.windowTitle
            self.setWindowTitle(
                iniText + " | not working - press Ctrl+X to start")
            timer.stop()
        else:
            self.lowerPlotUpdate()
            self.programRunning = True
            print 'Start'
            iniText = self.windowTitle
            self.setWindowTitle(iniText.replace(
                " | not working - press Ctrl+X to start",''))
            timer.start()

    def setCalibrationFunct(self):
        print 'set calibration'
        DialogWindow = QtGui.QInputDialog(self)
        Text, ok = DialogWindow.getText(
            self, "Energy calibration",
            "keV/channel", 
            QtGui.QLineEdit.Normal,
            "0.5 default")
        if ok and len(Text):
            self.energyCalibAxis = float(Text)
            self.energyAxisUpper.setScale(self.energyCalibAxis)
            self.energyAxisLower.setScale(self.energyCalibAxis)
        else:
            print 'canceled or input error'

    def peakFindFunct(self):
        print 'pf start'
        self.peakFindUpper()
        try:
            self.peakFindLower()
        except AttributeError:
            print "no gates found - gated spectrum doesn't exist"
        
    def peakFindUpper(self):#executed on pf activation and on pf params change
        print 'upper PF'
        for label in self.peaksLabelsUpper:
            self.vbUpper.removeItem(label) 
        self.peaksListUpper = find_peaks_cwt(
            self.matrixProjectionX, 
            np.arange(self.minPeakWidth,self.maxPeakWidth), 
            noise_perc=self.noisePeakWidth)           
        for peak0 in self.peaksListUpper:
            peak = peak0
            self.roiLabel = pg.TextItem(
                text = str(peak*self.energyCalibAxis),
                color=(200, 200, 200), angle=0)
            self.roiLabel.setZValue(20)
            top = self.matrixProjectionX[peak]
            position = peak
            self.roiLabel.setPos(position, top)
            self.vbUpper.addItem(self.roiLabel)
            self.peaksLabelsUpper.append(self.roiLabel)

    #executed together with self.peakFindUpper()
    #and on lower plot refresh if self.peakFindActive = True
    #auto peak find disabled, freezes everyting
    def peakFindLower(self):
        print 'lower PF'
        #clear labels
        for label in self.peaksLabelsLower:
            self.vbLower.removeItem(label)
            
        #create peak list
        self.peaksListLower = find_peaks_cwt(
            self.dataToPlot, 
            np.arange(self.minPeakWidth,self.maxPeakWidth), 
            noise_perc=self.noisePeakWidth)
        for peak0 in self.peaksListLower:
            peak = peak0
            self.roiLabel = pg.TextItem(
                text = str(peak*self.energyCalibAxis),
                color=(200, 200, 200), angle=0)
            self.roiLabel.setZValue(20)
            top = self.dataToPlot[peak]
            position = peak
            self.roiLabel.setPos(position, top)
            self.vbLower.addItem(self.roiLabel)
            self.peaksLabelsLower.append(self.roiLabel)

    def peakFindParamsFunct(self):
        print 'pf params change'
        self.pfWindow = pfParamsWindow()
        self.pfWindow.show()

    def transposeMatrixFunct(self):
        print 'transpose matrix'
        self.matrix = self.matrix.transpose()
        self.showMatrix()
        if self.ifTranspose:
            self.ifTranspose = False
            self.transposeStatus.setText(' ')
        else:
            self.ifTranspose = True
            self.transposeStatus.setText('TRANSPOSED')

    def displayLegendFunct(self):
        if self.legendVisible:
            self.vbUpper.legend.hide()
            self.vbLower.legend.hide()
            self.legendVisible = False
        else:
            self.vbUpper.legend.show()
            self.vbLower.legend.show()
            self.legendVisible = True

    def calcErrSpe(self): #calculates error spectrum ^2
        backgroundSpe = np.zeros(4096)
        if(len(self.plusRoiList)):
            rawSpe = self.plusRoiList[0].sliceMatrix()
            suppresionUp = self.plusRoiList[0].askWidth()
            self.plusRoiList[0].updateRoiLabel()
            for roiPlus in self.plusRoiList[1:]:
                rawSpe += roiPlus.sliceMatrix()
                suppresionUp += roiPlus.askWidth()
            dataToPlot = rawSpe
        if(len(self.minusRoiList)): 
            suppresionDown = self.minusRoiList[0].askWidth()
            backgroundSpe = self.minusRoiList[0].sliceMatrix()
            self.minusRoiList[0].updateRoiLabel()
            for roiMinus in self.minusRoiList[1:]:
                backgroundSpe += roiMinus.sliceMatrix()
                suppresionDown += roiMinus.askWidth()
            suppresionFactor = float(suppresionUp)/float(suppresionDown)
            dataToPlot = rawSpe + suppresionFactor**2*backgroundSpe
        return dataToPlot        

    def caclGatedSpe(self): #calculates gated spectrum
        backgroundSpe = np.zeros(4096)
        if(len(self.plusRoiList)):
            rawSpe = self.plusRoiList[0].sliceMatrix()
            suppresionUp = self.plusRoiList[0].askWidth()
            self.plusRoiList[0].updateRoiLabel()
            for roiPlus in self.plusRoiList[1:]:
                rawSpe += roiPlus.sliceMatrix()
                suppresionUp += roiPlus.askWidth()
            dataToPlot = rawSpe
        if(len(self.minusRoiList)): 
            suppresionDown = self.minusRoiList[0].askWidth()
            backgroundSpe = self.minusRoiList[0].sliceMatrix()
            self.minusRoiList[0].updateRoiLabel()
            for roiMinus in self.minusRoiList[1:]:
                backgroundSpe += roiMinus.sliceMatrix()
                suppresionDown += roiMinus.askWidth()
            suppresionFactor = float(suppresionUp)/float(suppresionDown)
            dataToPlot = rawSpe - suppresionFactor*backgroundSpe
        return dataToPlot

    def calcGatedSpeGroups(self): #calculates gated spe with groups
        finalSpe = np.zeros(4096)
        for group in self.groupRoiList:
            plusRegions = [roi for roi in self.plusRoiList \
                if roi.isInRegion(group.roiRegion.getRegion())]
            minusRegions = [roi for roi in self.minusRoiList \
                if roi.isInRegion(group.roiRegion.getRegion())]
            plusSpe = np.zeros(4096)
            minusSpe = np.zeros(4096)
            upFactor = 0.
            downFactor = 0.            
            for roi in plusRegions:
                plusSpe += roi.sliceMatrix()
                upFactor += roi.askWidth()
            for roi in minusRegions:
                minusSpe += roi.sliceMatrix()
                downFactor += roi.askWidth()
            supFact = float(upFactor)/float(downFactor)
            finalSpe = finalSpe + plusSpe - supFact*minusSpe
        return finalSpe

    def calcErrSpeGroups(self): #calculates error spe^2 with groups
        finalSpe = np.zeros(4096)
        for group in self.groupRoiList:
            plusRegions = [roi for roi in self.plusRoiList \
                if roi.isInRegion(group.roiRegion.getRegion())]
            minusRegions = [roi for roi in self.minusRoiList \
                if roi.isInRegion(group.roiRegion.getRegion())]
            plusSpe = np.zeros(4096)
            minusSpe = np.zeros(4096)
            upFactor = 0.
            downFactor = 0.            
            for roi in plusRegions:
                plusSpe += roi.sliceMatrix()
                upFactor += roi.askWidth()
            for roi in minusRegions:
                minusSpe += roi.sliceMatrix()
                downFactor += roi.askWidth()
            supFact = float(upFactor)/float(downFactor)
            finalSpe = finalSpe + plusSpe + supFact**2*minusSpe
        return finalSpe

    #saves all rois to text file with ".rl" ext
    def saveRoiListToFileFunct(self, *args):
        if args:
            FileName0 = str(args[0])
        else:
            FileName0 = QtGui.QFileDialog.getSaveFileName(
                self, "Save ROI list", "", "Text file (*.rl)")
        FileName = FileName0 + str('.rl')
        ToSaveText = str('PlusRois: ') + str(len(self.plusRoiList)) + '\n' + \
        str('MinusRois: ') + str(len(self.minusRoiList)) + '\n' + \
        str('GroupRois: ') + str(len(self.groupRoiList)) + '\n'
        for roi in self.plusRoiList:
            ToSaveText += str(int(roi.roiRegion.getRegion()[0])) + ' ' + \
            str(int(roi.roiRegion.getRegion()[1])) + '\n'
        for roi in self.minusRoiList:
            ToSaveText += str(int(roi.roiRegion.getRegion()[0])) + ' ' + \
            str(int(roi.roiRegion.getRegion()[1])) + '\n'
        for roi in self.groupRoiList:
            ToSaveText += str(int(roi.roiRegion.getRegion()[0])) + ' ' + \
            str(int(roi.roiRegion.getRegion()[1])) + '\n'
        with open(FileName, 'w') as f:
            f.write(ToSaveText)

    def loadRoiListFunct(self): #lazy, but works
        fileName = QtGui.QFileDialog.getOpenFileName(
            self, "Open file","", "Roi list(*.rl);;any(*)")
        with open(fileName, 'r') as f:
            roiFromFile = f.read().splitlines()
        pCount = int((roiFromFile[0]).split()[-1]) 
        mCount = int((roiFromFile[1]).split()[-1]) 
        gCount = int((roiFromFile[2]).split()[-1]) 
        for roiLine in roiFromFile[3:3+pCount]:
            viewBoxRange = roiLine.split()
            roiCenter = (int(viewBoxRange[1]) + int(viewBoxRange[0]))//2
            roiWidth = int(viewBoxRange[1]) - int(viewBoxRange[0])
            color = (255,0,0,90)
            newRoi = roi(roiCenter, roiWidth, color, 1, 'plus')
            self.plusRoiList.append(newRoi)
        for roiLine in roiFromFile[3+pCount:3+pCount+mCount]:
            viewBoxRange = roiLine.split()
            roiCenter = (int(viewBoxRange[1]) + int(viewBoxRange[0]))//2
            roiWidth = int(viewBoxRange[1]) - int(viewBoxRange[0])
            color = (0,0,255,90)
            newRoi = roi(roiCenter, roiWidth, color, 1, 'minus')
            self.minusRoiList.append(newRoi)
        for roiLine in roiFromFile[3+pCount+mCount:3+pCount+mCount+gCount]:
            viewBoxRange = roiLine.split()
            roiCenter = (int(viewBoxRange[1]) + int(viewBoxRange[0]))//2
            roiWidth = int(viewBoxRange[1]) - int(viewBoxRange[0])
            color = (0,255,0,90)
            newRoi = roi(roiCenter, roiWidth, color, 1, 'group')
            self.groupRoiList.append(newRoi)

    #### Display additional spectrum ###
    def addSpectrumFunct(self):
        self.newWindow = SubWindow()
        self.displaySpectra.append(self.newWindow)
        self.newWindow.show()
        
    def removeSpectrumFunct(self):
        print 'remove spectrum'
        for spe in self.displaySpectra:
            spe.clearSpectra()

    def bgRoiFunct(self):
        try:
            self.vbLower.removeItem(self.bgRoi)
        except:
            0
        viewBoxRange = np.array(self.vbLower.getViewBox().viewRange()[0])
        roiCenter = (viewBoxRange[1] + viewBoxRange[0])//2
        roiWidth = viewBoxRange[1] - viewBoxRange[0]
        color = (255,255,255,40)
        self.bgRoi = pg.LinearRegionItem(
            [roiCenter-0.03*roiWidth//2, roiCenter+0.03*roiWidth//2],
            brush=color)
        self.bgRoi.setZValue(-10)
        self.vbLower.addItem(self.bgRoi)            

    def bgRoiRemoveFunct(self):
        try:
            self.vbLower.removeItem(self.bgRoi)
            del self.bgRoi
        except:
            0           

    def utiRoiFunct(self):
        try:
            for i in xrange(len(self.peaksMarks)):
                self.vbLower.removeItem(self.peaksMarks.pop())
        except:
            self.peaksMarks = []
        try:
            self.vbLower.removeItem(self.utiRoi)
        except:
            0
        viewBoxRange = np.array(self.vbLower.getViewBox().viewRange()[0])
        roiCenter = (viewBoxRange[1] + viewBoxRange[0])//2
        roiWidth = viewBoxRange[1] - viewBoxRange[0]
        color = (0,255,255,70)
        self.utiRoi = pg.LinearRegionItem(
            [roiCenter-0.03*roiWidth//2, roiCenter+0.03*roiWidth//2],
            brush=color)
        self.utiRoi.setZValue(-5)
        self.vbLower.addItem(self.utiRoi)

    def utiRoiRemoveFunct(self):
        try:
            self.vbLower.removeItem(self.utiRoi)
            self.vbLower.removeItem(self.fitBackground)
            self.vbLower.removeItem(self.fitGauss)
            self.vbLower.removeItem(self.fitLabel)
            self.vbLower.removeItem(self.newFitGauss)
            self.vbLower.removeItem(self.fitLabel2)
        except:
            print 'utiRoiRemovefail'

    # fits single gaussian peak        
    def fitPeakFunct(self):
        self.roiLimits = self.utiRoi.getRegion()
        try:
            self.bgLimits = self.bgRoi.getRegion()
        except:
            self.bgLimits = self.roiLimits
        self.bgLen = int(self.bgLimits[1]) - int(self.bgLimits[0]) + 1
        self.speRegion = self.dataToPlot[
            int(self.roiLimits[0]):int(self.roiLimits[1] + 1)]
        #level of background
        self.bgPoints = self.dataToPlot[
            int(self.bgLimits[0])],self.dataToPlot[int(self.bgLimits[1])]
        fitBgParams = np.polyfit([0,self.bgLen],self.bgPoints,1)
        self.bgSpe = np.arange(self.bgLen)
        self.bgSpe = self.bgSpe*fitBgParams[0] + fitBgParams[1]
        self.fitSpeAxis = np.arange(
            int(self.roiLimits[0]),int(self.roiLimits[1] + 1)+1,1)    
        self.fitBgAxis = np.arange(
            int(self.bgLimits[0]),int(self.bgLimits[1] + 1)+1,1)
        
        try:
            self.fitBackground.setData(self.fitBgAxis, self.bgSpe)
        except:
            self.fitBackground = pg.PlotCurveItem(
                self.fitBgAxis, self.bgSpe, stepMode = True)
            self.fitBackground.setPen(0,0,255)
            self.vbLower.addItem(self.fitBackground)

        def gaussToFit(x,a,mu,sig):
            return a*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))\
                +self.bgSpeCut[x]
            
        def guassToFitErr(p,x,y):
            return y - gaussToFit(x,*p)

        if int(self.roiLimits[0]) - int(self.bgLimits[0]) >= 0:
            bgLeft = int(self.roiLimits[0]) - int(self.bgLimits[0])
        else:
            print 'Left Peak ROI out of Background ROI'
            return
        if int(self.bgLimits[1]) - int(self.roiLimits[1]) >= 0:
            bgRight = int(self.bgLimits[1]) - int(self.roiLimits[1])
        else:
            print 'Right Peak ROI out of Background ROI'
            return
        
        if bgRight and bgLeft:
            self.bgSpeCut = self.bgSpe[bgLeft:-bgRight]
        elif bgLeft:
            self.bgSpeCut = self.bgSpe[bgLeft:]
        elif bgRight:
            self.bgSpeCut = self.bgSpe[:-bgRight]
        else:
            self.bgSpeCut = self.bgSpe[:]
            
        guess = np.max(self.speRegion), self.roiLimits[1]-self.roiLimits[0], 5.
        x = xrange(len(self.speRegion))   
        y = self.speRegion
                    
        out = leastsq(guassToFitErr, guess, args = (x,y))
        self.firstSigma = out[0][2] #for 2nd peak fitting

        self.gaussToPlot = gaussToFit(xrange(len(self.speRegion)),*out[0])
        try:
            self.fitGauss.setData(self.fitSpeAxis, self.gaussToPlot)
        except:
            self.fitGauss = pg.PlotCurveItem(
                self.fitSpeAxis, self.gaussToPlot, stepMode = True)
            self.fitGauss.setPen(255,0,0)
            self.fitGauss.setZValue(10)
            self.vbLower.addItem(self.fitGauss)      
           
        area = np.sum(self.gaussToPlot) - np.sum(self.bgSpeCut)
        centroid1 = (out[0][1] + int(self.roiLimits[0]))*self.energyCalibAxis
        fwhm = 2.355*abs(out[0][2])*self.energyCalibAxis
        peaktext = str('Area=%d \nE=%.1fkeV \nFWHM=%.2fkeV' 
            % (area, centroid1, fwhm))
        print '------\n 1st peak fit:'
        print peaktext
        try:
            self.fitLabel.setText(peaktext)
            top = self.dataToPlot[int(out[0][1])+ int(self.roiLimits[0])]
            position = centroid1/self.energyCalibAxis
            self.fitLabel.setPos(position, top)
        except:
            self.fitLabel = pg.TextItem(
                text = peaktext,
                color=(255, 255, 255), angle=0,
                anchor=(0, 1))
            top = self.dataToPlot[int(out[0][1])+ int(self.roiLimits[0])]
            position = centroid1/self.energyCalibAxis
            self.fitLabel.setZValue(10)
            self.vbLower.addItem(self.fitLabel)
            self.fitLabel.setPos(position, top)

    # fits another gaussian peak if possible
    def fitNextPeakFunct(self):
        newSpeToFit = self.speRegion - self.gaussToPlot + self.bgSpeCut
        
        def gaussToFit(x,a,mu,sig):
            return a*np.exp(-np.power(x - mu, 2.) \
                /(2 * np.power(sig, 2.)))+self.bgSpeCut[x]
            
        def guassToFitErr(p,x,y):
            penalize = 0.
            if abs(p[2]) > self.firstSigma + 0.7 \
                or abs(p[2]) < self.firstSigma - 0.7:
                if abs(p[2]) < 1:
                    penalize = 0#abs(p[2])**-5*np.sum(y - self.bgSpe)
                else:
                    penalize = 0#abs(p[2])**5*np.sum(y - self.bgSpe)
            penalize1 = 0.
            if p[0] < 0.:
                penalize1 = 0#abs(p[0])*np.sum(y - self.bgSpe)
            return y - gaussToFit(x,*p) + penalize + penalize1

        guess = [np.max(newSpeToFit), self.roiLimits[1] - self.roiLimits[0],
                 self.firstSigma]            
        x = xrange(len(self.speRegion))   
        y = newSpeToFit               
        out = leastsq(guassToFitErr, guess, args = (x,y))

        self.newGaussToPlot = gaussToFit(xrange(len(self.speRegion)),*out[0])
        try:
            self.newFitGauss.setData(self.fitSpeAxis, self.newGaussToPlot)
        except:
            self.newFitGauss = pg.PlotCurveItem(
                self.fitSpeAxis, self.newGaussToPlot, stepMode = True)
            self.newFitGauss.setPen(0,255,0)
            self.newFitGauss.setZValue(10)
            self.vbLower.addItem(self.newFitGauss)       

        area = np.sum(self.newGaussToPlot) - np.sum(self.bgSpeCut)
        centroid1 = (out[0][1] + int(self.roiLimits[0]))*self.energyCalibAxis
        fwhm = 2.355*abs(out[0][2])*self.energyCalibAxis
        print '------\n 2nd peak fit:'
        '''print 'Area=' + str(area) + ' energy=' + str(centroid1) + \
        ' FWHM=' + str(fwhm)'''
        peaktext = str('Area=%d \nE=%.1fkeV \nFWHM=%.2fkeV' \
            % (area, centroid1, fwhm))
        print peaktext

        try:
            self.fitLabel2.setText(peaktext)
            top = self.dataToPlot[int(out[0][1])+ int(self.roiLimits[0])]
            position = centroid1/self.energyCalibAxis
            self.fitLabel2.setPos(position, top)
        except:
            self.fitLabel2 = pg.TextItem(
                text = peaktext,
                color=(255, 255, 255), angle=0,
                anchor=(0,1))
            top = self.dataToPlot[int(out[0][1])+ int(self.roiLimits[0])]
            position = centroid1/self.energyCalibAxis
            self.fitLabel2.setZValue(10)
            self.vbLower.addItem(self.fitLabel2)
            self.fitLabel2.setPos(position, top)

    def areaUnderPeakFunct(self):
        print 'calc area'
        self.roiLimits = self.utiRoi.getRegion()
        self.speRegion = self.dataToPlot[
            int(self.roiLimits[0]):int(self.roiLimits[1] + 1)]
        #level of background
        self.bgPoints = self.dataToPlot[
            int(self.roiLimits[0])],self.dataToPlot[int(self.roiLimits[1])]
        fitBgParams = np.polyfit([0,len(self.speRegion)],self.bgPoints,1)
        self.bgSpe = np.arange(len(self.speRegion))
        self.bgSpe = self.bgSpe*fitBgParams[0] + fitBgParams[1]
        self.fitSpeAxis = np.arange(
            int(self.roiLimits[0]),int(self.roiLimits[1] + 1)+1,1)      
        try:
            self.fitBackground.setData(self.fitSpeAxis, self.bgSpe)
        except:
            self.fitBackground = pg.PlotCurveItem(
                self.fitSpeAxis, self.bgSpe, stepMode = True)
            self.fitBackground.setPen(0,0,255)
            self.vbLower.addItem(self.fitBackground)
            
        area = np.sum(self.speRegion - self.bgSpe)
        print 'Area = ' + str(area)

    #creates spectrum for A.A. Pasternak software
    def pasternakShape(self): 
        self.roiLimits = self.utiRoi.getRegion()
        #spe fragment with peak
        self.speRegion = self.dataToPlot[
            int(self.roiLimits[0]):int(self.roiLimits[1] + 1)]
        if len(self.groupRoiList) == 0:                
            errSpe = self.calcErrSpe()
        else:
            errSpe = self.calcErrSpeGroups()
        #err spe
        self.error = errSpe[int(self.roiLimits[0]):int(self.roiLimits[1] + 1)]
        n0 = self.roiLimits[0] #first chan number
        nk = len(self.speRegion) #number of chan

        bg = np.min(self.speRegion) #bg level
        spe = self.speRegion - bg
        err = self.error + bg
        err_spe = err - spe

        FileName = QtGui.QFileDialog.getSaveFileName(
            self, "Save spectrum", "", "Text file (*)")
        ToSaveText = str(int(n0)) + ' ' + str(nk) + '\n'        
        a = np.around(err)
        b = np.around(err_spe)
        with open(FileName, 'w') as f:
            f.write(ToSaveText)
            np.savetxt(f, (a,b), fmt='%d')
            
    # A.A. Pasternak spe for Singlsh
    def pasternakSinglsh(self): 
        print 'test2'
        self.roiLimits = self.utiRoi.getRegion()
        #peak spe
        self.speRegion = self.dataToPlot[
            int(self.roiLimits[0]):int(self.roiLimits[1] + 1)]   
        #bg spe, linear
        self.bgPoints = self.speRegion[0], self.speRegion[1]
        fitBgParams = np.polyfit([0,len(self.speRegion)],self.bgPoints,1)
        self.bgSpe = np.arange(len(self.speRegion))
        self.bgSpe = self.bgSpe*fitBgParams[0] + fitBgParams[1]
        FileName = QtGui.QFileDialog.getSaveFileName(
            self, "Save spectrum", "", "Text file (*)")
        n0 = self.roiLimits[0] 
        nk = len(self.speRegion) 
        ToSaveText = str(int(n0)) + ' ' + str(nk) + '\n'        
        a = np.around(self.speRegion)
        b = np.around(self.bgSpe)
        with open(FileName, 'w') as f:            
            f.write(ToSaveText)
            np.savetxt(f, (a,b), fmt='%d')        

    def refreshAllLabels(self):
        for roi in self.plusRoiList:
            roi.updateRoiLabel()
        for roi in self.minusRoiList:
            roi.updateRoiLabel()
        for roi in self.groupRoiList:
            roi.updateRoiLabel()

    def lowerPlotUpdate(self):
        self.refreshAllLabels()
        try:
            if len(self.groupRoiList) == 0:
                self.dataToPlot = self.caclGatedSpe()
            else:
                self.dataToPlot = self.calcGatedSpeGroups()
            self.lowerSpe.setData(
                np.arange(0,len(self.dataToPlot)+1), 
                self.dataToPlot)
        except:
            # does nothing
            0
                                      
        #auto-peak find disabled, too heavy for cpu
        '''if self.peakFindActive: 
            for label in self.peaksLabelsLower:
                self.vbLower.removeItem(label)
            threading.Thread(target=self.peakFindLower()).start()'''
                    
    def closeEvent(self, event):
        reply = QtGui.QMessageBox.question(
            self, 'Close MakeMyGate',
            "Are you sure to quit?", QtGui.QMessageBox.Yes | 
            QtGui.QMessageBox.No, QtGui.QMessageBox.No)

        if reply == QtGui.QMessageBox.Yes:
            event.accept()
            print 'MakeMyGate: "bye, bye"'
        else:
            event.ignore()   

          
### auto refreshing init
def refreshInit():
    timer.setInterval(window.refreshTime)
    timer.timeout.connect(window.lowerPlotUpdate)    
    timer.start()       
  
def run():
    # PySide fix: Check if QApplication already exists. 
    # Create QApplication if it doesn't exist
    app=QtGui.QApplication.instance()       
    if not app:
        app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    refreshInit()
    app.exec_()

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    timer = QtCore.QTimer()
    refreshInit()
    sys.exit(app.exec_())
