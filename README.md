# MakeMyGate
Program dedicated to analysing data stored as 2D histograms written by  
Tomasz Marchlewski (marchlewski.tomasz@gmail.com)  
and  
Roman Szenborn (roman.szenborn@gmail.com).

Some screenshots of MMG's  
First screenshot. Upper spectrum is projection of matrix, bottom shows gated spectrum created by 
"red"(peak) and "blue"(background) gates. Bottom spectrum shows first few states in rotational band.
![124Xe rotational band](http://image.prntscr.com/image/b2c36f6ea231475d958a9bb2138c105d.png "124Xe rotational band")

Same as above, but for 60Co decay lines.
![60Co decay lines](http://image.prntscr.com/image/426f7a633af64ce88bcb1711ee99fb3f.png "60Co decay lines")

Short gif showing possibilities of MMG
![60Co decay lines as gif](http://i.imgur.com/VZnE7QX.gif "gif")



## --Quick guide--
1. General info
2. MakeMyGate requirements
3. Installing python and modules
    1. Anaconda distribution
    2. Miniconda 
    3. Package manager
4. Known problems

### 1\. General info  
MakeMyGate is a program dedicated to visualise and
slice 2D coincidence matrices. It strongly relies on
pyqtgraph library (info here: http://www.pyqtgraph.org/).
MMG is designed to easily slice coincidence matrix with
background subtraction.

### 2\.MakeMyGate requires:  
- python2.7
- python modules: numpy, scipy, pyqtgraph*
    (they are not included by default,
    altought numpy and scipy are included in Anaconda)

### 3\. Installing python and modules  
3.i. Easiest way to obtain python2.7 with numpy and scipy
is to download and install Anaconda distribution,
which is available here https://www.continuum.io/downloads
Be sure to download version with python2.7.
Anaconda distribution has numpy and scipy by default,
still you need to install pyqtgraph by typing in teminal:
"[sudo] conda install pyqtgraph=0.9.10 scipy=0.15.1"
Conda can be used to install other required modules.

### 3.ii. If you don't want to download and install whole Anaconda 
package (for example if you don't have much disk space) you 
can go for Miniconda distribution, which is 
available here: https://conda.io/miniconda.html
After installing it type in terminal:
"[sudo] conda install numpy scipy=0.15.1 pyqtgraph=0.9.10"

### 3.iii. Another option is to use apt, yum or another 
package manager. First make sure that you have python 2.7
and pip. If you don't, install it with availible 
package manager. NumPy and SciPy also can be installed 
via manager. To install pyqtgraph type in terminal:
"[sudo] pip install pyqtgraph=0.9.10"
Pip can be used to install other required 
python modules.
   
\* - pyqtgraph needs at least one of these modules: 
PyQt4, PyQt5 or PySide. This requirement should be met
automatically when you install pyqtgraph. It happens that 
the user himself must take care to install at least 
one of these libraries using pip, conda or 
other package manager. 

### 4\. Known prolems  
We recommend to install olders version of pyqtgraph - 0.9.10 instead of 0.10.0.  
Also for better peakfind consider using scipy 0.15.1

Special thanks to Wouter for introducing us to pyqtgraph
