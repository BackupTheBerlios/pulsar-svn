"""
$plot.py

PULSAR Project:

 Copyright (C) 2006-2007 Jean-Paul Amoureux, Christian Fernandez
 JPA - Unite de Catalyse et Chimie du Solide, Lille, France.
 CF  - Laboratoire Catalyse et Spectrochimie, Caen, France.

LICENSE:

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License (GPL)
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 To read the license please visit http://www.gnu.org/copyleft/gpl.html

PURPOSE OF THIS FILE:

 Class for the Pulsar Script Graphics

"""
__author__ = "C. Fernandez <christian.fernandez@ensicaen.fr>"
__contributors__ =""

import wx

import matplotlib
import matplotlib.numerix as numpy
matplotlib.use('WXagg')
from matplotlib.backends.backend_wxagg import FigureCanvasWx as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as Toolbar
from matplotlib.backends.backend_wxagg import FigureManager
from matplotlib.figure import Figure
import string
import pylab 
import sys,os
from Numeric import *
from scipy.fftpack import *

from debug import *

false=No=no=FALSE=False
true=Yes=yes=TRUE=True

DEBUG=False

#============================================================================
class PlotFigure(wx.Window):
#============================================================================
    def __init__(self,parent,size=wx.DefaultSize, filename=None, options=None):

        wx.Window.__init__( self, parent, id=-1, style=wx.NO_BORDER)
        self.fig = Figure((21,21),60)
        self.canvas = FigureCanvas(self, -1, self.fig)
        self.toolbar = Toolbar(self.canvas)
        self.choice=PlotChoice(self)
        self.choice.mode="Standard"
        self.toolbar.Realize()
        self.parent=parent
        # On Windows, default frame size behaviour is incorrect
        # you don't need this under Linux
        tw, th = self.toolbar.GetSizeTuple()
        fw, fh = self.canvas.GetSizeTuple()
        self.toolbar.SetSize(wx.Size(fw, th))
        
        # Create a figure manager to manage things
        self.figmgr = FigureManager(self.canvas, 1, self)
        self.Fit()

        self.filename=filename
        if DEBUG :
            DEBUG_MSG("plot :"+filename)
        self.n=0
        self.ni=0
        self.vl=0
        self.vl1=0
        self.sw=0
        self.sw1=0
        self.flag_sum=False
        
        try:          
            [dirname,shortfilename]=os.path.split(filename)
            [basename,extension]=os.path.splitext(shortfilename)

            self.title=basename
            self.filename=dirname+"\\"+basename+".spe"            

            if not options:
                self.options={'units': 'ppm',
                              'titles': 'no',
                              'yoffset': '0.001',
                              'yaxis': 'yes'}
            self.grid=False
        except:
            pass

    def GetToolBar(self):
        # You will need to override GetToolBar if you are using an 
        # unmanaged toolbar in your frame
         
        return self.toolbar
    
    def read_spectra(self):
        #frequency=(i+1)*sw/n-sw/2.
        #ppm=frequency*1.e6/abs(vl)
        DEBUG_MSG("read_spectra : "+self.filename)
        
        try: #Handling filename error (Correction of a bug for NEw in file menu)
            f=open(self.filename,"rb")
        except:
            return False
        
        try: # Handling opening erros
            #first the header
            line=f.readline()    # PULSAR string
            DEBUG_MSG(rtn(line))
            line=f.readline()    # NP string
            DEBUG_MSG(rtn(line))
            items=line.split('=')
            n=int(items[1])     
            line=f.readline()    # SW string
            DEBUG_MSG(rtn(line))
            items=line.split('=')
            sw=float(items[1])
            line=f.readline()    # VL
            DEBUG_MSG(rtn(line))
            items=line.split('=')
            vl=float(items[1])
            line=f.readline()    # NI or LB string?
            DEBUG_MSG(rtn(line))
            items=line.split('=')
            ni=1 #by default
            if string.strip(items[0])=="NI":     
                ni=int(items[1])
                line=f.readline()    #SW1
                DEBUG_MSG(rtn(line))
                items=line.split('=')
                sw1=float(items[1])
                line=f.readline()    #VL1
                DEBUG_MSG(rtn(line))
                items=line.split('=')
                vl1=float(items[1])
                line=f.readline()    #SUM
                DEBUG_MSG(rtn(line))
                items=line.split('=')
                flag_sum=eval(string.upper(items[1]))
                line=f.readline()    # LB string
                DEBUG_MSG(rtn(line))
                items=line.split('=')     
            #now we have to read LB
            self.lb=[]
            self.lb.append(float(items[1]))
            #and GB
            line=f.readline()
            DEBUG_MSG(rtn(line))
            items=line.split('=')
            self.gb=[]
            self.gb.append(float(items[1]))
            #and CONCENTRATION
            line=f.readline() 
            DEBUG_MSG(rtn(line))
            items=line.split('=')
            self.concentration=[]
            self.concentration.append(float(items[1]))
            if ni>1:
                for j in range(ni-1):
                    line=f.readline() #LB
                    DEBUG_MSG(rtn(line))
                    items=line.split('=')
                    self.lb.append(float(items[1]))
                    line=f.readline() #GB
                    DEBUG_MSG(rtn(line))
                    items=line.split('=')
                    self.gb.append(float(items[1]))
                    line=f.readline() #concentration
                    DEBUG_MSG(rtn(line))
                    items=line.split('=')
                    self.concentration.append(float(items[1])) 

            self.ydata = []
            self.spec = []
            
            # read TYPE 
            line=f.readline()
            DEBUG_MSG(rtn(line))

            # read YDATA or XDATA
            line=f.readline()  
            DEBUG_MSG(rtn(line))

        except: # There was an error
            WRITE_STRING("***ERROR***: Badly formatted header")
            f.close()
            return false

        if string.find(line,"YDATA")>-1: # check if is a YDATA directive
            try:     
                self.label=string.join(string.split(line)[1:],' ')
                for j in range(ni):
                    line=f.readline()
                    self.ydata.append(float(line))
                line=f.readline()
                DEBUG_MSG(rtn(line))
            except:
                WRITE_STRING("***ERROR***: bad YDATA labels")
                f.close()
                return false
        try:
            for j in range(ni):
                DEBUG_MSG("\treading spectrum"+str(j+1))
                xy = []
                DEBUG_MSG("n %13i"%n)
                for i in range(n):
                    line=f.readline()
                    items=line.split(' ')
                    x=float(items[0])
                    y=float(items[1])
                    xy.append(complex(x,y)*self.concentration[j])
                DEBUG_MSG("end reading "+str(j+1))
                self.spec.append(array(xy))

            line=f.readline()       
            DEBUG_MSG(rtn(line))
            if string.find(line,"END")>-1:
                    DEBUG_MSG("reading file success!")
        except:
            WRITE_STRING("***ERROR***: problem when reading X,Y data in the *.spe file!")
            f.close()
            return false
            
        self.frequency, self.ppm = [], []
        for i in range(n):
            self.frequency.append((i+1)*sw/n-sw/2.)
            self.ppm.append(((i+1)*sw/n-sw/2.)*1.e6/abs(vl))

        #save some data as class properties
        self.n=n
        self.ni=ni
        self.vl=vl
        self.sw=sw
        if ni > 1:
            self.vl1=vl1
            self.sw1=sw1
            self.flag_sum=flag_sum
            
        #perform the fft and broadening
        addlb(ni, self.spec, self.sw, self.lb, self.gb, self.concentration)

        #sum of the spectra
        if self.flag_sum:
            p=[]
            for i in range (n):
                sum=complex(0.,0.)
                for j in range(ni):
                    sum=sum+self.spec[j][i]
                p.append(sum)      
            self.spec.append(array(p))
            self.ni=ni+1    
        return True
    
    #------------------------------------------------    
    def plot_data(self,mode=None, offset=None):
    #------------------------------------------------
        """ PLOT_SPECTRA:
        Plot one or more 1D spectra on the same figure
        """

        if not mode :
            mode='reset'
        else:
            mode= string.lower(mode)

        if mode == 'reset' :
            try:
                if not self.read_spectra(): return False
            except:
                return False
                
            self.maxi=0
            for j in range(self.ni):
                mxi=max(map(abs,self.spec[j]))
                self.maxi=max(mxi,self.maxi)
            self.offset=0

        mx=self.ni+1
        self.dlg = wx.ProgressDialog("Spectrum plotting",
                       "...wait...",
                       maximum = mx+1,
                       parent=self,
                        # style = wx.PD_CAN_ABORT
                        #| wx.PD_APP_MODAL
                        #| wx.PD_ELAPSED_TIME
                        #| wx.PD_ESTIMATED_TIME
                        #| wx.PD_REMAINING_TIME
                           )

        keepGoing = True
        count = 0
        self.fig.clear()
        self.plotdata = self.fig.add_subplot(111)
        if offset :
            self.offset=offset
        else:
            self.offset=0
        
        if mode == 'reset' or mode=="standard" :               
            if string.lower(self.options['units'])=="hz":
                for i in range(self.ni):
                    spec=[]
                    shift=(i+1)*self.offset*self.maxi           
                    spec=[shift + x.real for x in self.spec[i]]
                    self.plotdata.plot(self.frequency, array(spec))
                    count += 1
                    keepGoing = self.dlg.Update(count)
                self.plotdata.set_xlim(self.frequency[len(self.frequency)-1],self.frequency[0])
                self.plotdata.set_xlabel('Frequency (Hz)')
            else:
                for i in range(self.ni):
                    shift=(i+1)*self.offset*self.maxi           
                    spec=[shift + x.real for x in self.spec[i]]
                    self.plotdata.plot(self.ppm, array(spec))
                    count += 1
                    keepGoing = self.dlg.Update(count)
                self.plotdata.set_xlim(self.ppm[len(self.ppm)-1],self.ppm[0])
                self.plotdata.set_xlabel('Chemical shift (ppm)')
            if string.lower(self.options['yaxis'])!="no": 
                self.plotdata.set_ylabel('Intensity (a.u.)\n')
            else:    
                self.plotdata.set_yticks([])   
                self.plotdata.set_yticklabels([])         
            self.plotdata.set_title(self.title)
            self.plotdata.grid(self.grid)

        #------------------
        elif mode == "popt":

            if self.ni == 1:
                self.dlg.Destroy()
                return

            self.fig.clear()
            self.plotdat=None
            self.plotdata = self.fig.add_subplot(111)

            xticks=[]
            xtickslabel=[]
            spacing=10*(round(log10(self.ni))-1)
            DEBUG_MSG("%5i %5i"%(self.ni, spacing))
            if spacing<1: spacing=1
            for i in range(self.ni):
                X=[float(self.n-j-1+i*self.n)/self.n for j in range(self.n)]
                if self.ydata:
                    if i%spacing==0:
                        xticks.append(float(self.n/2+i*self.n)/self.n)
                        xtickslabel.append(string.strip(str('%10.3f'%self.ydata[i])))
                Y=[k.real for k in self.spec[i]]
                self.plotdata.plot(X, array(Y))
                count += 1
                keepGoing = self.dlg.Update(count)

            self.plotdata.set_xlabel('Index')   
            if self.ydata:
                self.plotdata.set_xticks(xticks)   
                self.plotdata.set_xticklabels(xtickslabel)
                self.plotdata.set_xlabel(self.label)
                
            self.plotdata.set_ylabel('Intensity (a.u.)\n\n')                   
            self.plotdata.set_title("POPT display, "+self.title)
            self.plotdata.grid(self.grid)

         #------------------
        elif mode == "popt 2d":

            if self.ni == 1:
                self.dlg.Destroy()
                return

            self.fig.clear()
            self.plotdata = self.fig.add_subplot(111)

            count += 1
            keepGoing = self.dlg.Update(count)

            X=[]
            Y=[]
            Z=[]
            X=self.ppm
            Y=range(self.ni)
            yticks=[]
            ytickslabel=[]
            spacing=10*(round(log10(self.ni))-1)
            if spacing<1: spacing=1
            maxi=0.
            mini=0.
            for i in range(self.ni):
                Zt=[k.real for k in self.spec[i]]
                mxi=max(map(abs,Zt))
                mni=min(map(abs,Zt))
                maxi=max(maxi,mxi)
                mini=min(mini,mni)
                Z.append(array(Zt))
                if self.ydata:
                    if i%spacing==0:
                        yticks.append(i)
                        ytickslabel.append(string.strip(str('%10.3f'%self.ydata[i])))
            amp=(maxi-mini)
            levels=arange(mini-amp*.05,maxi+amp*.05,amp/self.ni)

            cmap=pylab.cm.get_cmap('jet', len(levels)-1)
            
            cset = self.plotdata.contourf(X, Y, Z, levels, cmap=cmap,)

            #cset2 = self.plotdata.contour(X, Y, Z, cset.levels,
            #        colors = 'k',
            #        hold='on')
            self.plotdata.set_xlim(self.ppm[len(self.ppm)-1],self.ppm[0])
            self.plotdata.set_xlabel('Chemical shift (ppm)')
            self.plotdata.set_ylabel('Index\n')
            if self.ydata:
                self.plotdata.set_yticks(yticks)   
                self.plotdata.set_yticklabels(ytickslabel)
                self.plotdata.set_ylabel(self.label+'\n')

            self.plotdata.set_title("POPT 2D display, "+self.title)
##                clabels=[mini,0.,maxi]
##                self.fig.colorbar(cset,clabels=clabels)
##                #self.plotdata.grid(self.grid)

       #update everything and exit
        self.canvas.draw()
        self.toolbar.update()

        self.dlg.Destroy()
        return True

#============================================================================
class PlotChoice(wx.Panel):
#============================================================================    
    def __init__(self, parent):

        wx.Panel.__init__(self, parent, -1)
        self.parent=parent
        sampleList = ['Standard', 'Popt', 'Popt 2D']

        wx.StaticText(self, -1, "Display:", pos=(15, 10), size=(75, -1))
        self.ch = wx.Choice(self, -1, pos=(15, 25), choices = sampleList)
        self.Bind(wx.EVT_CHOICE, self.EvtChoice, self.ch)

        l1 = wx.StaticText(self, -1, "Offset (Std)", pos=(15, 50), size=(75, -1))
        t1 = wx.TextCtrl(self, -1, "0.0",pos=(15, 65), size=(75, -1))
        wx.CallAfter(t1.SetInsertionPoint, 0)
        self.tc1 = t1
        self.Bind(wx.EVT_TEXT, self.EvtText, t1)
        t1.Bind(wx.EVT_CHAR, self.EvtChar)
        
        self.mode="Standard"
        self.offset=0.0
        
    def EvtChoice(self, event):
        self.mode=event.GetString()
        self.parent.plot_data(mode=self.mode, offset=self.offset)

    def EvtText(self, event):
        self.offset=float(event.GetString())
    
    def EvtChar(self, event):
        try:
            if event.GetKeyCode()==13: 
                self.parent.plot_data(mode=self.mode, offset=self.offset)
        except:
            pass
        event.Skip()

#---------------------------------------------------------------------------------    
# Module functions
#---------------------------------------------------------------------------------    
def rtn(st):
#---------------------------------------------------------------------------------    
    """ remove trailing \\n in string read by readline()"""
    st=st.split("\n")[0]
    return st

#---------------------------------------------------------------------------------    
def addlb(ni, spec ,sw ,lb=0 ,gb=0, concentration=0):
#---------------------------------------------------------------------------------    
    """ Make line broadenings """

    spc=[]
    DEBUG_MSG("%5i %13f %s %s" % (ni, sw, str(lb), str(gb)))

    for j in range(ni):        
        npts=len(spec[j])
        fid=[]
        fid = fft(spec[j],None,-1,1)
        #broadening
        lbx=lb[j]
        gbx=gb[j]
        for i in range(npts/2):
             widthh=lbx+(i-1)*gbx*gbx*1.133/sw/2.
             fid[i]=fid[i]*exp(-pi*(i-1)*widthh/sw/2.0)/npts
             fid[npts-i-1]=fid[npts-i-1]*exp(-pi*(i-1)*widthh/sw/2.0)/npts
        spc.append(ifft(fid,None, -1,1))
    return spc
