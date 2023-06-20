#! /usr/bin/env python

import pylab
import Pmw
from Tkinter import *
import Tkinter as Tk
import tkFileDialog
import tkMessageBox
import pyfits
import os
import shutil
import shelve
import numpy as np
import subprocess


#import deca

#configdir = '.'
#sys.path.append(configdir)
PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0]
sys.path.append(PATH_TO_PACKAGE)
import deca_setup
import deca

filters = deca.filters()

Current_Dir = os.getcwd()


generalStatus = Tk.StringVar()
generalScaleName = Tk.StringVar()
generalM0Name = Tk.StringVar()
generalGainName = Tk.StringVar()
generalNcombineName = Tk.StringVar()
generalExptimeName = Tk.StringVar()
generalExtName = Tk.StringVar()
generalImageName = Tk.StringVar()
generalMaskName = Tk.StringVar()
generalWeightName = Tk.StringVar()
generalPSFName = Tk.StringVar()
generalRaName = Tk.DoubleVar()
generalDecName = Tk.DoubleVar()
generalNameName = Tk.StringVar()
generalRonName = Tk.StringVar()
generalFWHMName = Tk.StringVar()
generalSky =  Tk.StringVar()
SkySubtr = Tk.IntVar()
generalSampling = Tk.IntVar()
generalDistanceName = Tk.DoubleVar()
generalExtinctionName = Tk.DoubleVar()
generalKCorrName = Tk.DoubleVar()
CoordsCheck = Tk.IntVar()
Filter_index = Tk.IntVar()
#Filter = Tk.StringVar()

ChooseCode = StringVar()
ChooseCode.set(deca_setup.code)

ChooseLum = StringVar()
ChooseLum.set(deca_setup.lum_units)

ChooseGeom = StringVar()
ChooseGeom.set(deca_setup.geom_units)

ChooseSB = StringVar()
ChooseSB.set(deca_setup.SB_units)

nx1 = Tk.IntVar()
ny1 = Tk.IntVar()

#### TODO: REPLACE IT WITH THE SAME FUNCTION FROM DECA
def read_keyw(header,keyw):
    if keyw in header:
      return header[keyw]
    else:
      return 1.

def read_header(galaxy_image):
    hdulist = pyfits.open(galaxy_image)
    prihdr = hdulist[0].header

    scale = read_keyw(prihdr, 'SCALE')
    fwhm = read_keyw(prihdr, 'FWHM')
    gain = read_keyw(prihdr, 'GAIN')
    read_out_noise = read_keyw(prihdr, 'READOUT')
    ncombine = read_keyw(prihdr, 'NCOMBINE')
    exptime = read_keyw(prihdr, 'EXPTIME')
    m0 = read_keyw(prihdr, 'M0')

    #scale = 0.396 #### TODO
    #fwhm = 3. #### TODO

    hdulist.close()
    ADD_INFO = [scale,m0,gain,read_out_noise,ncombine,exptime,fwhm]
    return ADD_INFO
####



def getValuesFromAllFields():   
    inputParamsgeneral = {}
    inputParamsgeneral["image_fits"]=generalImageName.get()
    inputParamsgeneral["weight_fits"]=generalWeightName.get()
    inputParamsgeneral["psf_fits"]=generalPSFName.get()
    inputParamsgeneral["mask_fits"]=generalMaskName.get()

    inputParamsgeneral["magZP"]=generalM0Name.get()    
    inputParamsgeneral["scale"]=generalScaleName.get()
    inputParamsgeneral["GAIN"]=generalGainName.get()
    inputParamsgeneral["NCOMBINE"]=generalNcombineName.get()
    inputParamsgeneral["EXPTIME"]=generalExptimeName.get()
    inputParamsgeneral["READNOISE"]=generalRonName.get()
    inputParamsgeneral["FWHM"]=generalFWHMName.get()
    inputParamsgeneral["Sampling"]=generalSampling.get()
    inputParamsgeneral["Sky"]=generalSky.get()
    inputParamsgeneral["SkySubtr"]=SkySubtr.get()
   
    inputParamsgeneral["RA"]=generalRaName.get()
    inputParamsgeneral["DEC"]=generalDecName.get()
    inputParamsgeneral["NAME"]=generalNameName.get()
    inputParamsgeneral["Distance"]=generalDistanceName.get()
    inputParamsgeneral["Extinction"]=generalExtinctionName.get()
    inputParamsgeneral["Kcorr"]=generalKCorrName.get()
    inputParamsgeneral["CoordsCheck"]=CoordsCheck.get()
    inputParamsgeneral["Filter_index"]=filters.index(Filter.get())
    return inputParamsgeneral
    
    
    
def setValuesToAllFields(inputParamsgeneral):
    if inputParamsgeneral is None:
        return
    

    generalImageName.set(inputParamsgeneral["image_fits"])
    generalWeightName.set(inputParamsgeneral["weight_fits"])
    generalPSFName.set(inputParamsgeneral["psf_fits"])
    generalMaskName.set(inputParamsgeneral["mask_fits"])

    generalM0Name.set(inputParamsgeneral["magZP"])
    generalScaleName.set(inputParamsgeneral["scale"])
    generalGainName.set(inputParamsgeneral["GAIN"])
    generalNcombineName.set(inputParamsgeneral["NCOMBINE"])
    generalExptimeName.set(inputParamsgeneral["EXPTIME"])
    generalRonName.set(inputParamsgeneral["READNOISE"])
    generalFWHMName.set(inputParamsgeneral["FWHM"])
    generalSky.set(inputParamsgeneral["Sky"])
    SkySubtr.set(inputParamsgeneral["SkySubtr"])
    generalSampling.set(inputParamsgeneral["Sampling"])

    generalRaName.set(inputParamsgeneral["RA"])
    generalDecName.set(inputParamsgeneral["DEC"])
    generalNameName.set(inputParamsgeneral["NAME"])
    generalDistanceName.set(inputParamsgeneral["Distance"])   
    generalExtinctionName.set(inputParamsgeneral["Extinction"])
    generalKCorrName.set(inputParamsgeneral["Kcorr"])
    CoordsCheck.set(inputParamsgeneral["CoordsCheck"])
    Filter_index.set(inputParamsgeneral["Filter_index"])
    Filter.selectitem(filters[Filter_index.get()])
  
    

def saveParams(master, params):
    """Store all parameters of galaxy in a file"""
    fileName = tkFileDialog.asksaveasfilename(parent=master,
                                              filetypes=[("Data Base files", "*.inp")],
                                              title="Open file to save parameters")
    if not fileName:
        return
    try:
        os.remove(fileName)
    except OSError:
        pass

    dataBase = shelve.open(fileName)
    dataBase["inputParamsgeneral"] = params

    dataBase.close()

def loadParams(master):
    """Load prevoiusly saved parameters from file"""
    fileName = tkFileDialog.askopenfilename(parent=master,
                                            filetypes=[("Data Base files", "*.inp")],
                                            title="Open file to load parameters")
    if not fileName:
        return None
    dataBase = shelve.open(fileName)
    inputParamsgeneral = dataBase["inputParamsgeneral"]

    dataBase.close()
    return inputParamsgeneral

def change_setup():
  zz = 1

def settings(master):
      
      # Create window where you can choose the model
      window_set = Tk.Toplevel(master)
      window_set.title("Settings")
      window_set.geometry(("600x200"))


      Tk.Label(window_set, text="Choose the code:",anchor='w', font=("Helvetica", 10, "bold"), width=35).grid(column=0, row=2,sticky='w')
      
      GalfitButton = Radiobutton(window_set, text='GALFIT', fg='blue', variable=ChooseCode, value='GALFIT')
      GalfitButton.grid(column=1,row=2, sticky=Tk.W)

      ImfitButton = Radiobutton(window_set, text='IMFIT', fg='blue', variable=ChooseCode, value='IMFIT')
      ImfitButton.grid(column=1,row=3, sticky=Tk.W)

      SkirtButton = Radiobutton(window_set, text='SKIRT', fg='blue', variable=ChooseCode, value='SKIRT')
      SkirtButton.grid(column=1,row=4, sticky=Tk.W)      
      



      Tk.Label(window_set, text="Choose units for",anchor='w', font=("Helvetica", 10, "bold"), width=35).grid(column=0, row=5,sticky='w')
      
      Tk.Label(window_set, text="     - Geometrical parameters:",anchor='w', font=("Helvetica", 10, "bold"), width=35).grid(column=0, row=6,sticky='w')

      PixButton = Radiobutton(window_set, text='pix', fg='blue', variable=ChooseGeom, value='pix')
      PixButton.grid(column=1,row=6, sticky=Tk.W)        

      ArcButton = Radiobutton(window_set, text='arcsec', fg='blue', variable=ChooseGeom, value='arcsec')
      ArcButton.grid(column=2,row=6, sticky=Tk.W)   

      kpcButton = Radiobutton(window_set, text='kpc', fg='blue', variable=ChooseGeom, value='kpc')
      kpcButton.grid(column=3,row=6, sticky=Tk.W)    
      
      
      
      
      Tk.Label(window_set, text="     - Surface brightness:",anchor='w', font=("Helvetica", 10, "bold"), width=35).grid(column=0, row=7,sticky='w')

      SBADUButton = Radiobutton(window_set, text=u'ADU/pix\u00b2', fg='blue', variable=ChooseSB, value='ADU/pix2')
      SBADUButton.grid(column=1,row=7, sticky=Tk.W)        

      SBMagButton = Radiobutton(window_set, text=u'mag/\u25a1\"', fg='blue', variable=ChooseSB, value='mag/arcsec2')
      SBMagButton.grid(column=2,row=7, sticky=Tk.W)   

      LumButton = Radiobutton(window_set, text=u'L\u2609/pc\u00b2', fg='blue', variable=ChooseLum, value='Lsun/pc2')
      LumButton.grid(column=3,row=7, sticky=Tk.W)   




      Tk.Label(window_set, text="     - Luminosity:",anchor='w', font=("Helvetica", 10, "bold"), width=35).grid(column=0, row=8,sticky='w')

      ADUButton = Radiobutton(window_set, text='ADU', fg='blue', variable=ChooseLum, value='ADU')
      ADUButton.grid(column=1,row=8, sticky=Tk.W)        

      MagButton = Radiobutton(window_set, text='mag', fg='blue', variable=ChooseLum, value='mag')
      MagButton.grid(column=2,row=8, sticky=Tk.W)   

      LumButton = Radiobutton(window_set, text=u'L\u2609', fg='blue', variable=ChooseLum, value='Lsun')
      LumButton.grid(column=3,row=8, sticky=Tk.W)   
      
      
      
      ApplyFrame = Tk.Frame(window_set)
      ApplyFrame.grid(column=1, row=9)
      ApplyFrameButton = Tk.Button(window_set, text="Apply", state="normal",command=change_setup, bg="green",font=("Helvetica", 10, "bold"), width=10)
      ApplyFrameButton.grid(column=1, row=9)


    #-------------------------------------------------------------

def beginning(master):
  def loadIniData():
      fileName = tkFileDialog.askopenfilename(parent=master,
				  filetypes=[("DECA-TK input files", ".inp"), ("DECA input file", "deca_input.dat")],
				  title="Open data file")


  menubar = Tk.Menu(master)
  menubar.configure(background='white')
  
  
  fileMenu = Tk.Menu(menubar, tearoff=0)
  fileMenu.add_command(label="Load input data", 
		      command=lambda: setValuesToAllFields(loadParams(master)), 
		      state="normal")
  fileMenu.add_command(label="Save input data",
		      command=lambda: saveParams(master, getValuesFromAllFields()),
		      state="normal")
  
  
  menubar.add_cascade(label="File", menu=fileMenu)

  menubar.add_command(label="Settings", command=lambda: settings(master))

  menubar.add_command(label="Quit", command=master.quit)
  master.config(menu=menubar)

  ############################################################
  #                  End of menu                             #
  ############################################################

class HyperlinkManager:

    def __init__(self, text):

        self.text = text

        self.text.tag_config("hyper", foreground="blue", underline=1)

        self.text.tag_bind("hyper", "<Enter>", self._enter)
        self.text.tag_bind("hyper", "<Leave>", self._leave)
        self.text.tag_bind("hyper", "<Button-1>", self._click)

        self.reset()

    def reset(self):
        self.links = {}

    def add(self, action):
        # add an action to the manager.  returns tags to use in
        # associated text widget
        tag = "hyper-%d" % len(self.links)
        self.links[tag] = action
        return "hyper", tag

    def _enter(self, event):
        self.text.config(cursor="hand2")

    def _leave(self, event):
        self.text.config(cursor="")

    def _click(self, event):
        for tag in self.text.tag_names(CURRENT):
            if tag[:6] == "hyper-":
                self.links[tag]()
                return

def onFrameConfigure(canvas):
	canvas.configure(scrollregion=canvas.bbox("all"))

def help_file1():
  import webbrowser
  root = Tk.Tk()
  root.title("Help: input")
  scrollbar = Scrollbar(root)
  scrollbar.pack(side=RIGHT, fill=Y)

  
  text = Text(root, wrap=WORD, yscrollcommand=scrollbar.set)
  text.pack()
  scrollbar.config(command=text.yview)
  hyperlink = HyperlinkManager(text)

  def click2():
      webbrowser.open_new(r"http://www.mpe.mpg.de/~erwin/code/imfit/")

  def click1():
      webbrowser.open_new(r"https://users.obs.carnegiescience.edu/peng/work/galfit/galfit.html")

  def click3():
      webbrowser.open_new(r"deca.html")

  text.insert(INSERT, "The input values should be given in a way as they are required by the GALFIT and IMFIT codes.")
  text.insert(INSERT, "\nFor more information, see ")
  text.insert(INSERT, "GALFIT",hyperlink.add(click1))
  text.insert(INSERT, " and ")
  text.insert(INSERT, "IMFIT",hyperlink.add(click2))
  text.insert(INSERT, " websites.\nAlso read ")
  text.insert(INSERT, "DECA manual.\n\n",hyperlink.add(click3))

  text.insert(INSERT, u"Image:\t\t\tSpecify the full path to the fits-image.\nWeight:\t\t\tSpecify the full path to the weight-image or leave 'none'. \nPSF:\t\t\tSpecify the full path to the PSF-image or leave 'none'.\nSampling:\t\t\tInput the sampling factor for the given PSF-image (default: 1.0).\nMask:\t\t\tSpecify the full path to the mask-image.\n\nZero-point:\t\t\tDefined as \u03bc [mag/\u25a1\"] = -2.5*lg(ADU) + Zero-point + 5*lg(Scale).\n\t\t\tFor GALFIT input file: magZP = Zero-point - 2.5lg(Exptime).\nScale:\t\t\tPlate scale in arcsec/pixel.\nGain:\t\t\tIn e-/ADU for the whole exposition.\nNcombine:\t\t\tNumber of combined images.\nExptime:\t\t\tExposure time in sec.\nRead-out noise:\tIn e-.\nFWHM:\t\t\tThe full-width-half-maximum in arcsec.\nSky level:\t\t\tBackground sky level in ADU (check if background has been subtracted).\nFilter:\t\t\tIn nm.\n\nCoordinates:\t\t\tRA and DEC in decimal degrees OR x and y in pixels.\nObject name:\t\t\tWhatever you want, usually the first name from NED.\nDistance:\t\t\tIn Mpc.\nGalactic extinction:\t\t\tIn mag.\nK-correction:\t\t\tIn mag.\n")

  


def copy_file(fileName):
	# Copy the file if it is not in the current directory:
	if os.path.exists(fileName):
	  if os.path.dirname(fileName)==Current_Dir:
	    fileName = fileName.split('/')[-1]
	  else:
	    shutil.copy(fileName,fileName.split('/')[-1])
	    fileName = fileName.split('/')[-1]  
	return fileName


      
def inp(p1,input_par_file):
	############################################################
	#                       Right panel                        #
	############################################################

	rightPanel = Tk.Frame(p1)
	rightPanel.pack(side=Tk.LEFT)
	generalImageName.set("none")
	generalWeightName.set("none")
	generalPSFName.set("none")
	generalMaskName.set("none")

	if input_par_file!=None and os.path.exists(input_par_file):
	    dataBase = shelve.open(input_par_file)
	    inputParamsgeneral = dataBase["inputParamsgeneral"]
	    dataBase.close()
	    
	    generalImageName.set(inputParamsgeneral["image_fits"])
	    generalWeightName.set(inputParamsgeneral["weight_fits"])
	    generalPSFName.set(inputParamsgeneral["psf_fits"])
	    generalMaskName.set(inputParamsgeneral["mask_fits"])

	    generalM0Name.set(inputParamsgeneral["magZP"])
	    generalScaleName.set(inputParamsgeneral["scale"])
	    generalGainName.set(inputParamsgeneral["GAIN"])
	    generalNcombineName.set(inputParamsgeneral["NCOMBINE"])
	    generalExptimeName.set(inputParamsgeneral["EXPTIME"])
	    generalRonName.set(inputParamsgeneral["READNOISE"])
	    generalFWHMName.set(inputParamsgeneral["FWHM"])
	    generalSky.set(inputParamsgeneral["Sky"])
	    SkySubtr.set(inputParamsgeneral["SkySubtr"])
	    generalSampling.set(inputParamsgeneral["Sampling"])

	    generalRaName.set(inputParamsgeneral["RA"])
	    generalDecName.set(inputParamsgeneral["DEC"])
	    generalNameName.set(inputParamsgeneral["NAME"])
	    generalDistanceName.set(inputParamsgeneral["Distance"])   
	    generalExtinctionName.set(inputParamsgeneral["Extinction"])
	    generalKCorrName.set(inputParamsgeneral["Kcorr"])
	    CoordsCheck.set(inputParamsgeneral["CoordsCheck"])
	    Filter_index.set(inputParamsgeneral["Filter_index"])
	    #Filter.selectitem(filters[Filter_index.get()])  
	  

	
	def loadImage():
		fileName = tkFileDialog.askopenfilename(parent=rightPanel,
				    filetypes=[("Fits files", "*.fit*"),("All files", "*.*")],
				    title="Open data file")

		try:
		    ADD_INFO = read_header(fileName)
		    [scale,m0,gain,read_out_noise,ncombine,exptime,fwhm] = ADD_INFO

		    generalM0Name.set(str(m0))
		    generalScaleName.set(str(scale))
		    generalGainName.set(str(gain))
		    generalNcombineName.set(str(ncombine))
		    generalExptimeName.set(str(exptime))
		    generalRonName.set(str(read_out_noise))
		    generalFWHMName.set(str(fwhm))
		    
		    generalSky.set(str(0.)) #### TODO: Find center internally in deca

		    hdulist = pyfits.open(fileName) #### TODO: Find center internally in deca
		    scidata = hdulist[0].data
		    ny,nx = np.shape(scidata)
		    xc = nx/2.
		    yc = ny/2.
		    
		    generalRaName.set(str(xc))
		    generalDecName.set(str(yc))
		except:
		  zz = 1
		
		fileName = copy_file(fileName)

		if fileName!='':
		  generalImageName.set(fileName)
		else:
		  generalImageName.set('none')

		return fileName

		
	  
	  
	def loadWeight():
		fileName = tkFileDialog.askopenfilename(parent=rightPanel,
				  filetypes=[("Fits files", "*.fit*"),("All files", "*.*")],
				  title="Open data file")
		  
		fileName = copy_file(fileName)  
		if fileName!='':
		  generalWeightName.set(fileName)
		else:
		  generalWeightName.set('none')

		return fileName	

	def loadPSF():
		fileName = tkFileDialog.askopenfilename(parent=rightPanel,
				  filetypes=[("Fits files", "*.fit*"),("All files", "*.*")],
				  title="Open data file")
		  
		fileName = copy_file(fileName)  
		if fileName!='':
		  generalPSFName.set(fileName)
		else:
		  generalPSFName.set('none')
		return fileName

	def loadMask():
		fileName = tkFileDialog.askopenfilename(parent=rightPanel,
				  filetypes=[("Fits files", "*.fit*"),("All files", "*.*")],
				  title="Open data file")
		  
		fileName = copy_file(fileName)  
		if fileName!='':
		  generalMaskName.set(fileName)
		else:
		  generalMaskName.set('none')
		return fileName

	def load_ds9():
	  if generalImageName.get()!='none' and generalImageName.get()!='':
	    if os.path.exists(generalImageName.get()):
	      try:
		ds9Proc = subprocess.Popen(["ds9", generalImageName.get(),"-scale", "log"])
	      except:
		tkMessageBox.showerror('Error','Could not launch DS9!')
		return 1

	generalStatus.set("disabled")




	#******************************************* GALAXY IMAGES *******************************************:
	generalPanel = Tk.Frame(rightPanel, pady=5)
	generalPanel.grid(column=0, row=1)

	Tk.Label(generalPanel, text="INPUT FITS-FILES:", fg="blue", font=("Helvetica", 13)).grid(column=0, row=0)


	# HELP PANEL
	helpPanel1 = Tk.Frame(generalPanel)
	helpPanel1.grid(column=6, row=0)
	helpPanelButton1 = Tk.Button(helpPanel1, text="Help", state="normal",command=help_file1, fg="green", font=("Helvetica", 12, "bold"), width=15)
	helpPanelButton1.grid(column=6, row=0)


	# GALAXY IMAGE
	ImagePanel = Tk.Frame(generalPanel)
	ImagePanel.grid(column=0, row=1)
	ImagePanelButton = Tk.Button(ImagePanel, text="Load Image", state="normal",command=loadImage, fg="black", font=("Helvetica", 10, "bold"), width=15)
	ImagePanelButton.grid(column=0, row=1)

	generalImageEntry = Tk.Entry(generalPanel, 
		                             textvariable=generalImageName, 
		                             width=20, 
		                             bg="white")

	generalImageEntry.grid(column=1, row=1, sticky=Tk.W)



	# Show DS9
	CreateDS9 = Tk.Frame(generalPanel)
	CreateDS9.grid(column=2, row=1)
	CreateDS9Button = Tk.Button(CreateDS9, text="View in DS9", state="normal",command= load_ds9, fg="green", font=("Helvetica", 10, "bold"))
	CreateDS9Button.grid(column=2, row=1)


	# WEIGHT IMAGE
	WeightPanel = Tk.Frame(generalPanel)
	WeightPanel.grid(column=0, row=2)
	WeightPanelButton = Tk.Button(WeightPanel, text="Load Weight*", state="normal",command=loadWeight, fg="black", font=("Helvetica", 10), width=15)
	WeightPanelButton.grid(column=0, row=2)
	
	generalWeightEntry = Tk.Entry(generalPanel, 
		                             textvariable=generalWeightName, 
		                             width=20, 
		                             bg="white")
	generalWeightEntry.grid(column=1, row=2, sticky=Tk.W)


	# PSF IMAGE
	PSFPanel = Tk.Frame(generalPanel)
	PSFPanel.grid(column=0, row=3)
	PSFPanelButton = Tk.Button(PSFPanel, text="Load PSF*", state="normal",command=loadPSF, fg="black", font=("Helvetica", 10), width=15)
	PSFPanelButton.grid(column=0, row=3)
	
	generalPSFEntry = Tk.Entry(generalPanel, 
		                             textvariable=generalPSFName, 
		                             width=20, 
		                             bg="white")
	generalPSFEntry.grid(column=1, row=3, sticky=Tk.W)


	if input_par_file==None:
	  Filter_index.set(0)
	  generalSampling.set(1)
	  generalM0Name.set("25.0")
	  generalScaleName.set("1.0")
	  generalGainName.set("1.0")
	  generalNcombineName.set("1")
	  generalExptimeName.set("1.0")
	  generalRonName.set("1.0")
	  generalFWHMName.set("1.0")
	  generalSky.set("0.0")
	  SkySubtr.set(1)
	  generalRaName.set(0.)
	  generalDecName.set(0.)
	  CoordsCheck.set(1)
	  generalDistanceName.set(0.)
	  generalExtinctionName.set(0.)
	  generalKCorrName.set(0.)	  
		    
	  
	  

	#Sampling for PSF
	Tk.Label(generalPanel, text="Sampling").grid(column=2, row=3,sticky='e')
	#generalSampling.set(1)
	generalSamplingEntry = Tk.Entry(generalPanel, 
		                             textvariable=generalSampling, 
		                             width=10, 
		                             bg="white")
	generalSamplingEntry.grid(column=3, row=3, sticky=Tk.E)



	# MASK IMAGE
	MaskPanel = Tk.Frame(generalPanel)
	MaskPanel.grid(column=0, row=4)
	MaskPanelButton = Tk.Button(MaskPanel, text="Load Mask*", state="normal",command=loadMask, fg="black", font=("Helvetica", 10), width=15)
	MaskPanelButton.grid(column=0, row=4)
	
	generalMaskEntry = Tk.Entry(generalPanel, 
		                             textvariable=generalMaskName, 
		                             width=20, 
		                             bg="white")
	generalMaskEntry.grid(column=1, row=4, sticky=Tk.W)



	#******************************************* CCD INFORMATION *******************************************:
	Tk.Label(generalPanel, text="OBSERVATION INFO:", fg="blue", font=("Helvetica", 13)).grid(column=0, row=5,sticky='w')

	
	if generalImageName.get()=='none':
		nx1, ny1 =0, 0
	else:
	  hdulist = pyfits.open(generalImageName.get())
	  prihdr = hdulist[0].header
	  inframe = hdulist[0].data
	  try:
		  nx1.set(prihdr['NAXIS1'])
		  ny1.set(prihdr['NAXIS2'])
	  except:
		  nx1, ny1 =inframe.shape[1], inframe.shape[0]


	#Zero-point
	Tk.Label(generalPanel, text=u"Zero-point [mag/\u25a1\"]", font=("Helvetica", 10, "bold")).grid(column=0, row=6,sticky='w')
	#generalM0Name.set("25.0")
	generalM0Entry = Tk.Entry(generalPanel, 
		                             textvariable=generalM0Name, 
		                             width=10, 
		                             bg="white")
	generalM0Entry.grid(column=1, row=6, sticky=Tk.W)


	#Scale
	Tk.Label(generalPanel, text="Scale [arcsec/pix]", font=("Helvetica", 10, "bold")).grid(column=0, row=7,sticky='w')
	#generalScaleName.set("1.0")
	generalScaleEntry = Tk.Entry(generalPanel, 
		                             textvariable=generalScaleName, 
		                             width=10, 
		                             bg="white")
	generalScaleEntry.grid(column=1, row=7, sticky=Tk.W)

	#Gain
	Tk.Label(generalPanel, text="Gain [e-/ADU]", font=("Helvetica", 10, "bold")).grid(column=0, row=8,sticky='w')
	#generalGainName.set("1.0")
	generalGainEntry = Tk.Entry(generalPanel, 
		                             textvariable=generalGainName, 
		                             width=10, 
		                             bg="white")
	generalGainEntry.grid(column=1, row=8, sticky=Tk.W)

	#Number of combined images
	Tk.Label(generalPanel, text="Ncombine", font=("Helvetica", 10, "bold")).grid(column=0, row=9,sticky='w')
	#generalNcombineName.set("1")
	generalNcombineEntry = Tk.Entry(generalPanel, 
		                             textvariable=generalNcombineName, 
		                             width=10, 
		                             bg="white")
	generalNcombineEntry.grid(column=1, row=9, sticky=Tk.W)

	#Exposition time(s)
	Tk.Label(generalPanel, text="Exptime [sec]", font=("Helvetica", 10, "bold")).grid(column=0, row=10,sticky='w')
	#generalExptimeName.set("1.0")
	generalExptimeEntry = Tk.Entry(generalPanel, 
		                             textvariable=generalExptimeName, 
		                             width=10, 
		                             bg="white")
	generalExptimeEntry.grid(column=1, row=10, sticky=Tk.W)

	#Read-out noise
	Tk.Label(generalPanel, text="Read-out noise [e-]*").grid(column=0, row=11,sticky='w')
	#generalRonName.set("1.0")
	generalRonEntry = Tk.Entry(generalPanel, 
		                             textvariable=generalRonName, 
		                             width=10, 
		                             bg="white")
	generalRonEntry.grid(column=1, row=11, sticky=Tk.W)

	#FWHM
	Tk.Label(generalPanel, text="FWHM [arcsec]*").grid(column=0, row=12,sticky='w')
	#generalFWHMName.set("1.0")
	generalFWHMEntry = Tk.Entry(generalPanel, 
		                             textvariable=generalFWHMName, 
		                             width=10, 
		                             bg="white")
	generalFWHMEntry.grid(column=1, row=12, sticky=Tk.W)

	#SKY LEVEL
	Tk.Label(generalPanel, text="Sky level [ADU]", font=("Helvetica", 10, "bold")).grid(column=0, row=13,sticky='w')
	#generalSky.set("0.0")
	generalSkyEntry = Tk.Entry(generalPanel, 
		                             textvariable=generalSky, 
		                             width=10, 
		                             bg="white")
	generalSkyEntry.grid(column=1, row=13, sticky=Tk.W)

	#SkySubtr.set(1)
	SkySubtrButton = Checkbutton(generalPanel, text='Subtr', fg='black', variable=SkySubtr, offvalue=0, onvalue=1)
	SkySubtrButton.grid(column=2,row=13, sticky=Tk.W)

	#Filter
	Tk.Label(generalPanel, text="Filter [nm]*").grid(column=0, row=14,sticky='w')
        
	global Filter
        Filter = Pmw.ComboBox(generalPanel,
                scrolledlist_items = filters,
                dropdown = 1,
                entry_width=10) 
	Filter.grid(column=1, row=14, sticky=Tk.W)

        Filter.selectitem(filters[Filter_index.get()])


	Tk.Label(generalPanel, text="____________").grid(column=0, row=18,sticky='w')
	Tk.Label(generalPanel, text="* Optional", font=("Helvetica", 10, "italic")).grid(column=0, row=19,sticky='w')
	
	
	

	
	
	#******************************************* OBJECT INFORMATION *******************************************:
	Tk.Label(generalPanel, text="OBJECT INFO:", fg="blue", font=("Helvetica", 13)).grid(column=4, row=5,sticky='w')	
	

	#Coordinates: RA
	Tk.Label(generalPanel, text="Coordinates:            RA/x", font=("Helvetica", 10, "bold")).grid(column=4, row=6,sticky='w')
	#generalRaName.set(0.)
	generalRaEntry = Tk.Entry(generalPanel, 
		                             textvariable=generalRaName, 
		                             width=10, 
		                             bg="white")
	generalRaEntry.grid(column=5, row=6, sticky=Tk.W)

	#Coordinates: DEC
	Tk.Label(generalPanel, text="DEC/y", font=("Helvetica", 10, "bold")).grid(column=6, row=6,sticky='e')
	#generalDecName.set(0.)
	generalDecEntry = Tk.Entry(generalPanel, 
		                             textvariable=generalDecName, 
		                             width=10, 
		                             bg="white")
	generalDecEntry.grid(column=7, row=6, sticky=Tk.E)

	
	#CoordsCheck.set(1)
	CoordsButton = Checkbutton(generalPanel, text='Pix', fg='black', variable=CoordsCheck, offvalue=0, onvalue=1)
	CoordsButton.grid(column=8,row=6, sticky=Tk.W)

	#Galaxy name
	Tk.Label(generalPanel, text="Object name*").grid(column=4, row=7,sticky='w')
	generalNameEntry = Tk.Entry(generalPanel, 
		                             textvariable=generalNameName, 
		                             width=10, 
		                             bg="white")
	generalNameEntry.grid(column=5, row=7, sticky=Tk.W)


	#Distance
	#generalDistanceName.set(0.)
	Tk.Label(generalPanel, text="Distance [Mpc]*").grid(column=4, row=8,sticky='w')
	generalDistanceEntry = Tk.Entry(generalPanel, 
		                             textvariable=generalDistanceName, 
		                             width=10, 
		                             bg="white")
	generalDistanceEntry.grid(column=5, row=8, sticky=Tk.W)
	
	
	#Extinction
	generalExtinctionName.set(0.)
	Tk.Label(generalPanel, text="Galactic extinction [mag]*").grid(column=4, row=9,sticky='w')
	generalExtinctionEntry = Tk.Entry(generalPanel, 
		                             textvariable=generalExtinctionName, 
		                             width=10, 
		                             bg="white")
	generalExtinctionEntry.grid(column=5, row=9, sticky=Tk.W)
	
	
	#K-correction
	generalKCorrName.set(0.)
	Tk.Label(generalPanel, text="K-correction [mag]*").grid(column=4, row=10,sticky='w')
	generalKCorrEntry = Tk.Entry(generalPanel, 
		                             textvariable=generalKCorrName, 
		                             width=10, 
		                             bg="white")
	generalKCorrEntry.grid(column=5, row=10, sticky=Tk.W)


	input_files = [generalImageName,generalWeightName,generalPSFName,generalMaskName,generalSampling]
	observation_info = [nx1,ny1,generalM0Name,generalScaleName,generalGainName,generalNcombineName,generalExptimeName,generalRonName,generalFWHMName,generalSky,SkySubtr,Filter]
	object_info = [generalRaName,generalDecName,CoordsCheck,generalNameName,generalDistanceName,generalExtinctionName,generalKCorrName]
	settings = [ChooseCode,ChooseGeom,ChooseSB,ChooseLum]
	
	return input_files,observation_info,object_info,settings

	#return nx1,ny1,generalImageName,generalWeightName,generalPSFName,generalMaskName,generalM0Name,generalScaleName,generalGainName,generalNcombineName,generalExptimeName, generalRaName,generalDecName,generalNameName,generalRonName,generalFWHMName,generalSampling,generalSky,CoordsCheck

