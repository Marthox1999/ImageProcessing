#!/usr/local/bin/python
# This Python file uses the following encoding: utf-8
import os ,sys
import pydicom
import numpy
import math
import tkinter as tk
import tkinter.filedialog
import imagefunc

#funciones:
#Busca la imagen haciendo uso de un filedialog
def searchImage():
	imageDir = tk.filedialog.askopenfilename()
	global dicomImage
	dicomImage = pydicom.dcmread(imageDir)
#Inserta la informacion de la imagen dentro del cuadro de texto de la interfaz
def loadImage():
	title.insert('end', imagefunc.dicomInfo(ds))
	
#Direccion del archivo Dicom
filename = "MRI Images/MRI01.dcm"
#Lee el archivo .dcm
ds = pydicom.dcmread(filename)
dpa = ds.pixel_array
da = ds.pixel_array

#Aqui comienza la interfaz
interfaz=tk.Tk()
panel = tk.Frame(interfaz, width=800, height=600)
title = tk.Text(panel, width=50, height=20)
searchImageButton = tk.Button(panel, text="Search Image", command = searchImage)
showImageButton = tk.Button(panel, text="Show Image", command = lambda: imagefunc.showImage(da))
showImageInfoButton = tk.Button(panel, text="Show Image Info", command = loadImage)
showHistogramButton = tk.Button(panel, text="Show histogram", command = lambda: imagefunc.createHistogram(ds))

#Empaquetando dentro de la interfaz todos los elementos creados previamente
panel.pack()
title.pack(padx=5,pady=5)
searchImageButton.pack()
showImageButton.pack()
showImageInfoButton.pack()
showHistogramButton.pack()
interfaz.mainloop()