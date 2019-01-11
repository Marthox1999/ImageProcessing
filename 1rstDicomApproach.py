import os
import pydicom
import numpy
import math
import tkinter as tk
import ImageFunc

#funciones:
#Busca la imagen haciendo uso de un filedialog
def searchImage():
	imageDir=tk.filedialog.askopenfilename()
	global dicomImage
	dicomImage = pydicom.dcmread(imageDir)
#Direccion del archivo Dicom
filename = "MRI Images/MRI01.dcm"
#Lee el archivo .dcm
ds = pydicom.dcmread(filename)
dpa = ds.pixel_array
da = ds.pixel_array
imageInfo=ImageFunc.dicomInfo(ds)
interfaz=tk.Tk()
panel = tk.Frame(interfaz, width=800, height=600)
title = tk.Text(panel, width=50, height=20)
searchImageButton=tk.Button(panel, text="Search Image", command=searchImage)
showImage=tk.Button(panel, text="Show Image", command=ImageFunc.showImage)
showImageInfoButton = tk.Button(panel, text="Show Image Info", command=ImageFunc.loadImage(title, imageInfo))
showHistogramButton = tk.Button(panel, text="Show histogram", command= ImageFunc.createHistogram(ds))
panel.pack()
title.pack(padx=5,pady=5)
searchImageButton.pack()
showImage.pack()
showImageInfoButton.pack()
showHistogramButton.pack()
interfaz.mainloop()