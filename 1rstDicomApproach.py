import os
import pydicom
import numpy
from matplotlib import pyplot, cm
import tkinter as tk

#Direccion del archivo Dicom
filename = "MRI Images/MRI01.dcm"
#Lee el archivo .dcm
ds = pydicom.dcmread(filename)

#funciones:

def dicomInfo():
	fullInfo = ""
	try:
		fullInfo += "Largest Image Pixel " + str(ds.LargestImagePixelValue) + "\n"
	except AttributeError:
		fullInfo += "Este atributo no esta disponible" + "\n"
	try:
		fullInfo += "Largest Image Pixel " + str(ds.SmallestImagePixelValue) +"\n"
	except AttributeError:
		fullInfo += "Este atributo no esta disponible" + "\n"
	try:
		fullInfo += "Manufacturer " + str(ds.Manufacturer) + "\n"
	except AttributeError:
		fullInfo += "Este atributo no esta disponible" + "\n"
	try:
		fullInfo += "Rows " + str(ds.Rows) + "\n"
	except AttributeError:
		fullInfo += "Este atributo no esta disponible" + "\n"
	try:
		fullInfo += "Columns" + str(ds.Columns) + "\n"
	except AttributeError:
		fullInfo += "Este atributo no esta disponible" + "\n"
	try:
		fullInfo += "Patient ID " + str(ds.PatientID) + "\n"
	except AttributeError:
		fullInfo += "Este atributo no esta disponible" + "\n"
	try:
		fullInfo += "Series Number " + str(ds.SeriesNumber) + "\n"
	except AttributeError:
		fullInfo += "Este atributo no esta disponible" + "\n"
	try:
		fullInfo += "BitsAllocated " + str(ds.BitsAllocated) + "\n"
	except AttributeError:
		fullInfo += "Este atributo no esta disponible" + "\n"
	try:
		fullInfo += "BitsStored " + str(ds.BitsStored) + "\n"
	except AttributeError:
		fullInfo += "Este atributo no esta disponible" + "\n"
	try:
		fullInfo += "HighBit " + str(ds.HighBit) + "\n"
	except AttributeError:
		fullInfo += "Este atributo no esta disponible" + "\n"
	try:
		fullInfo += "Frecuency " + str(ds.ImagingFrequency) + "\n"
	except AttributeError:
		fullInfo += "Este atributo no esta disponible" + "\n"
	try:
		fullInfo += "Pixel Bandwidth " + str(ds.PixelBandwidth) + "\n"
	except AttributeError:
		fullInfo += "Este atributo no esta disponible" + "\n"
	try:
		fullInfo += "Pixel Spacing " + str(ds.PixelSpacing) + "\n"
	except AttributeError:
		fullInfo += "Este atributo no esta disponible" + "\n"
	try:
		fullInfo += "Slice Thickness " + str(ds.SliceThickness) + "\n"
	except AttributeError:
		fullInfo += "Este atributo no esta disponible" + "\n"
	try:
		fullInfo += "Spacing Between Slices " + str(ds.SpacingBetweenSlices) + "\n"
	except AttributeError:
		fullInfo += "Este atributo no esta disponible" + "\n"
	return fullInfo
def showImage():
	pyplot.clf()
	pyplot.imshow(numpy.flipud(ds.pixel_array),cmap=pyplot.cm.gray)
	pyplot.show()
def searchImage():
	imageDir=tk.filedialog.askopenfilename()
	global ds
	ds = pydicom.dcmread(imageDir)
def loadImage():
	title.insert('end', imageInfo)
def createHistogram():
	da = ds.pixel_array
	try:
		histogram=[0]*ds.LargestImagePixelValue
	except AttributeError:
		histogram=[0]*65536
	for i in range (0,ds.Rows-1):
		for j in range (0, ds.Columns-1):
			index = da[i][j]
			histogram[index] += 1
			#histogram.insert(da[i][j], histogram[da[i][j]]+1)
	pyplot.clf()
	pyplot.plot(histogram)
	pyplot.show()


da = ds.pixel_array
imageInfo=dicomInfo()
interfaz=tk.Tk()
panel = tk.Frame(interfaz, width=800, height=600)
title = tk.Text(panel, width=50, height=20)
searchImageButton=tk.Button(panel, text="Search Image", command=searchImage)
showImage=tk.Button(panel, text="Show Image", command=showImage)
showImageInfoButton = tk.Button(panel, text="Show Image Info", command=loadImage)
showHistogramButton = tk.Button(panel, text="Show histogram", command= createHistogram)
panel.pack()
title.pack(padx=5,pady=5)
searchImageButton.pack()
showImage.pack()
showImageInfoButton.pack()
showHistogramButton.pack()
interfaz.mainloop()