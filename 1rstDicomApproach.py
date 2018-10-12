import os
import pydicom
import numpy
import math
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

#Refleja la imagen para poder hacer la operacion de combolucion en las filas
#y columnas de la imagen original que normamente no se le podria aplicar (bordes)

def combolutionMiror (matriz, kernel):
	return 0	

#Se reduce el tamaño de la imagen en los bordes a los que no se les pudo aplicar 
#la operacion combolucion
def combolutionReduccion(matriz, kernel):
	#la diferencia entre este y el combolutionIgnore es que este en vez de ignorar los bordes que no
	#puede afectar los cambia a 0
	#Prevent the manipulation of the original data
	numOfRows=len(matriz[0])
	numOfColumns=len(matriz)
	ksize=len(kernel)
	if(ksize%2 == 0):
		return EOFError()	
	neighbors = math.floor(len(kernel)/2)
	newMatrix=[[0 for _ in range(numOfColumns-(2*neighbors))] for _ in range(numOfRows-(2*neighbors))]
	#go through the matrix avoiding the borders necesaries depending on the neighbors
	for i in range (neighbors,numOfColumns-neighbors):
		for j in range (neighbors, numOfRows-neighbors):
			sumKernel=0
			for k in range(i-neighbors, i+neighbors):
				for l in range(j-neighbors, j+neighbors):
					sumKernel += matriz[k][l]*kernel[k-i-neighbors][l-j-neighbors]
			newMatrix[i][j]=sumKernel
	return newMatrix
#No se hace ningun cambio sobre los bordes de la imagen
def combolutionIgnore(matriz, kernel):
	#Prevent the manipulation of the original data
	newMatrix=matriz
	numOfRows=len(matriz[0])
	numOfColumns=len(matriz)
	ksize=len(kernel)
	if(ksize%2 == 0):
		return EOFError()	
	neighbors = math.floor(len(kernel)/2)
	#go through the matrix avoiding the borders necesaries depending on the neighbors
	for i in range (neighbors,numOfColumns-neighbors):
		for j in range (neighbors, numOfRows-neighbors):
			sumKernel=0
			for k in range(i-neighbors, i+neighbors):
				for l in range(j-neighbors, j+neighbors):
					sumKernel += matriz[k][l]*kernel[k-i-neighbors][l-j-neighbors]
			newMatrix[i][j]=sumKernel
	return newMatrix


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