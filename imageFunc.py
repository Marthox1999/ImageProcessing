#!/usr/local/bin/python
# This Python file uses the following encoding: utf-8
import os, sys
import pydicom
import numpy
import math
from matplotlib import pyplot, cm


#Recive a dicom file and return the deseable values of the header

def dicomInfo(dicomImage):
	fullInfo = ""
	try:
		fullInfo += "Largest Image Pixel Value:" + str(dicomImage.LargestImagePixelValue) + "\n"
	except AttributeError:
		fullInfo += "LargestImagePixelValue no esta disponible" + "\n"
	try:
		fullInfo += "Largest Image Pixel: " + str(dicomImage.SmallestImagePixelValue) +"\n"
	except AttributeError:
		fullInfo += "Largest Image Pixel no esta disponible" + "\n"
	try:
		fullInfo += "Manufacturer: " + str(dicomImage.Manufacturer) + "\n"
	except AttributeError:
		fullInfo += "Manufacturer no esta disponible" + "\n"
	try:
		fullInfo += "Rows: " + str(dicomImage.Rows) + "\n"
	except AttributeError:
		fullInfo += "Rows no esta disponible" + "\n"
	try:
		fullInfo += "Columns: " + str(dicomImage.Columns) + "\n"
	except AttributeError:
		fullInfo += "Columns no esta disponible" + "\n"
	try:
		fullInfo += "Patient ID: " + str(dicomImage.PatientID) + "\n"
	except AttributeError:
		fullInfo += "Patiente ID no esta disponible" + "\n"
	try:
		fullInfo += "Series Number: " + str(dicomImage.SeriesNumber) + "\n"
	except AttributeError:
		fullInfo += "Series Number no esta disponible" + "\n"
	try:
		fullInfo += "BitsAllocated: " + str(dicomImage.BitsAllocated) + "\n"
	except AttributeError:
		fullInfo += "BitsAllocated no esta disponible" + "\n"
	try:
		fullInfo += "BitsStored: " + str(dicomImage.BitsStored) + "\n"
	except AttributeError:
		fullInfo += "BitsStored no esta disponible" + "\n"
	try:
		fullInfo += "HighBit " + str(dicomImage.HighBit) + "\n"
	except AttributeError:
		fullInfo += "HighBit no esta disponible" + "\n"
	try:
		fullInfo += "Frecuency: " + str(dicomImage.ImagingFrequency) + "\n"
	except AttributeError:
		fullInfo += "Frecuency no esta disponible" + "\n"
	try:
		fullInfo += "Pixel Bandwidth: " + str(dicomImage.PixelBandwidth) + "\n"
	except AttributeError:
		fullInfo += "Pixel Bandwidth no esta disponible" + "\n"
	try:
		fullInfo += "Pixel Spacing: " + str(dicomImage.PixelSpacing) + "\n"
	except AttributeError:
		fullInfo += "Pixel Spacing no esta disponible" + "\n"
	try:
		fullInfo += "Slice Thickness: " + str(dicomImage.SliceThickness) + "\n"
	except AttributeError:
		fullInfo += "Slice Thickness no esta disponible" + "\n"
	try:
		fullInfo += "Spacing Between Slices: " + str(dicomImage.SpacingBetweenSlices) + "\n"
	except AttributeError:
		fullInfo += "Spacing Between Slices no esta disponible" + "\n"
	return fullInfo
	
#Filters

#Refleja la imagen para poder hacer la operacion de combolucion en las filas
#y columnas de la imagen original que normamente no se le podria aplicar (bordes)
def convolutionMirror(matrix,kernel):
	neighbors = int(math.floor(len(kernel)/2))
	rows, columns = matrix.shape
	mirrorMatrix = numpy.pad(matrix,neighbors,'symmetric')
	newMatrix = numpy.array([[0 for _ in range (columns)] for _ in range (rows)])
	rows, columns = mirrorMatrix.shape
	for i in range (neighbors, rows-neighbors):
		for j in range (neighbors, columns-neighbors):
			newMatrix[i-neighbors][j-neighbors]=numpy.sum(mirrorMatrix[i-neighbors:i+neighbors+1,j-neighbors:j+neighbors+1]*kernel[:,:])
	return newMatrix
def convolutionReduccion(matrix,kernel):
	neighbors = int(math.floor(len(kernel)/2))
	rows, columns = matrix.shape
	newMatrix = numpy.array([[0 for _ in range (columns)] for _ in range (rows)])
	for i in range (neighbors, rows-neighbors-1):
		for j in range (neighbors, columns-neighbors-1):
			newMatrix[i][j]=numpy.sum(matrix[i-neighbors:i+neighbors+1,j-neighbors:j+neighbors+1]*kernel[:,:])
	return newMatrix
def convolutionIgnore(matrix,kernel):
	neighbors = int(math.floor(len(kernel)/2))
	rows, columns = matrix.shape
	newMatrix = numpy.array([[0 for _ in range (columns)] for _ in range (rows)])
	for i in range (neighbors, rows-neighbors-1):
		for j in range (neighbors, columns-neighbors-1):
			newMatrix[i][j]=numpy.sum(matrix[i-neighbors:i+neighbors+1,j-neighbors:j+neighbors+1]*kernel[:,:])
	return newMatrix
#Sobel Filter
def sobelFilter(matrix):
	sobelRigthKernel = numpy.array([[-1,0,1], [-2,0,2], [-1,0,1]])
	sobelDownKernel = numpy.array([[-1,-2,-1], [0,0,0], [1,2,1]])
	neighbors = int(math.floor(len(sobelRigthKernel)/2))
	rows, columns = numpy.shape(matrix)
	newMatrixValues = numpy.array([[0 for _ in range (columns)] for _ in range (rows)])
	newMatrixAngles = numpy.array([[0.0 for _ in range (columns)] for _ in range (rows)])
	horizontalMatrix = numpy.array([[0 for _ in range (columns)] for _ in range (rows)])
	verticalMatrix = numpy.array([[0 for _ in range (columns)] for _ in range (rows)])
	for i in range (neighbors, rows-neighbors-1):
		for j in range (neighbors, columns-neighbors-1):
			horizontalMatrix[i][j]=numpy.sum(matrix[i-neighbors:i+neighbors+1,j-neighbors:j+neighbors+1]*sobelRigthKernel[:,:])
			verticalMatrix[i][j]=numpy.sum(matrix[i-neighbors:i+neighbors+1,j-neighbors:j+neighbors+1]*sobelDownKernel[:,:])
			newMatrixValues[i][j]=abs(horizontalMatrix[i][j])+abs(verticalMatrix[i][j])
			try:
				newMatrixAngles[i][j]=math.atan(verticalMatrix[i][j]/horizontalMatrix[i][j])
			except ZeroDivisionError:
				newMatrixValues[i][j]=0
	horizontalMatrix=None
	verticalMatrix=None
	umbral=thresholdingOtsu(newMatrixValues)
	return newMatrixValues, newMatrixAngles
#Laplacial Filter
def laplacialFilter(matrix):
	laplacianKernel = numpy.array([[-1,-1,-1], [-1,8,-1], [-1,-1,-1]])
	neighbors = int(math.floor(len(laplacianKernel)/2))
	rows, columns = matrix.shape
	newMatrix = numpy.array([[0 for _ in range (columns)] for _ in range (rows)])
	for i in range (neighbors, rows-neighbors-1):
		for j in range (neighbors, columns-neighbors-1):
			newMatrix[i][j]=numpy.sum(matrix[i-neighbors:i+neighbors+1,j-neighbors:j+neighbors+1]*laplacianKernel[:,:])
	return newMatrix
#Crea el histograma de una imagen en formato.dicom
def createHistogram(matrix):
	print(numpy.amax(matrix))
	dicomTotalPixels = len(matrix)*len(matrix[0])
	rows, columns = matrix.shape
	#histogram=[0]*max(map(max, matrix))+1
	#histogram=np.zeros(65536)
	histogram=[0]*(numpy.amax(matrix)+1)
	for i in range (0,rows):
		for j in range (0, columns):
			index = matrix[i][j]
			#print(index)
			histogram[index] += 1
	for i in range (0, len(histogram)-1):
		histogram[i]=histogram[i]/dicomTotalPixels
	return histogram
def ShowHistogram(dicomImage):
	dicomPixelArray = dicomImage.pixel_array
	dicomTotalPixels = len(dicomPixelArray)*len(dicomPixelArray[0])
	try:
		histogram=[0]*dicomImage.LargestImagePixelValue
	except AttributeError:
		histogram=[0]*65536
	for i in range (0,dicomImage.Rows-1):
		for j in range (0, dicomImage.Columns-1):
			index = dicomPixelArray[i][j]
			histogram[index] += 1
	for i in range (0, len(histogram)):
		histogram[i]=histogram[i]/dicomTotalPixels
	pyplot.clf()
	pyplot.plot(histogram)
	pyplot.show()
def thresholdingOtsu(matrix):
	dicomTotalPixels = len(matrix)*len(matrix[0])
	histogram = createHistogram(matrix)
	maxt,q1,q2,m1,m2=0,0,0,0,0
	tArray=[]
	print(len(histogram))
	for t in range (len(histogram)):
		for i in range (0,t):
			q1+=histogram[i]
		for i in range (t+1,len(histogram)):
			q2+=histogram[i]
		for i in range (0,t):
			try:
				m1+=i*histogram[i]/q1
			except ZeroDivisionError:
				m1=0
		for i in range (t+1,len(histogram)):
			try:
				m2+=i*histogram[i]/q1
			except ZeroDivisionError:
				m2=0
		t=math.sqrt((q1*q1)*((m1-m2)**2))
		tArray.append(t)
	maxt=max(tArray)
#Permite seleccionar el filtro que se desa ejecutar
def applyFilter(kernel, size, filter, dicomPixelArray):
	if(filter == "Reduccion"):
		if(kernel == "Promedio" and size == "3x3"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[0][0])
			return filteredDicomPixelArray
		if(kernel == "Promedio" and size == "5x5"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[0][1])
			return filteredDicomPixelArray
		if(kernel == "Promedio" and size == "7x7"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[0][2])
			return filteredDicomPixelArray
		if(kernel == "Gaussiano" and size == "3x3"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[1][0])
			return filteredDicomPixelArray
		if(kernel == "Gaussiano" and size == "5x5"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[1][1])
			return filteredDicomPixelArray
		if(kernel == "Gaussiano" and size == "7x7"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[1][2])
			return filteredDicomPixelArray
		if(kernel == "Rayleigh" and size == "3x3"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[2][0])
			return filteredDicomPixelArray
		if(kernel == "Rayleigh" and size == "5x5"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[2][1])
			return filteredDicomPixelArray
		if(kernel == "Rayleigh" and size == "7x7"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[2][2])
			return filteredDicomPixelArray
		if(kernel == "Mediano" and size == "3x3"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[3][0])
			return filteredDicomPixelArray
		if(kernel == "Mediano" and size == "5x5"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[3][1])
			return filteredDicomPixelArray
		if(kernel == "Mediano" and size == "7x7"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[3][2])
			return filteredDicomPixelArray
	if(filter == "Ignorar"):
		if(kernel == "Promedio" and size == "3x3"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[0][0])
			return filteredDicomPixelArray
		if(kernel == "Promedio" and size == "5x5"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[0][1])
			return filteredDicomPixelArray
		if(kernel == "Promedio" and size == "7x7"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[0][2])
			return filteredDicomPixelArray
		if(kernel == "Gaussiano" and size == "3x3"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[1][0])
			return filteredDicomPixelArray
		if(kernel == "Gaussiano" and size == "5x5"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[1][1])
			return filteredDicomPixelArray
		if(kernel == "Gaussiano" and size == "7x7"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[1][2])
			return filteredDicomPixelArray
		if(kernel == "Rayleigh" and size == "3x3"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[2][0])
			return filteredDicomPixelArray
		if(kernel == "Rayleigh" and size == "5x5"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[2][1])
			return filteredDicomPixelArray
		if(kernel == "Rayleigh" and size == "7x7"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[2][2])
			return filteredDicomPixelArray
		if(kernel == "Mediano" and size == "3x3"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[3][0])
			return filteredDicomPixelArray
		if(kernel == "Mediano" and size == "5x5"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[3][1])
			return filteredDicomPixelArray
		if(kernel == "Mediano" and size == "7x7"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[3][2])
			return filteredDicomPixelArray
	if(filter == "Espejo"):
		if(kernel == "Promedio" and size == "3x3"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[0][0])
			return filteredDicomPixelArray
		if(kernel == "Promedio" and size == "5x5"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[0][1])
			return filteredDicomPixelArray
		if(kernel == "Promedio" and size == "7x7"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[0][2])
			return filteredDicomPixelArray
		if(kernel == "Gaussiano" and size == "3x3"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[1][0])
			return filteredDicomPixelArray
		if(kernel == "Gaussiano" and size == "5x5"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[1][1])
			return filteredDicomPixelArray
		if(kernel == "Gaussiano" and size == "7x7"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[1][2])
			return filteredDicomPixelArray
		if(kernel == "Rayleigh" and size == "3x3"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[2][0])
			return filteredDicomPixelArray
		if(kernel == "Rayleigh" and size == "5x5"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[2][1])
			return filteredDicomPixelArray
		if(kernel == "Rayleigh" and size == "7x7"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[2][2])
			return filteredDicomPixelArray
		if(kernel == "Mediano" and size == "3x3"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[3][0])
			return filteredDicomPixelArray
		if(kernel == "Mediano" and size == "5x5"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[3][1])
			return filteredDicomPixelArray
		if(kernel == "Mediano" and size == "7x7"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[3][2])
			return filteredDicomPixelArray
	else:
		filteredDicomPixelArray = dicomPixelArray
	return filteredDicomPixelArray

#Kernels previamente calculados para la aplicación de filtros

averageKernel3x3=numpy.array([[1,1,1],[1,1,1],[1,1,1]])
averageKernel5x5=numpy.array([[1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1]])
averageKernel7x7=numpy.array([[1,1,1,1,1,1,1],[1,1,1,1,1,1,1],[1,1,1,1,1,1,1],[1,1,1,1,1,1,1],[1,1,1,1,1,1,1],[1,1,1,1,1,1,1],[1,1,1,1,1,1,1]])

gaussianKernel3x3=numpy.array([[1,2,1],[2,4,2],[1,2,1]])
gaussianKernel5x5=numpy.array([[1,4,7,4,1],[4,16,26,16,4],[7,26,41,26,7],[4,16,26,16,4],[1,4,7,4,1]])
gaussianKernel7x7=numpy.array([[0,0,1,2,1,0,0],[0,3,13,22,13,3,0],[1,13,59,97,59,13,1],[2,22,97,159,97,22,2],[1,13,59,97,59,13,1],[0,3,13,22,13,3,0],[0,0,1,2,1,0,0]])

rayleighKernel3x3=numpy.array([[0,0,0],[0,367879,164169],[0,164169,73262]])
rayleighKernel5x5=numpy.array([[0,0,0,0,0],[0,367879,164169,20213,813],[0,164169,73262,9020,363],[0,20213,9020,1110,44],[0,813,363,44,1]])

kernelarray = [
	[averageKernel3x3, averageKernel5x5, averageKernel7x7], #posicion 0 tamaños kernel promedio
	[gaussianKernel3x3, gaussianKernel5x5, gaussianKernel7x7], #posición 1 tamaños kernel gaussiano
	[rayleighKernel3x3,rayleighKernel5x5,None], #posicion 2 tamaños kernel Rayleigh
	[None, None, None], #posicion 3 tamaños kernel Mediano (comming soon)
	]
