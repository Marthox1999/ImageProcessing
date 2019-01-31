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
		fullInfo += "El atributo 'LargestImagePixelValue' no esta disponible" + "\n"
	try:
		fullInfo += "Largest Image Pixel: " + str(dicomImage.SmallestImagePixelValue) +"\n"
	except AttributeError:
		fullInfo += "El atributo Largest Image Pixel no esta disponible" + "\n"
	try:
		fullInfo += "Manufacturer: " + str(dicomImage.Manufacturer) + "\n"
	except AttributeError:
		fullInfo += "El atributo Manufacturer no esta disponible" + "\n"
	try:
		fullInfo += "Rows: " + str(dicomImage.Rows) + "\n"
	except AttributeError:
		fullInfo += "El atributo Rows no esta disponible" + "\n"
	try:
		fullInfo += "Columns: " + str(dicomImage.Columns) + "\n"
	except AttributeError:
		fullInfo += "El atributo Columns no esta disponible" + "\n"
	try:
		fullInfo += "Patient ID: " + str(dicomImage.PatientID) + "\n"
	except AttributeError:
		fullInfo += "El atributo Patiente ID no esta disponible" + "\n"
	try:
		fullInfo += "Series Number: " + str(dicomImage.SeriesNumber) + "\n"
	except AttributeError:
		fullInfo += "El atributo Series Number no esta disponible" + "\n"
	try:
		fullInfo += "BitsAllocated: " + str(dicomImage.BitsAllocated) + "\n"
	except AttributeError:
		fullInfo += "El atributo BitsAllocated no esta disponible" + "\n"
	try:
		fullInfo += "BitsStored: " + str(dicomImage.BitsStored) + "\n"
	except AttributeError:
		fullInfo += "El atributo BitsStored no esta disponible" + "\n"
	try:
		fullInfo += "HighBit " + str(dicomImage.HighBit) + "\n"
	except AttributeError:
		fullInfo += "El atributo HighBit no esta disponible" + "\n"
	try:
		fullInfo += "Frecuency: " + str(dicomImage.ImagingFrequency) + "\n"
	except AttributeError:
		fullInfo += "El atributo Frecuency no esta disponible" + "\n"
	try:
		fullInfo += "Pixel Bandwidth: " + str(dicomImage.PixelBandwidth) + "\n"
	except AttributeError:
		fullInfo += "El atributo Pixel Bandwidth no esta disponible" + "\n"
	try:
		fullInfo += "Pixel Spacing: " + str(dicomImage.PixelSpacing) + "\n"
	except AttributeError:
		fullInfo += "El atributo Pixel Spacing no esta disponible" + "\n"
	try:
		fullInfo += "Slice Thickness: " + str(dicomImage.SliceThickness) + "\n"
	except AttributeError:
		fullInfo += "El atributo Slice Thickness no esta disponible" + "\n"
	try:
		fullInfo += "Spacing Between Slices: " + str(dicomImage.SpacingBetweenSlices) + "\n"
	except AttributeError:
		fullInfo += "El atributo Spacing Between Slices no esta disponible" + "\n"
	return fullInfo

#Kernels

def addNewKernel(kernel):
	kernelarray.append(kernel)
	
#Filters

#Refleja la imagen para poder hacer la operacion de combolucion en las filas
#y columnas de la imagen original que normamente no se le podria aplicar (bordes)
def convolutionMirror (matriz, kernel):
	#esto es una prueba
	return 0	
#Se reduce el tamaño de la imagen en los bordes a los que no se les pudo aplicar 
#la operacion combolucion
def convolutionReduccion(matriz, kernel):
	#la diferencia entre este y el convolutionIgnore es que este en vez de ignorar los bordes que no
	#puede afectar los cambia a 0
	#Prevent the manipulation of the original dicomPixelArraydata
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
def convolutionIgnore(matriz, kernel):
	#Prevent the manipulation of the original dicomPixelArraydata
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
			sumMatrixKernel=0
			for k in range(i-neighbors, i+neighbors):
				for l in range(j-neighbors, j+neighbors):
					sumKernel += kernel[k-i-neighbors][l-j-neighbors]
					sumMatrixKernel += matriz[k][l]*kernel[k-i-neighbors][l-j-neighbors]
			newMatrix[i][j]=sumMatrixKernel/sumKernel
	return newMatrix
#Usando la libreria pyplot muestra la imagen .dicom
def showImage(dicomPixelArray):
	pyplot.clf()
	pyplot.imshow(dicomPixelArray, cmap=pyplot.cm.bone)
	pyplot.show()
#Crea el histograma de una imagen en formato.dicom
def createHistogram(dicomImage):
	dicomPixelArray = dicomImage.pixel_array
	try:
		histogram=[0]*dicomImage.LargestImagePixelValue
	except AttributeError:
		histogram=[0]*65536
	for i in range (0,dicomImage.Rows-1):
		for j in range (0, dicomImage.Columns-1):
			index = dicomPixelArray[i][j]
			histogram[index] += 1
			#histogram.insert(dicomPixelArray[i][j], histogram[dicomPixelArray[i][j]]+1)
	pyplot.clf()
	pyplot.plot(histogram)
	pyplot.show()
#Permite seleccionar el filtro que se desa ejecutar
def applyFilter(kernel, size, filter, dicomPixelArray):
	if(filter == "Reduccion"):
		if(kernel == "Promedio" and size == "3x3"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[0][0])
		if(kernel == "Promedio" and size == "5x5"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[0][1])
		if(kernel == "Promedio" and size == "7x7"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[0][2])
		if(kernel == "Gaussiano" and size == "3x3"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[1][0])
		if(kernel == "Gaussiano" and size == "5x5"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[1][1])
		if(kernel == "Gaussiano" and size == "7x7"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[1][2])
		if(kernel == "Medio" and size == "3x3"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[2][0])
		if(kernel == "Medio" and size == "5x5"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[2][1])
		if(kernel == "Medio" and size == "7x7"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[2][2])
		if(kernel == "Mediano" and size == "3x3"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[3][0])
		if(kernel == "Mediano" and size == "5x5"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[3][1])
		if(kernel == "Mediano" and size == "7x7"):
			filteredDicomPixelArray = convolutionReduccion(dicomPixelArray, kernelarray[3][2])
	if(filter == "Ignorar"):
		if(kernel == "Promedio" and size == "3x3"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[0][0])
		if(kernel == "Promedio" and size == "5x5"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[0][1])
		if(kernel == "Promedio" and size == "7x7"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[0][2])
		if(kernel == "Gaussiano" and size == "3x3"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[1][0])
		if(kernel == "Gaussiano" and size == "5x5"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[1][1])
		if(kernel == "Gaussiano" and size == "7x7"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[1][2])
		if(kernel == "Medio" and size == "3x3"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[2][0])
		if(kernel == "Medio" and size == "5x5"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[2][1])
		if(kernel == "Medio" and size == "7x7"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[2][2])
		if(kernel == "Mediano" and size == "3x3"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[3][0])
		if(kernel == "Mediano" and size == "5x5"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[3][1])
		if(kernel == "Mediano" and size == "7x7"):
			filteredDicomPixelArray = convolutionIgnore(dicomPixelArray, kernelarray[3][2])
	if(filter == "Espejo"):
		if(kernel == "Promedio" and size == "3x3"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[0][0])
		if(kernel == "Promedio" and size == "5x5"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[0][1])
		if(kernel == "Promedio" and size == "7x7"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[0][2])
		if(kernel == "Gaussiano" and size == "3x3"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[1][0])
		if(kernel == "Gaussiano" and size == "5x5"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[1][1])
		if(kernel == "Gaussiano" and size == "7x7"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[1][2])
		if(kernel == "Medio" and size == "3x3"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[2][0])
		if(kernel == "Medio" and size == "5x5"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[2][1])
		if(kernel == "Medio" and size == "7x7"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[2][2])
		if(kernel == "Mediano" and size == "3x3"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[3][0])
		if(kernel == "Mediano" and size == "5x5"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[3][1])
		if(kernel == "Mediano" and size == "7x7"):
			filteredDicomPixelArray = convolutionMirror(dicomPixelArray, kernelarray[3][2])
	else:
		filteredDicomPixelArray = dicomPixelArray
	return filteredDicomPixelArray
def undoFilter(dicomPixelArray, filteredDicomPixelArray):
	filteredDicomPixelArray = dicomPixelArray.pixel_array
	return dicomPixelArray

#Kernels previamente calculados para la aplicación de filtros

averageKernel3x3=[[1,1,1],[1,1,1],[1,1,1]]
averageKernel5x5=[[1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1]]
averageKernel7x7=[[1,1,1,1,1,1,1],[1,1,1,1,1,1,1],[1,1,1,1,1,1,1],[1,1,1,1,1,1,1],[1,1,1,1,1,1,1],[1,1,1,1,1,1,1],[1,1,1,1,1,1,1]]

gaussianKernel3x3=[[1,2,1],[2,4,2],[1,2,1]]
gaussianKernel5x5=[[1,4,7,4,1],[4,16,26,16,4],[7,26,41,26,7],[4,16,26,16,4],[1,4,7,4,1]]
gaussianKernel7x7=[[0,0,1,2,1,0,0],[0,3,13,22,13,3,0],[1,13,59,97,59,13,1],[2,22,97,159,97,22,2],[1,13,59,97,59,13,1],[0,3,13,22,13,3,0],[0,0,1,2,1,0,0]]

kernelarray = [
	[averageKernel3x3, averageKernel5x5, averageKernel7x7], #posicion 0 tamaños kernel promedio
	[gaussianKernel3x3, gaussianKernel5x5, gaussianKernel7x7], #posición 1 tamaños kernel gaussiano
	[], #posicion 2 tamaños kernel Medio (comming soon)
	[], #posicion 3 tamaños kernel Mediano (comming soon)
	[] #posicion 4 tamaños kernel Rayleigh (comming soon)
	]
