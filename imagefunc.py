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

#Filters

#Refleja la imagen para poder hacer la operacion de combolucion en las filas
#y columnas de la imagen original que normamente no se le podria aplicar (bordes)

def combolutionMiror (matriz, kernel):
	#esto es una prueba
	return 0	

#Se reduce el tamaño de la imagen en los bordes a los que no se les pudo aplicar 
#la operacion combolucion
def combolutionReduccion(matriz, kernel):
	#la diferencia entre este y el combolutionIgnore es que este en vez de ignorar los bordes que no
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
def combolutionIgnore(matriz, kernel):
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
	pyplot.imshow(dicomPixelArray)
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
def applyFilter(kernel, filter, dicomPixelArray):
	if(filter == 1):
		filteredDicomPixelArray = combolutionIgnore(dicomPixelArray, kernel)
	if(filter == 2):
		filteredDicomPixelArray = combolutionReduccion(dicomPixelArray, kernel)
	if(filter == 3):
		filteredDicomPixelArray = combolutionMiror(dicomPixelArray, kernel)
	else:
		filteredDicomPixelArray = dicomPixelArray
	return filteredDicomPixelArray
def undoFilter(dicomPixelArray, filteredDicomPixelArray):
	filteredDicomPixelArray = dicomPixelArray.pixel_array
	return dicomPixelArray