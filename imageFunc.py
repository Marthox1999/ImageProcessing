#!/usr/local/bin/python
# This Python file uses the following encoding: utf-8
import os, sys
import pydicom
import numpy
import math
from matplotlib import pyplot, cm
import createfilters
import copy

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
	print("valor maximo imagen filtrada = " + str(numpy.amax(newMatrix)))
	return newMatrix
def convolutionReduccion(matrix,kernel):
	neighbors = int(math.floor(len(kernel)/2))
	rows, columns = matrix.shape
	newMatrix = numpy.array([[0 for _ in range (columns)] for _ in range (rows)])
	for i in range (neighbors, rows-neighbors-1):
		for j in range (neighbors, columns-neighbors-1):
			newMatrix[i][j]=numpy.sum(matrix[i-neighbors:i+neighbors+1,j-neighbors:j+neighbors+1]*kernel[:,:])
	print("valor maximo imagen filtrada = " + str(numpy.amax(newMatrix)))
	return newMatrix
def convolutionIgnore(matrix,kernel):
	neighbors = int(math.floor(len(kernel)/2))
	rows, columns = matrix.shape
	newMatrix = numpy.array([[0 for _ in range (columns)] for _ in range (rows)])
	for i in range (neighbors, rows-neighbors-1):
		for j in range (neighbors, columns-neighbors-1):
			newMatrix[i][j]=numpy.sum(matrix[i-neighbors:i+neighbors+1,j-neighbors:j+neighbors+1]*kernel[:,:])
	print("valor maximo imagen filtrada = " + str(numpy.amax(newMatrix)))
	return newMatrix
#Sobel Filter
def sobelFilter(matrix):
	sobelRigthKernel = numpy.array([[-1,0,1], [-2,0,2], [-1,0,1]])
	sobelDownKernel = numpy.array([[-1,-2,-1], [0,0,0], [1,2,1]])
	rows, columns = numpy.shape(matrix)
	newMatrixValues = numpy.array([[0 for _ in range (columns)] for _ in range (rows)])
	newMatrixAngles = numpy.array([[0 for _ in range (columns)] for _ in range (rows)])
	print("valor maximo antes de sobel = " + str(numpy.amax(matrix)))
	verticalMatrix = convolutionMirror(matrix, sobelDownKernel)
	horizontalMatrix = convolutionMirror(matrix, sobelRigthKernel)
	newMatrixValues=numpy.absolute(verticalMatrix)+numpy.absolute(horizontalMatrix)
	newMatrixValues=newMatrixValues.astype(numpy.int64)
	print("valor maximo despues de sobel = " + str(numpy.amax(newMatrixValues)))
	horizontalMatrix=None
	verticalMatrix=None
	newMatrixValues=thresholdingOtsu(newMatrixValues)
	return newMatrixValues, newMatrixAngles
#Laplacial Filter
def laplacialFilter(matrix):
	laplacianKernel = numpy.array([[-1,-1,-1], [-1,8,-1], [-1,-1,-1]])
	neighbors = int(math.floor(len(laplacianKernel)/2))
	rows, columns = matrix.shape
	newMatrix = numpy.array([[0 for _ in range (columns)] for _ in range (rows)])
	for i in range (neighbors, rows-neighbors-1):
		for j in range (neighbors, columns-neighbors-1):
			newMatrix[i][j]=numpy.sum(matrix[i-neighbors:i+neighbors+1,j-neighbors:j+neighbors+1]*laplacianKernel[:,:])/numpy.sum(laplacialFilter)
	return newMatrix
#Crea el histograma de una imagen en formato.dicom
def createHistogram(matrix):
	print("Creando histograma")
	print("Valor maximo gradiente = " + str(numpy.amax(matrix)))
	#dicomTotalPixels = len(matrix)*len(matrix[0])
	rows, columns = matrix.shape
	#histogram=[0]*max(map(max, matrix))+1
	#histogram=np.zeros(65536)
	histogram=numpy.zeros((numpy.amax(matrix))+1)
	for i in range (0,rows):
		for j in range (0, columns):
			index = matrix[i][j]
			#print(index)
			histogram[index] += 1
	#for i in range (0, len(histogram)-1):
	#	histogram[i]=histogram[i]/dicomTotalPixels
	print("valor maximo histograma = ", str(numpy.amax(histogram)))
	print("tamaño histograma = ", str(len(histogram)))
	#pyplot.clf()
	#pyplot.plot(histogram)
	#pyplot.show()
	return histogram
def ShowHistogram(dicomPixelArray):
	dicomTotalPixels = len(dicomPixelArray)*len(dicomPixelArray[0])
	histogram=[0]*(numpy.amax(dicomPixelArray)+1)
	for i in range (0,len(dicomPixelArray)):
		for j in range (0, len(dicomPixelArray[0])):
			index = dicomPixelArray[i][j]
			histogram[index] += 1
	for i in range (0, len(histogram)):
		histogram[i]=histogram[i]/dicomTotalPixels
	pyplot.clf()
	pyplot.plot(histogram)
	pyplot.show()
def thresholdingOtsu(matrix):
	print("Calculando umbral \n")
	totalPixels = matrix.size
	print(totalPixels)
	histogram = createHistogram(matrix)
	print("Histograma creado")
	sumatory, threshold = 0,0
	sumB, maxt, wB, wF = 0,0,0,0
	for t in range (0,numpy.amax(matrix)):
		sumatory += t * histogram[t]
	for t in range (0,numpy.amax(matrix)):
		wB += histogram[t]
		if(wB == 0):
			continue
		wF = totalPixels - wB
		if(wF == 0):
			break
		sumB += t*histogram[t]
		mB = sumB/wB
		mF = (sumatory-sumB)/wF
		varBetween = wB*wF*(mB - mF)*(mB - mF)
		if(maxt < varBetween):
			maxt = varBetween
			threshold = t
	print("t maximo es " + str(threshold))
	row, column = matrix.shape
	for i in range (0,row):
		for j in range (0,column):
			if matrix[i][j] >= threshold:
				matrix[i][j]=1
			elif matrix[i][j] < threshold:
				matrix[i][j]=0
	return matrix
def dilatation(matrix, struct):
    neighbor = math.floor(int(len(struct))/2)
    rowO, columnO=matrix.shape
    row, column = struct.shape
    newImage = numpy.zeros((rowO,columnO))
    change = numpy.zeros((row,column))
    for i in range (neighbor,rowO-neighbor):
        for j in range (neighbor,columnO-neighbor):
            firsti = i - neighbor
            firstj = j - neighbor  
            endi = i + neighbor + 1
            endj = j + neighbor + 1
            #erosion    
            if ( not numpy.array_equiv(matrix[firsti:endi,firstj:endj], struct[:,:])):
                #pero debo tener en cuenta si en la nueva ya ahi unos marcados en ese sector
                if (numpy.sum(newImage[firsti:endi,firstj:endj]) == 0):
                    newImage[firsti:endi,firstj:endj] = change[:,:]
            else:
                #no debo quitar marcas realizadas
                newImage[i,j] = 1
    return newImage

def erosion(matrix, struct):
    neighbor = math.floor(int(len(struct))/2)
    rowO, columnO = matrix.shape
    row, column = struct.shape
    newImage = numpy.zeros((rowO,columnO))
    change = numpy.ones((row,column))
    for i in range (neighbor,rowO-neighbor):
        for j in range (neighbor,columnO-neighbor):
            firsti = i - neighbor
            firstj = j - neighbor  
            endi = i + neighbor + 1
            endj = j + neighbor + 1
            #dilatation   
            # se esta suponiendo origen centro debo tener en cuenta opcion de centro 
            if ( matrix[i,j] == 1 ):
                newImage[firsti:endi,firstj:endj] = change[:,:]
    return newImage

def dilatacionErosion(matrix, struct):
	newMatrix = copy.copy(matrix)
	newMatrix = thresholdingOtsu(newMatrix)
	matrixStruct = numpy.array([[0]*3]*3)
	if(struct == "+"):
		matrixStruct=structarray[0]
	elif(struct == "x"):
		matrixStruct=structarray[1]
	elif(struct == "-"):
		matrixStruct=structarray[2]
	elif(struct == "|"):
		matrixStruct=structarray[3]
	elif(struct == "□"):
		matrixStruct=structarray[4]
	result =dilatation(newMatrix, matrixStruct)-erosion(newMatrix, matrixStruct)
	return (result)

def rgb(a,b,c):
	return([a,b,c])

def calculateCentroids(matrix, centroids):
	print("Calculando Centroides")
	rows, columns = numpy.shape(matrix)
	#Se almacenaran los pixeles que correspondan a cada centroide
	groups = [[] for i in range (len(centroids))]
	paintGroups = numpy.array([[0 for _ in range (columns)] for _ in range (rows)])
	for i in range (0, rows):
		for j in range (0, columns):
			distance = list(map(lambda x: numpy.absolute(x-matrix[i,j]), centroids))
			minDistaneIndex = distance.index(min(distance))
			groups[minDistaneIndex].append(matrix[i,j])
			paintGroups[i][j]=minDistaneIndex
	newCentroids = [0]*len(centroids)
	for i in range (0, len(centroids)):
		try:
			newCentroids[i] = int(round(numpy.sum(groups[i])/ len(groups[i])))
		except(ZeroDivisionError):
			newCentroids[i] = int(round(sum(groups[i])))
	return centroids==newCentroids, newCentroids, paintGroups

def kmeans(matrix, baseCentroids):
	print("Empezo K-means")
	shouldFinish, newCentroids, paintGroups = calculateCentroids(matrix, baseCentroids)
	while not shouldFinish:
		shouldFinish, newCentroids, paintGroups = calculateCentroids(matrix, newCentroids)
		print(newCentroids)
	return newCentroids, paintGroups

def kmeansIntoImage (matrix, baseCentroids):
	rows, columns = numpy.shape(matrix)
	newCentroids, paintGroups=kmeans(matrix, baseCentroids)
	newMatrix=numpy.zeros(shape=(rows, columns, 3))
	colors=[rgb(0,0,0),rgb(255,255,255),rgb(0,0,255),rgb(0,255,250),rgb(255,0,0),
	rgb(0,255,0),rgb(255,255,0),rgb(255,0,255)]
	print("Pintando Imagen")
	for i in range(rows):
		for j in range(columns):
			for k in range (len(newCentroids)):
				if(paintGroups[i][j]==newCentroids.index(newCentroids[k])):
					newMatrix[i][j]=colors[k]
	return newMatrix 

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

averageKernel3x3=numpy.array(
		[[1,1,1],
		[1,1,1],
		[1,1,1]])
averageKernel5x5=numpy.array(
		[[1,1,1,1,1],
		[1,1,1,1,1],
		[1,1,1,1,1],
		[1,1,1,1,1],
		[1,1,1,1,1]])
averageKernel7x7=numpy.array(
		[[1,1,1,1,1,1,1],
		[1,1,1,1,1,1,1],
		[1,1,1,1,1,1,1],
		[1,1,1,1,1,1,1],
		[1,1,1,1,1,1,1],
		[1,1,1,1,1,1,1],
		[1,1,1,1,1,1,1]])

gaussianKernel3x3, _ = createfilters.get_gaussian_filter(1,1)
gaussianKernel5x5, _ = createfilters.get_gaussian_filter(2,1)
gaussianKernel7x7, _ = createfilters.get_gaussian_filter(3,1)

rayleighKernel3x3, _ = createfilters.get_rayleigh_filter(1,1)
rayleighKernel5x5, _ = createfilters.get_rayleigh_filter(2,1)
rayleighKernel7x7, _ = createfilters.get_rayleigh_filter(3,1)

kernelarray = [
	[averageKernel3x3, averageKernel5x5, averageKernel7x7], #posicion 0 tamaños kernel promedio
	[gaussianKernel3x3, gaussianKernel5x5, gaussianKernel7x7], #posición 1 tamaños kernel gaussiano
	[rayleighKernel3x3,rayleighKernel5x5,None], #posicion 2 tamaños kernel Rayleigh
	[None, None, None], #posicion 3 tamaños kernel Mediano (comming soon)
	]

structarray=[
	numpy.array([[0,1,0],[1,1,1],[0,1,0]]),
	numpy.array([[1,0,1],[0,1,0],[1,0,1]]),
	numpy.array([[0,0,0],[1,1,1],[0,0,0]]),
	numpy.array([[0,1,0],[0,1,0],[0,1,0]]),
	numpy.array([[1,1,1],[1,1,1],[1,1,1]])
]
