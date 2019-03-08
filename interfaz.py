#!/usr/local/bin/python
# This Python file uses the following encoding: utf-8
import os, sys
import tkinter as tk
from tkinter import ttk
import tkinter.filedialog
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from PIL import Image
import pydicom
import numpy
import math
import imageFunc
import copy

fileName = "/home/mateo/Documents/Programing/ImageProcessing/prove_images/cameraM.png"
dicomImage = Image.open(fileName)
dicomPixelArray = numpy.array(dicomImage.convert(mode='L').getdata()).reshape(dicomImage.size[0], dicomImage.size[1])
dicomFilteredPixelArray = None
dicomThresholdingPixelArray = None

def searchImage():
    global fileName, dicomImage, dicomPixelArray
    fileName = tk.filedialog.askopenfilename()
    dicomImage = pydicom.dcmread(fileName)
    dicomPixelArray = dicomImage.pixel_array
    dicomFilteredPixelArray = None
    dicomThresholdingPixelArray = None

def showInformation():
    imageInfoLabel['text'] = imageFunc.dicomInfo(dicomImage)

def applyFilter():
    kernel=selectKernelCBox.get()
    size=selectKernelSizeCBox.get()
    filter=selectFilterCBox.get()
    localDicomPixelArray = dicomPixelArray
    global dicomFilteredPixelArray
    dicomFilteredPixelArray=copy.copy(imageFunc.applyFilter(kernel, size, filter, localDicomPixelArray))
    showFilteredImage()

def thresholding():
    global dicomThresholdingPixelArray
    global dicomAnglesMatrix
    if(dicomFilteredPixelArray is None):
        dicomThresholdingPixelArray, dicomAnglesMatrix = copy.copy(imageFunc.sobelFilter(dicomPixelArray))
    else:
        dicomThresholdingPixelArray,dicomAnglesMatrix=copy.copy(imageFunc.sobelFilter(dicomFilteredPixelArray))
    showThresholdImage()

def showImage():
    for widget in imageFrame.winfo_children():
        widget.destroy()
    figure = plt.Figure()
    subPlot = figure.add_subplot(111)
    subPlot.imshow(dicomPixelArray, cmap=plt.cm.gray)
    imagesTemp = FigureCanvasTkAgg(figure, master=imageFrame)
    imagesTemp.draw()
    imagesTemp.get_tk_widget().pack(padx=5, pady=15)

def showFilteredImage():
    for widget in imageFrame.winfo_children():
        widget.destroy()
    figure = plt.Figure()
    subPlot = figure.add_subplot(111)
    subPlot.imshow(dicomFilteredPixelArray, cmap=plt.cm.gray)
    imagesTemp = FigureCanvasTkAgg(figure, master=imageFrame)
    imagesTemp.draw()
    imagesTemp.get_tk_widget().pack(padx=5, pady=15)

def showThresholdImage():
    for widget in imageFrame.winfo_children():
        widget.destroy()
    figure = plt.Figure()
    subPlot = figure.add_subplot(111)
    subPlot.imshow(dicomThresholdingPixelArray, cmap=plt.cm.gray)
    imagesTemp = FigureCanvasTkAgg(figure, master=imageFrame)
    imagesTemp.draw()
    imagesTemp.get_tk_widget().pack(padx=5, pady=15)

filterList=["Sin Filtro","Reduccion","Ignorar","Espejo"]
kernelList=["Promedio","Gaussiano","Rayleigh","Mediano"]
kernelSize=["3x3","5x5","7x7"]

root = tk.Tk()
root.title("Procesamiento de imagenes")
root.configure(bg="black")
root.geometry("1000x500")
root.resizable(0,0)

imageInfoFrame = tk.Frame(root, bg="black")
imageInfoFrame.pack(side=tk.LEFT, padx= 10, pady=10)
buttonFrame = tk.Frame(root, bg="black")
buttonFrame.pack(side=tk.RIGHT, padx=10, pady=10)

searchImageButton = tk.Button(buttonFrame, text="Search Image", bg="gray30", fg="white", height=2, width=15, command=searchImage)
searchImageButton.pack(padx=5, pady=5)
showImageButton = tk.Button(buttonFrame, text="Show Image", bg="gray30", fg="white", height=2, width=15, command=showImage)
showImageButton.pack(padx=5, pady=5)
showImageInfoButton = tk.Button(buttonFrame, text="Show Image Info", bg="gray30", fg="white", height=2, width=15, command=showInformation)
showImageInfoButton.pack(padx=5, pady=5)
showHistogramButton = tk.Button(buttonFrame, text="Show Histogram", bg="gray30", fg="white", height=2, width=15, command = lambda: imageFunc.ShowHistogram(dicomPixelArray))
showHistogramButton.pack(padx=5, pady=5)
applyFilterButton = tk.Button(buttonFrame, text="Apply Filter", bg="gray30",fg="white", height=2, width=15, command=applyFilter)
applyFilterButton.pack(padx=5, pady=5)
applyBordersButton = tk.Button(buttonFrame, text="Thresholding", bg="gray30",fg="white", height=2, width=15, command=thresholding)
applyBordersButton.pack(padx=5, pady=5)
selectFilterLabel = tk.Label(buttonFrame, text="Seleccione un filtro", bg="black", fg="white")
selectFilterLabel.pack(padx=5, pady=5)
selectFilterCBox = ttk.Combobox(buttonFrame, values=filterList, state="readonly")
selectFilterCBox.current(0)
selectFilterCBox.pack(padx=5, pady=5)
selectKernelLabel = tk.Label(buttonFrame, text="Seleccione un kernel", bg="black", fg="white")
selectKernelLabel.pack(padx=5, pady=5)
selectKernelCBox = ttk.Combobox(buttonFrame, values=kernelList, state="readonly", width=12)
selectKernelCBox.current(0)
selectKernelCBox.pack(side=tk.LEFT, padx=5, pady=5)
selectKernelSizeCBox = ttk.Combobox(buttonFrame, values=kernelSize, state="readonly", width=3)
selectKernelSizeCBox.current(0)
selectKernelSizeCBox.pack(side=tk.RIGHT, padx=5, pady=5)

titleLabel = tk.Label(imageInfoFrame, text="Procesamiento de imagenes", bg="black", fg="white", height=2, width=25, font="Arial 12 bold")
titleLabel.pack(side=tk.TOP, padx=5, pady=5)
imageInfoLabel = tk.Label(imageInfoFrame, text="Aqui se desplegara la informacion "+"\n"+"de la imagen seleccionada", bg="gray", fg="white", height=30, width=35)
imageInfoLabel.pack(padx=5, pady=5)

imageFrame=tk.Frame(root, bg="black")
imageFrame.pack()

root.mainloop()
