#!/usr/local/bin/python
# This Python file uses the following encoding: utf-8
import os, sys
import tkinter as tk
from tkinter import ttk
import tkinter.filedialog
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pydicom
import numpy
import math
import imageFunc

def searchImage():
    global fileName, dicomImage, dicomPixelArray
    fileName = tk.filedialog.askopenfilename()
    dicomImage = pydicom.dcmread(fileName)
    dicomPixelArray = dicomImage.pixel_array

def showInformation():
    imageInfoLabel['text'] = imageFunc.dicomInfo(dicomImage)

def applyFilter():
    kernel=selectKernelCBox.get()
    size=selectKernelSizeCBox.get()
    filter=selectFilterCBox.get()
    global dicomFilteredPixelArray
    dicomFilteredPixelArray=imageFunc.applyFilter(kernel, size, filter, dicomPixelArray)
    showFilteredImage()

def showImage():
    for widget in imageFrame.winfo_children():
        widget.destroy()
    figure = plt.Figure()
    subPlot = figure.add_subplot(111)
    subPlot.imshow(dicomPixelArray, cmap=plt.cm.gray)
    imagesTemp = FigureCanvasTkAgg(figure, master=imageFrame)
    imagesTemp.draw()
    imagesTemp.get_tk_widget().pack()

def showFilteredImage():
    for widget in imageFrame.winfo_children():
        widget.destroy()
    figure = plt.Figure()
    subPlot = figure.add_subplot(111)
    subPlot.imshow(dicomFilteredPixelArray, cmap=plt.cm.gray)
    imagesTemp = FigureCanvasTkAgg(figure, master=imageFrame)
    imagesTemp.draw()
    imagesTemp.get_tk_widget().pack()


#fileName = ""
#dicomImage = pydicom.dcmread(fileName)
#dicomPixelArray = dicomImage.pixel_array

filterList=["Sin Filtro","Reduccion","Ignorar","Espejo"]
kernelList=["Promedio","Gaussiano","Medio","Mediano"]
kernelSize=["3x3","5x5","7x7"]

root = tk.Tk()
root.title("Procesamiento de imagenes")
root.configure(bg="black")
root.geometry("1000x500")
root.resizable(0,0)

imageInfoFrame = tk.Frame(root, bg="gray")
imageInfoFrame.pack(side=tk.LEFT, padx= 10, pady=10)
buttonFrame = tk.Frame(root, bg="black")
buttonFrame.pack(side=tk.RIGHT, padx=10, pady=10)

searchImageButton = tk.Button(buttonFrame, text="Search Image", bg="gray", fg="white", height=3, width=15, command=searchImage)
searchImageButton.pack(padx=5, pady=5)
showImageButton = tk.Button(buttonFrame, text="Show Image", bg="gray", fg="white", height=3, width=15, command=showImage)
showImageButton.pack(padx=5, pady=5)
showImageInfoButton = tk.Button(buttonFrame, text="Show Image Info", bg="gray", fg="white", height=3, width=15, command=showInformation)
showImageInfoButton.pack(padx=5, pady=5)
showHistogramButton = tk.Button(buttonFrame, text="Show Histogram", bg="gray", fg="white", height=3, width=15, command = lambda: imageFunc.createHistogram(dicomImage))
showHistogramButton.pack(padx=5, pady=5)
applyFilterButton = tk.Button(buttonFrame, text="Apply Filter", bg="gray",fg="white", height=3, width=15, command=applyFilter)
applyFilterButton.pack(padx=5, pady=5)
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

imageInfoLabel = tk.Label(imageInfoFrame, text="Aqui se desplegara la informacion "+"\n"+"de la imagen seleccionada", bg="gray", fg="white", height=30, width=50)
imageInfoLabel.pack(padx=5, pady=5)

textFrame=tk.Frame(root, bg="black")
textFrame.pack(side=tk.TOP, padx=5, pady=15)

imageFrame=tk.Frame(root, bg="black")
imageFrame.pack()

root.mainloop()
