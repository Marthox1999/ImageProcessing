#!/usr/local/bin/python
# This Python file uses the following encoding: utf-8
import os, sys
import tkinter as tk
from tkinter import ttk
import tkinter.filedialog
import pydicom
import numpy
import math
import imagefunc

def searchImage():
    global fileName, dicomImage, dicomPixelArray
    fileName = tk.filedialog.askopenfilename()
    dicomImage = pydicom.dcmread(fileName)
    dicomPixelArray = dicomImage.pixel_array

def showInformation():
    imageInfoLabel['text'] = imagefunc.dicomInfo(dicomImage)

fileName = "D:\Programación\procesamientoDeImagenes\MRI Images\MRI01.dcm"
dicomImage = pydicom.dcmread(fileName)
dicomPixelArray = dicomImage.pixel_array

filterList=["Reducción","Ignorar","Espejo"]
kernelList=[]

root = tk.Tk()
root.title("Procesamiento de imagenes")
root.configure(bg="black")
root.geometry("800x500")
root.resizable(0,0)

imageInfoFrame = tk.Frame(root, bg="gray")
imageInfoFrame.pack(side=tk.LEFT, padx= 10, pady=10)
buttonFrame = tk.Frame(root, bg="black")
buttonFrame.pack(side=tk.RIGHT, padx=10, pady=10)

searchImageButton = tk.Button(buttonFrame, text="Search Image", bg="black", fg="white", height=3, width=15, command=searchImage)
searchImageButton.pack(padx=5, pady=5)
showImageButton = tk.Button(buttonFrame, text="Show Image", bg="black", fg="white", height=3, width=15, command = lambda: imagefunc.showImage(dicomPixelArray))
showImageButton.pack(padx=5, pady=5)
showImageInfoButton = tk.Button(buttonFrame, text="Show Image Info", bg="black", fg="white", height=3, width=15, command=showInformation)
showImageInfoButton.pack(padx=5, pady=5)
showHistogramButton = tk.Button(buttonFrame, text="Show Histogram", bg="black", fg="white", height=3, width=15, command = lambda: imagefunc.createHistogram(dicomImage))
showHistogramButton.pack(padx=5, pady=5)
applyFilterButton = tk.Button(buttonFrame, text="Apply Filter", bg="black",fg="white", height=3, width=15)
applyFilterButton.pack(padx=5, pady=5)
selectFilterLabel = tk.Label(buttonFrame, text="Seleccione un filtro", bg="black", fg="white")
selectFilterLabel.pack(padx=5, pady=5)
selectFilterCBox = ttk.Combobox(buttonFrame, values= filterList, state="readonly")
selectFilterCBox.pack(padx=5, pady=5)
selectKernelLabel = tk.Label(buttonFrame, text="Seleccione un kernel", bg="black", fg="white")
selectKernelLabel.pack(padx=5, pady=5)
selectKernelCBox = ttk.Combobox(buttonFrame, values= kernelList, state="readonly")
selectKernelCBox.pack(padx=5, pady=5)

imageInfoLabel = tk.Label(imageInfoFrame, text="Aqui se desplegara la informacion de la imagen seleccionada", bg="gray", fg="white", height=30, width=80)
imageInfoLabel.pack(padx=5, pady=5)

root.mainloop()
