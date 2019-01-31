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
import imagefunc

root=tk.Tk()
root.title("Procesamiento de imagenes")
root.configure(bg="black")
root.geometry("700x450")
root.resizable(0,0)

textFrame=tk.Frame(root, bg="black")
textFrame.pack(side=tk.TOP, padx=5, pady=15)

imageLabel=tk.Label(textFrame, text="Imagen original", bg="black", fg="White", font=10)
imageLabel.pack(side=tk.LEFT, padx=75, pady=5)
filteredImageLabel=tk.Label(textFrame, text="Imagen filtrada", bg="black", fg="White", font=10)
filteredImageLabel.pack(side=tk.RIGHT, padx=75, pady=5)

imageFrame=tk.Frame(root, bg="black")
imageFrame.pack()
filteredImageCanvas=tk.Canvas(imageFrame, bg="gray", height=300, width=300)
filteredImageCanvas.pack(side=tk.RIGHT, padx=5, pady=5)
imageCanvas=tk.Canvas(imageFrame, bg="gray", height=300, width=300)
imageCanvas.pack(side=tk.LEFT, padx=5, pady=5)

backButton=tk.Button(text="Volver", padx=15)
backButton.pack(side=tk.LEFT, padx=40, pady=5)

def initWindow():
    root.mainloop()
initWindow()