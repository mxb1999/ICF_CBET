import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import tkinter as tk
from tkinter import *

window = tk.Tk()
implSim = h5py.File('src/Output/implSim.hdf', 'r')
templist = list(implSim.keys())
datalist = []
for i in templist:
    datalist.append(i)
var1 = StringVar(window)
var2 = StringVar(window)
var3 = StringVar(window)
var1.set("") # default value
var2.set("") # default value
var3.set("") # default value
w = OptionMenu(window, var1, *datalist)
v = OptionMenu(window, var2, *datalist)
y = OptionMenu(window, var3, *datalist)
namer = tk.Entry(
    foreground = "white",
    background = "grey",
    width = 15
)
entryLabel = tk.Label(
    text = "Plot Name",
    foreground="white",
    background="grey",
    width=15,
    height=1
)
window.title("Implosion Simulation Plotting")
#can use hexidecimal values as well
#optionMenu = OptionMenu(*datalist)
dropLabel1 = tk.Label(
    text = "Select a dataset",
    foreground="white",
    background="grey",
    width=15,
    height=1
)
dropLabel2 = tk.Label(
    text = "Select a dataset",
    foreground="white",
    background="grey",
    width=15,
    height=1
)
dropLabel3 = tk.Label(
    text = "Select a dataset",
    foreground="white",
    background="grey",
    width=15,
    height=1
)
dim2Button = tk.Button(
    text="2DPlot",
    width=7,
    height=1,
    bg="gray",
    fg="white"
)
button = tk.Button(
    text="3DPlot",
    width=7,
    height=1,
    bg="gray",
    fg="white"
)
def handle_keypress3D(event):
    path1 = '/'+var1.get()
    path2 = '/'+var2.get()
    path3 = '/'+var3.get()
    set1 = np.array(implSim[path1][:])
    set2 = np.array(implSim[path2][:])
    #set1, set2 = np.meshgrid(set1, set2)
    temp = np.array(implSim[path3][:])
    plt.contourf(set1,set2,temp, 100, cmap='plasma')
    plt.title(namer.get())
    plt.xlabel(var1.get())
    plt.ylabel(var2.get())
    plt.colorbar()
    plt.show()
    #ax.plot_surface(set1,set2,temp,cmap=cm.coolwarm,linewidth=0, antialiased=False)
    #plt.show()


def handle_keypress2D(event):
    path1 = '/'+var1.get()
    path2 = '/'+var2.get()
    set1 = np.array(implSim[path1][:])
    set2 = np.array(implSim[path2][:])
    #set1, set2 = np.meshgrid(set1, set2)
    print ("Hello")
    plt.plot(set1,set2)
    plt.title(namer.get())
    plt.xlabel(var1.get())
    plt.ylabel(var2.get())
    plt.show()
    #ax.plot_surface(set1,set2,temp,cmap=cm.coolwarm,linewidth=0, antialiased=False)
    #plt.show()


dim2Button.bind("<Button-1>", handle_keypress2D)

button.bind("<Button-1>", handle_keypress3D)
#optionMenu.pack()
dropLabel1.pack()
w.pack()
dropLabel2.pack()
v.pack()
dropLabel3.pack()
y.pack()
namer.pack()
entryLabel.pack()
dim2Button.pack()
button.pack()

window.mainloop()
