import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import tkinter as tk
from tkinter import *


implSim = h5py.File('output/implSim.hdf', 'r')#Opens the output hdf file in read mode

##Initial output arrays, the 'highlights' if you will: electric field, electron density gradient, machnum gradient, density perturbation
xarr = np.array(implSim['/x'][:]);
zarr = np.array(implSim['/z'][:]);
eden = np.array(implSim['/eden_ncrit'][:]);
field = np.array(implSim['/Original_Field'][:]);
mach = np.array(implSim['/machnum'][:]);
perturb = np.array(implSim['/density_perturbation'][:]);
#convert to microns
xarr*=10000;
zarr*=10000;
fig, axs = plt.subplots(2, 2)
gs = fig.add_gridspec(2, 2, hspace=10, wspace=10000)
cbe = axs[0, 0].contourf(xarr,zarr,eden, 100, cmap='nipy_spectral')
axs[0, 0].set_title('Density Profile')
fig.colorbar(cbe, ax=axs[0,0]);

cbf = axs[0, 1].contourf(xarr,zarr,field, 100, cmap='nipy_spectral')
axs[0, 1].set_title('Field Amplitude')
fig.colorbar(cbf, ax=axs[0,1]);

cbm = axs[1, 0].contourf(xarr,zarr,mach, 100, cmap='nipy_spectral')
axs[1, 0].set_title('Mach Number')
fig.colorbar(cbm, ax = axs[1,0]);

cbp = axs[1, 1].contourf(xarr,zarr,perturb, 100, cmap='nipy_spectral')
axs[1, 1].set_title('Density Perturbation')
fig.colorbar(cbp, ax = axs[1,1]);
#Set the axes labels
for ax in axs.flat:
    ax.set(xlabel='x (\u03BCm)', ylabel='z (\u03BCm)')
fig.tight_layout()
plt.show();
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

#open a Tkinter window
window = tk.Tk()
templist = list(implSim.keys())#Fill a list with the titles of the quantities being plotted
datalist = []
print("Initialized");
for i in templist:
    if i != 'x' and i != 'z':
        datalist.append(i)
#Was used for a previous layout
#window = tk.Frame(master=window,padx=5,pady=5)
#window.grid(row=0,column=0)
#window = tk.Frame(master=window,padx=5,pady=5)
#window.grid(row=0,column=1)

var1 = StringVar(window)
var2 = StringVar(window)
var3 = StringVar(window)
var1.set("") # default value
var2.set("") # default value
var3.set("") # default value
w = OptionMenu(window, var1, *datalist)
v = OptionMenu(window, var2, *datalist)
y = OptionMenu(window, var3, *datalist)

#lblCmpr = tk.Label(
#    master=window,
#    text = "Compare Simulation Results",
#    foreground="white",
#    background="teal",
#    width=30,
#    height=2
#)
lblPlt = tk.Label(
    master=window,
    text = "Plot from C++ Simulation",
    foreground="white",
    background="teal",
    width=30,
    height=2
)
namer = tk.Entry(
    master=window,
    text="Enter Plot Title",
    foreground = "white",
    background = "grey",
    width = 15
)
entryLabel = tk.Label(
    master=window,
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
    master=window,
    text = "Select a dataset",
    foreground="white",
    background="grey",
    width=15,
    height=1
)
dropLabel2 = tk.Label(
    master=window,
    text = "Select a dataset",
    foreground="white",
    background="grey",
    width=15,
    height=1
)
dropLabel3 = tk.Label(
    master=window,
    text = "Select a dataset",
    foreground="white",
    background="grey",
    width=15,
    height=1
)
dim2Button = tk.Button(
    master=window,
    text="2DPlot",
    height=1,
    width=15,
    bg="gray",
    fg="white"
)
#compButton1 = tk.Button(
#    master=window,
#    text="Compare C++ to Matlab",
#    width=20,
#    height=1,
#    bg="gray",
#    fg="white"
#)
#compButton2 = tk.Button(
#    master=window,
#    text="Compare C++ to Yorick",
#    width=20,
#    height=1,
#    bg="gray",
#    fg="white"
#)
#compButton3 = tk.Button(
#    master=window,
#    text="Compare Yorick to Matlab",
#    width=20,
#    height=1,
#    bg="gray",
#    fg="white"
#)
#yorBtn = tk.Button(
#    master=window,
#    text="Plot Yorick Field",
#    width=20,
#    height=1,
#    bg="gray",
#    fg="white"
#)
#matBtn = tk.Button(
#    master=window,
#    text="Plot Matlab Field",
#    width=20,
#    height=1,
#    bg="gray",
#    fg="white"
#)
button = tk.Button(
    master=window,
    text="3DPlot",
    width=15,
    height=1,
    bg="gray",
    fg="white"
)
#compares C++ to Matlab: NOT IN USE
def handle_keypressComp1(event):
    plt.figure();
    path1 = '/z'#+var1.get()
    path2 = '/x'#+var2.get()
    path3 = '/New_Field'
    path4 = '/fieldAmp'
    set1 = np.array(implSim[path1][:])
    set2 = np.array(implSim[path2][:])
    set3 = np.array(implSim[path3][:])
    set4 = np.array(matComp[path4][:])
    for i in range(set3.shape[0]):
        for j in range(set3.shape[1]):
            temp = set4[i][j]
            set4[i][j]= set4[j][i]
            set4[j][i]= temp
    nonzero = np.array(set4);
    nonzero = abs(set4)
    nonzero[nonzero < 1e-6] = 1;
    set5 = 100*(set4-set3)/nonzero
    plt.contourf(set1,set2,set5, 100, cmap='nipy_spectral')
    plt.title("Matlab Field - C++ Field (%)")
    plt.xlabel("z (cm)")
    plt.ylabel("x (cm)")
    plt.colorbar()
    plt.show()
    #ax.plot_surface(set1,set2,temp,cmap=cm.coolwarm,linewidth=0, antialiased=False)
    #plt.show()
#compares C++ to yorick: NOT IN USE
def handle_keypressComp2(event):
    plt.figure();
    path1 = '/z'#+var1.get()
    path2 = '/x'#+var2.get()
    path3 = '/New_Field'
    path4 = '/field'
    set1 = np.array(implSim[path1][:])
    set2 = np.array(implSim[path2][:])
    #set1, set2 = np.meshgrid(set1, set2)
    set3 = np.array(implSim[path3][:])
    set4 = np.array(yorickComp[path4][:])
    nonzero = np.array(set4);
    nonzero = abs(set4)
    nonzero[nonzero < 1e-6] = 1;
    set5 = 100*(set4-set3)/nonzero
    plt.contourf(set1,set2,set5, 100, cmap='plasma')
    plt.title("Yorick Field - C++ Field (%)")
    plt.xlabel("z (cm)")
    plt.ylabel("x (cm)")
    plt.colorbar()
    plt.show()
    #ax.plot_surface(set1,set2,temp,cmap=cm.coolwarm,linewidth=0, antialiased=False)
    #plt.show()
    #compared C++ to matlab: NOT IN USE
def handle_keypressComp3(event):
    plt.figure();
    path1 = '/z'#+var1.get()
    path2 = '/x'#+var2.get()
    path3 = '/field'
    path4 = '/fieldAmp'
    set1 = np.array(implSim[path1][:])
    set2 = np.array(implSim[path2][:])
    #set1, set2 = np.meshgrid(set1, set2)
    set3 = np.array(yorickComp[path3][:])
    set4 = np.array(matComp[path4][:])
    nonzero = abs(set4)
    nonzero[nonzero < 1e-6] = 1;
    set6 = 100*(set4-set3)/nonzero
    plt.contourf(set1,set2,set6, 100, cmap='plasma')
    plt.title("Matlab Field - Yorick Field (%)")
    plt.xlabel("z (cm)")
    plt.ylabel("x (cm)")
    plt.colorbar()
    plt.show()
    #ax.plot_surface(set1,set2,temp,cmap=cm.coolwarm,linewidth=0, antialiased=False)
    #plt.show()
    #plot matlab : NOT IN USE
def pltMatlab(event):
    plt.figure();
    path1 = '/z'#+var1.get()
    path2 = '/x'#+var2.get()
    path3 = '/fieldAmp'
    set1 = np.array(implSim[path1][:])
    set2 = np.array(implSim[path2][:])
    set3 = np.array(matComp[path3][:])
    for i in range(set3.shape[0]):
        for j in range(set3.shape[1]):
            temp = set3[i,j]
            set3[i,j]= set3[j,i]
            set3[j,i]= temp
    plt.contourf(set1,set2,set3, 100, cmap='nipy_spectral')
    plt.title("Matlab Electric Field")
    plt.xlabel("z (cm)")
    plt.ylabel("x (cm)")
    plt.colorbar()
    plt.show()
    #ax.plot_surface(set1,set2,temp,cmap=cm.coolwarm,linewidth=0, antialiased=False)
    #plt.show()
    #plot yorick: NOT IN USE
def pltYorick(event):
    plt.figure();
    path1 = '/z'#+var1.get()
    path2 = '/x'#+var2.get()
    path3 = '/field'
    set1 = np.array(implSim[path1][:])
    set2 = np.array(implSim[path2][:])
    set3 = np.array(yorickComp[path3][:])
    plt.contourf(set1,set2,set3, 100, cmap='plasma')
    plt.title("Yorick Electric Field")
    plt.xlabel("z (cm)")
    plt.ylabel("x (cm)")
    plt.colorbar()
    plt.show()
    #ax.plot_surface(set1,set2,temp,cmap=cm.coolwarm,linewidth=0, antialiased=False)
    #plt.show()

    #plot countour of desired array
def handle_keypress3D(event):
    plt.figure();
    path1 = '/z'
    path2 = '/x'
    path3 = '/'+var1.get()
    set1 = np.array(implSim[path1][:])
    set2 = np.array(implSim[path2][:])
    #set1, set2 = np.meshgrid(set1, set2)
    temp = np.array(implSim[path3][:])
    plt.contourf(set1,set2,temp, 100, cmap='nipy_spectral')
    plt.title(namer.get())
    plt.xlabel("z (cm)")
    plt.ylabel("x (cm)")
    plt.colorbar()
    plt.show()
    #ax.plot_surface(set1,set2,temp,cmap=cm.coolwarm,linewidth=0, antialiased=False)
    #plt.show()

    #plot a line graph of the desired array
def handle_keypress2D(event):
    plt.figure();
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

#Binding buttons to functions
dim2Button.bind("<Button-1>", handle_keypress2D)
button.bind("<Button-1>", handle_keypress3D)
#optionMenu.pack()
lblPlt.pack();
dropLabel1.pack()
w.pack()
namer.pack()
namer.insert(0,"Enter Plot Title")
#dim2Button.pack()
button.pack()

window.mainloop()
