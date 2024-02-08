from numpy.lib.function_base import average
from numpy.lib.shape_base import get_array_prepare
import tables
import numpy as np
import matplotlib.pyplot as plt
import os
import glob


hashattrmap = {'inten': 'intensity', 'width_high_inten': 'width',
               'width_low_inten': 'width', 'order_high_inten': 'order',
               'order': 'order'}


class EntryHolder:
    def __init__(self, intensity, width, order, hashattr):
        self.intensity = intensity
        self.width = width
        self.order = order
        self.maxentries = []
        self.avgentries = []
        self.hashattr = hashattr

    def add_entry(self, maxval, avg):
        self.maxentries.append(maxval)
        self.avgentries.append(avg)

    def get_data(self) -> tuple:
        maxarr = np.array(self.maxentries)
        avgarr = np.array(self.avgentries)
        return (np.average(maxarr), np.std(maxarr),
                np.average(avgarr), np.std(avgarr))

    def __eq__(self, o: object) -> bool:
        return (self.name == o.name and
                getattr(self, self.hashattr) == getattr(o, o.hashattr))

    def __hash__(self) -> int:
        return hash(getattr(self, self.hashattr))

    def __str__(self) -> str:
        return 'Entry: {} = {}'.format(self.hashattr, getattr(self, self.hashattr))

    def __repr__(self) -> str:
        return 'Entry: {} = {}'.format(self.hashattr, getattr(self, self.hashattr))


def compare(filename_other: str,
            fieldname: str,
            filename1: str = None,
            hfile: tables.File = None,
            plot: bool = False):
    if hfile is None:
        hfile1 = tables.open_file(filename1, mode='r')
    else:
        hfile1 = hfile
    hfile2 = tables.open_file(filename_other, mode='r')
    if('/' + fieldname not in hfile1 or '/' + fieldname not in hfile2):
        hfile1.close()
        hfile2.close()
        return

    set1 = np.array(hfile1.root[fieldname])
    set2 = np.array(hfile2.root[fieldname])
    sh1 = np.shape(set1)
    sh2 = np.shape(set2)
    indices = np.zeros_like(sh1)
    for ind1, ind2, i in zip(sh1, sh2, range(len(sh2))):
        indices[i] = min(ind1, ind2)
    comparearr = np.zeros(indices)
    maxval = 0.0
    for i in range(min(sh1[0], sh2[0])):
        for j in range(min(sh1[1], sh2[1])):
            if set1[i][j] != set1[i][j] or set2[i][j] != set2[i][j]:
                break
            denom = set2[i][j]

            compare = min(abs(set1[i][j] - set2[i][j])/denom*100, 5)
            #  if compare > 1:
            #    compare = 0.0;
            comparearr[i, j] = compare
            #  print(compare)
            maxval = max(compare, maxval)
    maxdiff = np.max(comparearr)
    average = np.average(comparearr)
    print('\t{}:{}'.format(maxdiff, average))
    hfile2.close()
    if plot:
        shp = comparearr.shape
        fig1 = plt.figure()
        plt.contourf(range(shp[1]), range(shp[0]), comparearr, 100)
        plt.colorbar()
        plt.title("Relative Intensity Difference % |C++ - MATLAB|/MATLAB")
        plt.ylabel("Ray Index")
        plt.xlabel("Crossing Index")
        fig1.show()
        fig2 = plt.figure()
        plt.contourf(range(sh1[1]), range(sh1[0]), set1, 100)
        plt.colorbar()
        plt.title("C++ Intensity (1e14 W/cm^2)")
        plt.ylabel("Ray Index")
        plt.xlabel("Crossing Index")
        fig2.show()
        fig3 = plt.figure()
        plt.contourf(range(sh2[1]), range(sh2[0]), set2, 100)
        plt.colorbar()
        plt.title("MATLAB Intensity (1e14 W/cm^2)")
        plt.ylabel("Ray Index")
        plt.xlabel("Crossing Index")
        plt.show()
        # input()


def analyze_data(filename):
    data_dict = {} #  Store EntryHolder objects
    name = 'intensity'
    entry = None
    currdict = None
    with open(filename, mode='r') as f:
        all_lines = f.readlines()
        for line in all_lines:
            if line.find('suff=') != -1:
                name = line[5:-1]
                data_dict[name] = {}
                currdict = data_dict[name]
            elif line.find(',') != -1:
                paramlist = line.split(',')
                inten = float(paramlist[0])
                width = float(paramlist[1])
                order = float(paramlist[2])
                entry = EntryHolder(inten, width, order, hashattrmap[name])
                currdict[getattr(entry, entry.hashattr)] = entry
            else:
                paramlist = line.split(':')
                entry.add_entry(float(paramlist[0]), float(paramlist[1]))
    for name in data_dict:
        average_val_list = []
        average_std_list = []
        max_val_list = []
        max_std_list = []
        x_axis = []
        listorder = [max_val_list, max_std_list, average_val_list, average_std_list]
        for x_val in data_dict[name]:
            x_axis.append(x_val)
            entry = data_dict[name][x_val]
            data = entry.get_data()
            for (ls, val) in zip(listorder, data):
                ls.append(val)
        fig1 = plt.figure()
        plt.title('Max CBET Intensity Deviation (%) vs. ' + name)
        if name == 'inten':
            plt.title('Max CBET Intensity Deviation (%) vs. Peak Intensity')
            plt.xlabel('Peak Intensity (W/cm^2)')
        elif name == 'width':
            plt.xlabel('Width: I=1/e*I_peak')
        else:
            plt.xlabel('Beam Profile Order n: exp(-2*(|r/w|^n))')
        plt.ylabel("Max CBET Intensity Deviation (%)")
        plt.errorbar(x_axis, max_val_list, yerr=max_std_list)
        fig1.show()
        fig2 = plt.figure()
        plt.title('Average CBET Intensity Deviation (%) vs. ' + name)
        if name == 'inten':
            plt.title('Average CBET Intensity Deviation (%) vs. Peak Intensity')
            plt.xlabel('Peak Intensity (W/cm^2)')
        elif name == 'width':
            plt.xlabel('Width: I=1/e*I_peak')
        else:
            plt.xlabel('Beam Profile Order n: exp(-2*(|r/w|^n))')
        plt.ylabel("Average CBET Intensity Deviation (%)")
        plt.errorbar(x_axis, average_val_list, yerr=average_std_list)
        fig2.show()
    plt.show()




if __name__ == '__main__':
    """numvars = 10
    names = ['inten', 'width_high_inten', 'width_low_inten',
             'order_high_inten', 'order']
    folder = '/home/matt/Documents/csc/projects/matlab_clean/comparison_output'
    for suffix in names:
        print('suff='+suffix)
        for i in range(numvars):
            files = glob.glob(folder + os.sep + 'matlabcbet_'+suffix +
                              '_' + str(i) + '_*.h5')
            if len(files) == 0:
                continue
            basefile = files.pop(0)
            hfile = tables.open_file(basefile, mode='r')
            intensity = hfile.root['intensity'][0][0]
            width = hfile.root['width'][0][0]
            order = hfile.root['order'][0][0]
            print('{},{},{}'.format(intensity, width, order))
            for ofile in files:
                compare(filename_other=ofile, fieldname='new_field', hfile=hfile)
            hfile.close()"""
    #analyze_data('analysis_output.txt')
    compare(filename_other='/home/matt/Documents/csc/projects/matlab_3b/matlabcbet_3beam.h5', fieldname='new_field', filename1='output/implSim.h5', plot=True)

