from numpy.core.function_base import linspace
import tables
import numpy as np
import matplotlib.pyplot as plt
from tables.table import create_indexes_descr
def compare(filename1: str, filename2: str, fieldname):
    hfile1 = tables.open_file(filename1, mode='r')
    hfile2 = tables.open_file(filename2, mode='r')
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
    print(sh1)
    print(sh2)
    for i in range(min(sh1[0], sh2[0])):
        for j in range(min(sh1[1], sh2[1])):
            if set1[i][j] != set1[i][j] or set2[i][j] != set2[i][j]:
                break
            denom = set2[i][j]
            #if abs(set1[i][j]) < 1e-10:
            #    set1[i][j] = 1
            if denom < 1:
                denom = 1
            compare = abs(set1[i][j]- set2[i][j])
            #if compare > 1:
            #    compare = 0.0;
            comparearr[i, j] = compare
            #print(compare)
            maxval = max(compare, maxval)
    shp = comparearr.shape
    fig1 = plt.figure()#
    print(np.max(comparearr))
    plt.contourf(range(shp[1]), range(shp[0]), comparearr, 100)
    plt.colorbar()
    plt.title("Absolute Intensity Difference ABS(MATLAB1 - MATLAB2)")
    plt.ylabel("Ray Index")
    plt.xlabel("Crossing Index")
    fig1.show()
    fig2 = plt.figure();
    plt.contourf(range(sh1[1]), range(sh1[0]), set1 , 100)
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

    hfile1.close()
    hfile2.close()

if __name__ == '__main__':
    compare('/home/matt/Documents/csc/projects/matlab_clean/matlabcbet2.h5',
            '/home/matt/Documents/csc/projects/matlab_clean/matlabcbet.h5',
            'new_field')