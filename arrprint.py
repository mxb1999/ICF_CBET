import csv
import numpy as np

beamindices = np.zeros((2,16,11), dtype=int)
cnt = 0
with open("rayIndices.csv") as f:
    reader = csv.reader(f, dialect=csv.excel)
    for r in reader:
        r1 = int(r[0])
        r2 = int(r[1])
        l1 = int(r[2])
        beam1 = int((r1-1 )/ 16)
        beam2 = int((r2-1) / 16)
        if r1 == 32:
            print("{} {} {}".format(r1, r2, l1))
        r1 = r1 - beam1*16 - 1

        beamindices[beam1][r1][l1-1] = r2 - beam2*16
print("{", end='')
for i in range(2):
    for j in range(16):
        for q in range(11):
            print("{},".format(beamindices[i][j][q]), end='')
        if j < 15 or i != 1:
            print()
print("};\n", end='')