
import sys

with open('massnew.dat','r') as f:
    lines = f.readlines()

    bp_ndx = [[] for s in range(78)]
    for i in range(len(lines)):
        bpi, massi = lines[i].split()
        bpi = int(bpi)
        massi = float(massi)
        if massi > 2:
            old = bpi+1
            if bpi<=39:
                bpi = bpi
                bp_ndx[bpi-1] += [i]
            elif 108<=bpi<=146:
                bpi = bpi-108+39
                bp_ndx[bpi-1] += [i]
            elif 147<=bpi<=185:
                bpi = bpi-147+77
                bp_ndx[bpi-1] += [i]
            elif bpi>=254:
                bpi = bpi-254+39
                bp_ndx[bpi-1] += [i]

    for i in range(78):
        for j in bp_ndx[i]:
            sys.stdout.write("%5d %5d\n"%(i,j))
