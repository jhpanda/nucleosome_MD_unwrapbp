
bbDNA = 'C5\''
bbDNA = 'P'

Core = {
           #"C" : range(64,135),
           "C" : range(39,135),
           "D" : range(25,102),
           "E" : range(17,117),
           "F" : range(37,126),
           "G" : range(39,135),
           "H" : range(25,102),
           "I" : range(17,117),
           "J" : range(37,126)
       }
h3h4 = "C D G H".split()
h21 = "E F".split()
h22 = "I J".split()

idxheavyDNA = []
idxbbDNA = []
idxheavyCore = []
idxCaCore = []
idxH3H4Core = []

idxCoreDNA = []


def write_ndx(idxs,fidx):
    p = open('Core.ndx','w')
    for i in range(len(idxs)):
        p.write("[ %s ]\n"%(fidx[i]))
        for j in range(len(idxs[i])):
            p.write("%4d "%(idxs[i][j]))
            if (j+1)%15==0:
                p.write("\n")
        p.write("\n")
    p.close()

p = open('conf0-em2.pdb','r')
lines = p.readlines()
p.close()


atid = 0
nh3 = 0
nh4 = 0
for i in range(len(lines)):
    if lines[i].startswith("ATOM  "):
        atid += 1
        atname = lines[i][11:16].strip()
        resname = lines[i][17:21].strip()
        chain = lines[i][21]
        resi = int(lines[i][22:27])

        if chain in ['A','B']:
            if atname == bbDNA:
                idxbbDNA += [atid]
                if resi>-21 and resi<21:
                    idxCoreDNA += [atid]
            if not atname.startswith('H') or atname[0].isdigit():
                idxheavyDNA += [atid]

        else:
            if resi in Core[chain]:
                if atname == 'CA':
                    idxCaCore += [atid]
                    if chain in h3h4:
                        idxH3H4Core += [atid]
                        if chain in ['C']:
                            nh3 += 1
                        if chain in ['G']:
                            nh4 += 1
                if not atname.startswith('H') or atname[0].isdigit():
                    idxheavyCore += [atid]

#idxs = [idxheavyDNA,idxbbDNA,idxheavyCore,idxCaCore,idxH3H4Core]
#fidx = ['heavyDNA','bbDNA','heavyCore','CaCore','H3H4Core']

idxs = [idxheavyDNA,idxbbDNA,idxheavyCore,idxCaCore,idxH3H4Core,idxCoreDNA]
fidx = ['heavyDNA','bbDNA','heavyCore','CaCore','H3H4Core','CoreDNA']
write_ndx(idxs,fidx)
print(nh3, nh4)
