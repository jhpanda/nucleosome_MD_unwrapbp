
p = open('core.ndx','r')
lines = p.readlines()
p.close()

index = []
for line in lines:
    index += [int(s) for s in line.split()]

for idx in index:
    print("%8d"%idx)
