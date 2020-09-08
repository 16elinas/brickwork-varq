import matplotlib.pyplot as plt
import numpy as np
import sys

# plot singular values for each of the iterations
all_svals = []

filename = sys.argv[1]

with open(filename) as f:
    # get the length of the MPS (one side of circuit)
    line = f.readline()
    n = int(line[4:])
    # get the number of iterations (other side of circuit)
    m = int(f.readline()[4:])
    q = int(f.readline()[4:])
    err = float(f.readline()[15:])
    for i in range(m):
        while 'singular' not in line and line != '':
            line = f.readline()
        line = f.readline() # this is the data
        svals = line.split(',')[:-1]
        svals = sorted([float(v) for v in svals])[::-1]
        if len(svals) > 0:
            all_svals.append(svals)

# print(all_svals)

fig, ax = plt.subplots(1, len(all_svals), sharey='all')
for i in range(len(all_svals)):
    ax[i].plot(np.log(all_svals[i]))

fig.suptitle("Singular Values - " + str(n) + "x" + str(m) + " q=" + str(q))
plt.show()
