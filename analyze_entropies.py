import matplotlib.pyplot as plt
import numpy as np
import sys

# plot avg entropy as a function of iteration
avg_ents = []
middle_ents = []

filename = sys.argv[1]

with open(filename) as f:
    # get the length of the MPS (one side of circuit)
    line = f.readline()
    n = int(line[4:])
    # get the number of iterations (other side of circuit)
    m = int(f.readline()[4:])
    q = int(f.readline()[4:])
    for i in range(m):
        while 'iteration' not in line and line != '':
            line = f.readline()
        f.readline() # this says "entropies"
        line = f.readline() # this is the data
        ents = line.split(',')[:-1]
        ents = [float(e) for e in ents]
        # print(ents)
        if len(ents) > 0:
            avg_ents.append(np.average(ents))
            middle_ents.append(ents[n // 2])

plt.plot(middle_ents)
plt.xlabel("iteration")
plt.ylabel("half-chain entanglement")
plt.title(str(n) + "x" + str(m) + " q=" + str(q))
plt.show()
