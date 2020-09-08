import matplotlib.pyplot as plt
import numpy as np
import sys

# plot avg entropy as a function of iteration
avg_ents = []
avg_ents_before = []
middle_ents = []
middle_ents_before = []

qs = [3,4,5,6,7]
d = 20
l = 18


for q in qs:
    filename = 'o18x10_q'+str(q)+'_d'+str(d)+'.txt'
    with open(filename) as f:
        middle_ents = []
        middle_ents_before = []
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
                middle_ents.append(ents[n // 2])

            while "before truncation" not in line and line != '':
                line = f.readline()
            line = f.readline() # this is the data
            ents = line.split(',')[:-1]
            ents = [float(e) for e in ents]
            # print(ents)
            if len(ents) > 0:
                middle_ents_before.append(ents[n // 2])

        avg_ents.append(np.average(middle_ents[2:]))
        avg_ents_before.append(np.average(middle_ents_before[2:]))

delta_ents = [avg_ents_before[i] - avg_ents[i] for i in range(len(avg_ents))]
plt.plot(qs, delta_ents)
plt.xlabel("q")
plt.ylabel("entanglement lost by truncating ()")
# ax[0].title(str(n) + "x" + str(m) + " q=" + str(q))
plt.show()
