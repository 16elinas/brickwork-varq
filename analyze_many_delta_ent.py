import matplotlib.pyplot as plt
import numpy as np
import sys

# plot avg entropy as a function of iteration
avg_ents = []
avg_ents_before = []
all_ents = []
middle_ents = []
middle_ents_before = []
all_ents_before = []

es = [2,5,8]
q = 5
l = 18


for e in es:
    # print(q)
    filename = 'o18x10_q'+str(q)+'_e-'+str(e)+'.txt'
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
            while 'entropies' not in line and line != '':
                line = f.readline()
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

    all_ents.append(middle_ents)
    avg_ents.append(np.nanmean(middle_ents[2:]))

    all_ents_before.append(middle_ents_before)
    avg_ents_before.append(np.nanmean(middle_ents_before[2:]))

# delta_ents = [avg_ents_before[i] - avg_ents[i] for i in range(len(avg_ents))]
for i in range(len(es)):
    plt.plot(all_ents_before[i], label="1e-" + str(es[i]))
# plt.plot(es, np.nanmean(middle_ents[2:]))
plt.xlabel("iteration")
plt.ylabel("entanglement")
# ax[0].title(str(n) + "x" + str(m) + " q=" + str(q))
plt.legend()
plt.show()
