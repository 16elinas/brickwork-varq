# get sum(p_z Tr(rho_z^2)) / sum(p_z)
import matplotlib.pyplot as plt
import numpy as np
import sys

qs = [2, 3, 4, 5]
ls = [20, 25, 30, 35, 40]
zs = []
# number of circuit iterations, sample size
N = 5

# plot avg entropy as a function of iteration
ents_2d = []


for q in qs:
    temp_ents = []
    for L in ls:
        for i in range(N):
            filename = 'o' + str(L) + 'x10_q' + str(q) + '_e-13.txt'

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
                        before_mid_ents.append(ents[(n // 2) - 1])
                        after_mid_ents.append(ents[(n // 2) + 1])

            all_ents = [[middle_ents[i], after_mid_ents[i]] for i in range(len(middle_ents))]
            temp_ents.append(np.average(middle_ents[4:]))
            print(temp_ents)
    ents_2d.append(temp_ents)
print(ents_2d)

for i in range(len(qs)):
    plt.plot(ls, ents_2d[i], label=str(qs[i]))

plt.legend()
plt.show()
