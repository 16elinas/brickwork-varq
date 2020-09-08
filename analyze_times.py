import matplotlib.pyplot as plt
import numpy as np
import sys


# filenames are tL1xL2_qq_out.txt
# ls = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]

# set q
# q = 3

# set L
l = 18

qs = [2, 3, 4, 5, 6]

# whether to do the log or not
log = True

cube_root = False

scaleq = False

# the array to hold the times
it_times = []
trunc_times = []

avg_it_times = []
avg_trunc_times = []

std_it_times = []
std_trunc_times = []

avg_ents = []

for q in qs:
    filename = 'o' + str(l) + 'x10_q' + str(q) + '_e-13.txt'
    with open(filename) as f:
        it_times = []
        trunc_times = []
        middle_ents = []
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

            while 'it time' not in line and line != '':
                line = f.readline()

            it_times.append(float(f.readline()))
            line = f.readline() # says truncation time
            trunc_times.append(float(f.readline()))


        if log:
            std_it_times.append(np.std(np.log(it_times)))
            std_trunc_times.append(np.std(np.log(trunc_times)))
            avg_it_times.append(np.log(np.average(it_times)))
            avg_trunc_times.append(np.log(np.average(trunc_times)))
        elif cube_root:
            std_it_times.append(np.std(np.float_power(it_times, 1/6)))
            std_trunc_times.append(np.std(np.float_power(trunc_times, 1/6)))
            avg_it_times.append(np.float_power(np.average(it_times), 1/6))
            avg_trunc_times.append(np.float_power(np.average(trunc_times), 1/6))
        else:
            std_it_times.append(np.std(it_times))
            std_trunc_times.append(np.std(trunc_times))
            avg_it_times.append(np.average(it_times))
            avg_trunc_times.append(np.average(trunc_times))

        if scaleq:
            avg_ents.append(np.average(middle_ents)/q)
        else:
            avg_ents.append(np.average(middle_ents))


fig, ax = plt.subplots(2)
ax[0].errorbar(avg_ents, avg_it_times, yerr=std_it_times, fmt='.')
ax[0].set_xlabel("half-chain entropy")
ax[0].set_ylabel("iteration time (s)")
# ax[0].set_title('Iteration time vs. q, L=' + str(l))
# plt.show()

ax[1].errorbar(avg_ents, avg_trunc_times, yerr=std_trunc_times, fmt='.')
ax[1].set_xlabel('half-chain entropy')
ax[1].set_ylabel('truncation time (s)')
# plt.title('Truncation time vs. q, L=' + str(l))
# plt.show()
plt.show()
