import matplotlib.pyplot as plt
import numpy as np
import sys


# filenames are tL1xL2_qq_out.txt
# ls = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]

# set q
# q = 3

# set L
l = 20

qs = [2, 3, 4, 5, 6, 7, 8]

# whether to do the log or not
log = True

# the array to hold the times
it_times = []
trunc_times = []

avg_it_times = []
avg_trunc_times = []

std_it_times = []
std_trunc_times = []

for q in qs:
    filename = 't' + str(l) + 'x10_q' + str(q) + '_d8_out.txt'
    with open(filename) as f:
        it_times = []
        trunc_times = []
        line = f.readline()
        n = int(line[4:])
        # get the number of iterations (other side of circuit)
        m = int(f.readline()[4:])
        q = int(f.readline()[4:])
        for i in range(m):
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
        else:
            std_it_times.append(np.std(it_times))
            std_trunc_times.append(np.std(trunc_times))
            avg_it_times.append(np.average(it_times))
            avg_trunc_times.append(np.average(trunc_times))


plt.errorbar(qs, avg_it_times, yerr=std_it_times, fmt='.')
plt.xlabel("q")
plt.ylabel("iteration time")
plt.title('Iteration time vs. q, L=' + str(l))
plt.show()

plt.errorbar(qs, avg_trunc_times, yerr=std_trunc_times, fmt='.')
plt.xlabel('q')
plt.ylabel('truncation time')
plt.title('Truncation time vs. q, L=' + str(l))
plt.show()
