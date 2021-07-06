import matplotlib.pyplot as plt
import numpy as np
import sys


# Get all the information, and plot anything vs. anything else as desired

# first dimension is q, second is error, third is iteration (if that's the format)
avg_ents = []
avg_ents_before = []
all_ents = []

middle_ents = []
middle_ents_before = []
all_ents_before = []

overlaps = []
avg_overlaps = []

all_svals = []
all_svals_e = []
all_svals_q = []

es = [.0001,.001,.01,.1,1,2,4,5]
# es = list(range(10, 710, 10))
# q = 7
qs = [2,3,4,5,6,7]
l = 18


for q in qs:
    all_ents_before_q = []
    overlaps_q = []
    avg_overlaps_q = []
    all_ents_q = []
    avg_ents_before_q = []
    avg_ents_q = []
    all_svals_q = []

    for e in es:
        # print(q)
        if e == 1:
            filename = 'o18x10_q'+str(q)+'_e-3.txt'
        elif e == 10:
            filename = 'o18x10_q'+str(q)+'_e-2.txt'
        elif e == .1:
            filename = 'o18x10_q'+str(q)+'_e-4.txt'
        elif e == .01:
            filename = 'o18x10_q'+str(q)+'_e-5.txt'
        elif e == .001:
            filename = 'o18x10_q'+str(q)+'_e-6.txt'
        elif e == .0001:
            filename = 'o18x10_q'+str(q)+'_e-7.txt'
        elif e > 100:
            filename = 'o18x10_q'+str(q)+'_e.'+str(int(e/10))+'.txt'
        else:
            filename = 'o18x10_q'+str(q)+'_e'+str(np.round(e*.001, 3))[1:]+'.txt'
        # filename = 'o18x10_q'+str(q)+'_e-'+str(e)+'.txt'
        with open(filename) as f:
            overlaps_temp = []
            middle_ents = []
            middle_ents_before = []
            svals = []
            all_svals_e = []
            # get the length of the MPS (one side of circuit)
            line = f.readline()
            n = int(line[4:])
            # get the number of iterations (other side of circuit)
            m = int(f.readline()[4:])
            q = int(f.readline()[4:])
            for i in range(m):
                while 'overlap' not in line and line != '':
                    line = f.readline()
                line = f.readline()
                if line == '':
                    break
                overlaps_temp.append(float(line))
                while 'entropies' not in line and line != '':
                    line = f.readline()
                line = f.readline() # this is the data
                ents = line.split(',')[:-1]
                ents = [float(e) for e in ents]
                # print(ents)
                if len(ents) > 0:
                    middle_ents.append(ents[n // 2])

                while 'singular' not in line and line != '':
                    line = f.readline()
                line = f.readline() # this is the data
                svals = line.split(',')[:-1]
                svals = sorted([float(v) for v in svals])[::-1]
                if len(svals) > 0:
                    all_svals_e.append(svals)

                while 'before truncation' not in line and line != '':
                    line = f.readline()
                line = f.readline() # this is the data
                ents = line.split(',')[:-1]
                ents = [float(e) for e in ents]
                # print(ents)
                if len(ents) > 0:
                    middle_ents_before.append(ents[n // 2])
        # all_ents_q.append(middle_ents)
        # avg_ents_q.append(np.nanmean(middle_ents[2:]))
        #
        # all_ents_before_q.append(middle_ents_before)
        # avg_ents_before_q.append(np.nanmean(middle_ents_before[2:]))
        #
        # overlaps_q.append(overlaps_temp)
        # avg_overlaps_q.append(np.nanmean(overlaps_temp[2:]))

        if e <= 500 and e > 100:
            for i in range(2,6):
                filename = 'o18x10_q'+str(q)+'_e.'+str(int(e/10))+'_' + str(i) + '.txt'
                with open(filename) as f:
                    # overlaps_temp = []
                    # middle_ents = []
                    # middle_ents_before = []
                    # get the length of the MPS (one side of circuit)
                    line = f.readline()
                    n = int(line[4:])
                    # get the number of iterations (other side of circuit)
                    m = int(f.readline()[4:])
                    q = int(f.readline()[4:])
                    for i in range(m):
                        while 'overlap' not in line and line != '':
                            line = f.readline()
                        line = f.readline()
                        if line == '':
                            break
                        overlaps_temp.append(float(line))
                        while 'entropies' not in line and line != '':
                            line = f.readline()
                        line = f.readline() # this is the data
                        ents = line.split(',')[:-1]
                        ents = [float(e) for e in ents]
                        # print(ents)
                        if len(ents) > 0:
                            middle_ents.append(ents[n // 2])

                        while 'before truncation' not in line and line != '':
                            line = f.readline()
                        line = f.readline() # this is the data
                        ents = line.split(',')[:-1]
                        ents = [float(e) for e in ents]
                        # print(ents)
                        if len(ents) > 0:
                            middle_ents_before.append(ents[n // 2])


        all_ents_q.append(middle_ents)
        avg_ents_q.append(np.nanmean(middle_ents[2:]))
        # print(np.nanmean(middle_ents[2:]))

        all_ents_before_q.append(middle_ents_before)
        avg_ents_before_q.append(np.nanmean(middle_ents_before[2:]))

        overlaps_q.append(overlaps_temp)
        avg_overlaps_q.append(np.nanmean(overlaps_temp[2:]))
        all_svals_q.append(all_svals_e)

    all_ents.append(all_ents_q)
    avg_ents.append(avg_ents_q)

    all_ents_before.append(all_ents_before_q)
    avg_ents_before.append(avg_ents_before_q)

    overlaps.append(overlaps_q)
    avg_overlaps.append(avg_overlaps_q)

    all_svals.append(all_svals_q)

delta_ents = [[avg_ents_before[i][j] - avg_ents[i][j] for j in range(len(es))] for i in range(len(qs))]
print(np.shape(all_svals))
for i in range(len(qs)):
    plt.plot(np.log(all_svals[i][2][-1]), label="q=" + str(qs[i]))
# plt.plot(es, np.nanmean(middle_ents[2:]))
plt.xlabel("i")
plt.ylabel("log(lambda_i)")
# ax[0].title(str(n) + "x" + str(m) + " q=" + str(q))
plt.legend()
plt.show()
