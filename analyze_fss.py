import matplotlib.pyplot as plt
import numpy as np
import sys



# first dimension is q, second is error, third is iteration (if that's the format)
avg_ents = []
std_ents = []
avg_ents_before = []
all_ents = []
avg_last_ents = []
std_last_ents = []

middle_ents = []
middle_ents_before = []
all_ents_before = []

overlaps = []
avg_overlaps = []

all_svals = []
all_svals_q = []

# es = [.0001,.001,.01,.1,1,2,4,5]
# es += list(range(10, 710, 10))
# q = 7
qs = [2,3,4,5]
# ls = list(range(22, 70, 16))
ls = [22, 30, 38, 46]
ilist = list(range(1,51))
# ilist = [1]
# ls += [62]

for q in qs:
    all_ents_before_q = []
    # overlaps_q = []
    # avg_overlaps_q = []
    all_ents_q = []
    avg_ents_before_q = []
    avg_ents_q = []
    all_svals_q = []
    std_ents_q = []
    avg_last_ents_q = []
    std_last_ents_q = []

    for l in ls:
        ents_l = []
        last_ents_l = []
        for i in ilist:
            # for q=5, have 4, 1-3 and an unnamed one
            if q < 5 or i <= 5:
                filename = 'o' + str(l) + 'x20_q'+str(q)+'_e-6_' + str(i) +'.txt'
            # elif q==5 and i==4:
                # filename = 'o' + str(l) + 'x20_q'+str(q)+'_e-6.txt'
            else:
                continue

            # filename = 'o18x10_q'+str(q)+'_e-'+str(e)+'.txt'
            with open(filename) as f:
                overlaps_temp = []
                middle_ents = []
                middle_ents_before = []
                svals = []
                # get the length of the MPS (one side of circuit)
                line = f.readline()
                n = int(line[4:])
                # get the number of iterations (other side of circuit)
                m = int(f.readline()[4:])
                q = int(f.readline()[4:])
                for i in range(m):
                    # while 'overlap' not in line and line != '':
                        # line = f.readline()
                    # line = f.readline()
                    # if line == '':
                        # break
                    # overlaps_temp.append(float(line))
                    while 'iteration' not in line and line != '':
                        line = f.readline()
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
                        all_svals_q.append(svals)

                    # while 'before truncation' not in line and line != '':
                    #     line = f.readline()
                    # line = f.readline() # this is the data
                    # ents = line.split(',')[:-1]
                    # ents = [float(e) for e in ents]
                    # # print(ents)
                    # if len(ents) > 0:
                    #     middle_ents_before.append(ents[n // 2])
            ents_l += middle_ents[7:]
            last_ents_l.append(middle_ents[-1])

        all_ents_q.append(ents_l)
        avg_ents_q.append(np.nanmean(ents_l))
        avg_last_ents_q.append(np.nanmean(last_ents_l))
        std_ents_q.append(np.std(ents_l))
        std_last_ents_q.append(np.std(last_ents_l))


        # all_ents_before_q.append(middle_ents_before)
        # avg_ents_before_q.append(np.nanmean(middle_ents_before[2:]))

        # overlaps_q.append(overlaps_temp)
        # avg_overlaps_q.append(np.nanmean(overlaps_temp[2:]))

    all_ents.append(all_ents_q)
    avg_ents.append(avg_ents_q)
    std_ents.append(std_ents_q)

    avg_last_ents.append(avg_last_ents_q)
    std_last_ents.append(std_last_ents_q)

    # all_ents_before.append(all_ents_before_q)
    # avg_ents_before.append(avg_ents_before_q)

    # overlaps.append(overlaps_q)
    # avg_overlaps.append(avg_overlaps_q)

    all_svals.append(all_svals_q)



# delta_ents = [[avg_ents_before[i][j] - avg_ents[i][j] for j in range(len(ls))] for i in range(len(qs))]
print(np.transpose(avg_ents))
for i in range(len(ls)):
    plt.errorbar(qs, np.transpose(avg_ents)[i], yerr=np.transpose(std_ents)[i], label="L=" + str(ls[i]), capsize=2)
# plt.plot(es, np.nanmean(middle_ents[2:]))
plt.xlabel("q")
plt.ylabel("S")
# ax[0].title(str(n) + "x" + str(m) + " q=" + str(q))
plt.legend()
plt.show()
