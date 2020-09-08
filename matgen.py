from scipy.stats import unitary_group
import numpy as np
import sys

# first system argument is one side of the circuit
# second is m is length of other and the number of iterations
# number of unitaries should be m * (n + 1 + int(2 * (n / 8.0)))
# but for a cluster state it's m * (2 * n + 1)
# third is the length of one side of the matrix (should be q^2)
# fourth is the filename to write them to
n = int(sys.argv[1])
m = int(sys.argv[2])
# IMPORTANT: q2 is the SQUARE of the local qudit dimension
q2 = int(sys.argv[3])
filename = sys.argv[4]

# optional last argument to make all the matrices identities
# or to have maximally entangling matrices
# 'i' for identities
# 'e' for maximally entangling
id = False
ent = False

if len(sys.argv) > 5:
    if sys.argv[5] == 'i':
        id = True
    elif sys.argv[5] == 'e':
        ent = True

num = m * (2 * n + 1)
with open(filename,'w') as f:
    for i in range(num):
        if id:
            # make the real identity matrix
            xr = np.eye(q2)
            xi = np.zeros((q2, q2))
        elif ent:
            # make a q-dimensional GHZ state creator matrix (hopefully)
            xr = np.zeros((q2, q2))
            xi = np.zeros((q2, q2))
            q = int(np.sqrt(q2))
            negatives = [1] * q
            for j in range(q):
                if j > 0:
                    negatives[-j] = -1
                for k in range(q):
                    for m in range(q):
                        xr[k * q + (k + m) % q, j * q + m] = negatives[k] / np.sqrt(q)

        else:
            x = unitary_group.rvs(q2)
            xr = x.real
            xi = x.imag
        for j in range(q2):
            for k in range(q2):
                f.write("(%1.10f, %1.10f);" % (xr[j][k],xi[j][k]))
            f.write("\n")
