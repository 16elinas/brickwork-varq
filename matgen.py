import sys
from numpy import zeros, eye, sqrt
from scipy.stats import unitary_group

########################################################
## Generate all the random unitary gates (matrices)
## for a circuit
########################################################
# n and m are the same parameters as in sebd.cc (circuit params)
n = int(sys.argv[1])
m = int(sys.argv[2])
# IMPORTANT: q2 is the SQUARE of the local qudit dimension
q2 = int(sys.argv[3])
#  the filename to write the matrices to
filename = sys.argv[4]

# optional last argument (5th one) to make all the matrices identities
# 'i' for identities
# 'e' for maximally entangling (this doesn't work)
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
            xr = eye(q2)
            xi = zeros((q2, q2))
        elif ent:
            # make a q-dimensional GHZ state creator matrix (hopefully)
            # I don't think this worked
            xr = zeros((q2, q2))
            xi = zeros((q2, q2))
            q = int(sqrt(q2))
            negatives = [1] * q
            for j in range(q):
                if j > 0:
                    negatives[-j] = -1
                for k in range(q):
                    for m in range(q):
                        xr[k * q + (k + m) % q, j * q + m] = negatives[k] / sqrt(q)

        else:
            x = unitary_group.rvs(q2)
            xr = x.real
            xi = x.imag
        for j in range(q2):
            for k in range(q2):
                f.write("(%1.10f, %1.10f);" % (xr[j][k],xi[j][k]))
            f.write("\n")
