from scipy.stats import unitary_group
import numpy as np
import sys

# first system argument is one side of the circuit
# second is m is length of other and the number of iterations
# number of unitaries should be m * (n + 1 + int(2 * (n / 8.0)))
# third is the length of one side of the matrix (should be q^2)
# fourth is the filename to write them to
n = int(sys.argv[1])
m = int(sys.argv[2])
q2 = int(sys.argv[3])
filename = sys.argv[4]

# optional last argument to make all the matrices identities
# just put anything there
if len(sys.argv) > 5:
    id = True
else:
    id = False

num = m * (n + 1 + int(2 * (n / 8.0)))
with open(filename,'w') as f:
    for i in range(num):
        if id:
            xr = np.eye(q2)
            xi = np.zeros((q2, q2))
        else:
            x = unitary_group.rvs(q2)
            xr = x.real
            xi = x.imag
        for j in range(q2):
            for k in range(q2):
                f.write("(%1.10f, %1.10f);" % (xr[j][k],xi[j][k]))
            f.write("\n")
