import numpy as np
import sys

factor = float(sys.argv[1])
z, nz = np.loadtxt("NzBlue.txt", unpack=True)
np.savetxt("NzBlue_use.txt", np.transpose([z, factor * nz]))
z, nz = np.loadtxt("NzRed.txt", unpack=True)
np.savetxt("NzRed_use.txt", np.transpose([z, factor * nz]))
