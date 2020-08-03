import os
import sys
import numpy as np
def getStress (file):
    cmd = "grep 'in kB ' " + file + " | tail -1 | awk '{print $3, $4, $5, $7, $8, $6}'"
    list = os.popen(cmd).read().strip()
    list = np.array(list.split(),dtype=float)
    return (list)
f0 = sys.argv[1]
f1 = sys.argv[2]
ss = float(sys.argv[3])
l0 = getStress (f0)
l1 = getStress (f1)
CC = (l1-l0)/ss/20.

for x in CC:
    sys.stdout.write (" %8.2f" % (x))
sys.stdout.write ("\n")
