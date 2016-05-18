from __future__ import print_function
import numpy as np
from pysnptools.snpreader import Bed

data_dir = '/groups/price/hilary/ibd/data'
bedfile = data_dir+'/1000G.EUR.QC.22'
outfile = bedfile+'.f2snps'

bed = Bed(bedfile)
x = bed.read()
b = np.array([sum(x.val[:,i]) in [2,976] and 1 in x.val[:,i] for i in range(len(x.sid))])
f2snps = x.sid[b]
print('\n'.join(f2snps), file = open(outfile,'w'))
