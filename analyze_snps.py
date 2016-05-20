from __future__ import print_function, division
from pysnptools.snpreader import Bed
import numpy as np
import pandas as pd
import itertools as it
import argparse
import pyutils.iter as iter
import find_ibd_segment as fis

def main(args):
    print('reading seeed snps')
    seed_snps = pd.read_csv(args.seed_snps, header=None, names=['SNP'], index_col='SNP')
    seed_snps['ibs_length'] = 0
    seed_snps['ibd'] = 0

    print('reading typed snps')
    typed_snps = pd.read_csv(args.typed_snps, header=None, names=['SNP'])

    print('reading genotypes')
    data = Bed(args.bfile)
    X = data.read().val
    typed_snps_indices = np.sort(data.sid_to_index(typed_snps.SNP))
    typed_snps_bp = data.col_property[typed_snps_indices,2]

    print(len(seed_snps), 'snps in list')
    print(data.iid_count, data.sid_count, 'are dimensions of X')

    def analyze_snp(i):
        # find first typed snp after query snp
        snp_bp = data.col_property[i,2]
        v = np.where(typed_snps_bp > snp_bp)[0]
        if len(v) > 0:
            typed_i = v[0]
        else:
            typed_i = len(typed_snps_indices)-1

        n1, n2 = np.where(X[:,i] == 1)[0]
        if (X[n1,typed_snps_indices[typed_i]] - X[n2, typed_snps_indices[typed_i]])**2 == 4:
            return 0, 0

        typed_il, typed_ir = fis.find_boundaries(
                X[n1,typed_snps_indices],
                X[n2,typed_snps_indices],
                typed_i)
        typed_ir -= 1

        il = typed_snps_indices[typed_il]
        ir = typed_snps_indices[typed_ir]
        cM = data.col_property[ir, 1] - \
                data.col_property[il, 1]
        ibd = (np.mean(X[n1,il:ir] == X[n2,il:ir]) > 0.99)
        return cM, int(ibd)

    for (i, snp) in iter.show_progress(
            it.izip(data.sid_to_index(seed_snps.index), seed_snps.index),
            total=len(seed_snps)):
            # total=10):
        seed_snps.ix[snp, ['ibs_length', 'ibd']] = analyze_snp(i)

    print(seed_snps.iloc[:100])
    seed_snps.to_csv(args.outfile, sep='\t')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed-snps', required=True)
    parser.add_argument('--typed-snps', required=True)
    parser.add_argument('--bfile', required=True)
    parser.add_argument('--seq-err', type=float, default=0.01)
    parser.add_argument('--gen-err', type=float, default=0.01)
    parser.add_argument('--p-thresh', type=float, default=0.05)
    parser.add_argument('--outfile', required=True)
    args = parser.parse_args()

    main(args)
