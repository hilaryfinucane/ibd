from __future__ import print_function, division
from pysnptools.snpreader import Bed
import numpy as np
import pandas as pd
import itertools as it
import argparse

def main(args):
    print('reading seeed snps')
    seed_snps = pd.read_csv(args.seed_snps, header=None, names=['SNP'], index_col='SNP')
    seed_snps['ibs_length'] = 0
    seed_snps['ibd'] = 0

    print('reading typed snps')
    typed_snps = pd.read_csv(args.typed_snps, header=None, names=['SNP'])

    print('reading genotypes')
    data = Bed(args.bfile)
    X = data[:,:1000].read().val
    typed_snps_indices = np.sort(data.sid_to_index(typed_snps.SNP))
    typed_snps_bp = data.col_property[typed_snps_indices,2]

    print(len(seed_snps), 'snps in list')
    print(data.iid_count, data.sid_count, 'are dimensions of X')

    def analyze_snp(i):
        # find first typed snp after query snp
        snp_bp = data.col_property[i,2]
        typed_i = typed_snps_indices[np.where(typed_snps_bp > snp_bp)[0][0]]
        n1, n2 = np.where(X[:,i] == 1)[0]

        # extend list to left
        # extend list to right

        typed_il = typed_i - 1
        typed_ir = typed_i + 1
        ibd = 0

        cM = data.col_property[typed_ir, 1] - data.col_property[typed_il, 1]
        return cM, ibd

    for (i, snp) in it.izip(data.sid_to_index(seed_snps.index), seed_snps.index):
        print(i, snp)
        seed_snps.ix[snp, ['ibs_length', 'ibd']] = analyze_snp(i)
        if i > 100:
            break

    print(seed_snps.iloc[:10])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed-snps', required=True)
    parser.add_argument('--typed-snps', required=True)
    parser.add_argument('--bfile', required=True)
    parser.add_argument('--seq-err', type=float, default=0.01)
    parser.add_argument('--gen-err', type=float, default=0.01)
    parser.add_argument('--p-thresh', type=float, default=0.05)
    args = parser.parse_args()

    main(args)
