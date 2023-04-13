#!/usr/bin/env python
from scipy.stats import levene, ttest_ind
from outliers import smirnov_grubbs as grubbs
import numpy as np
import sys


def stat_var_vs_exp(in_var, in_exp, out_stat):
    print("Loading expressions")
    exp_db = {}
    with open(in_exp, 'r') as fin:
        for line in fin:
            data = line.strip().split()
            exp_db[data[0]] = float(data[1]) #np.log2(float(data[1])+1)
    
    print("Loading variants")
    with open(in_var, 'r') as fin:
        with open(out_stat, 'w') as fout:
            idx_db = {}
            for line in fin:
                data = line.strip().split()
                if line[0] == '#':
                    fout.write("%s\tLevene_pvalue\tTtest_pvalue\tValid_ref_count\tValid_alt_count\n"%line.strip())
                    smp_list = []
                    for i in range(4, len(data)):
                        idx_db[i] = data[i]
                        if idx_db[i] not in exp_db:
                            continue
                        smp_list.append(idx_db[i])
                else:
                    ref_list = []
                    alt_list = []
                    for i in range(4, len(data)):
                        if idx_db[i] not in exp_db:
                            continue
                        if data[i] == '0':
                            ref_list.append(exp_db[idx_db[i]])
                        elif data[i] == '1':
                            alt_list.append(exp_db[idx_db[i]])

                    if len(ref_list) == 1 or len(alt_list) == 1:
                        continue
                    # remove outliers with grubbs test
                    ref_list = grubbs.test(np.array(ref_list), alpha=0.05)
                    alt_list = grubbs.test(np.array(alt_list), alpha=0.05)
                    if len(ref_list) == 1 or len(alt_list) == 1:
                        continue
                    levene_val = levene(ref_list, alt_list).pvalue
                    if levene_val > 0.05:
                        equal_var = True
                    else:
                        equal_var = False
                    t = ttest_ind(ref_list, alt_list, equal_var=equal_var)
                    t_pval = t.pvalue
                    if t_pval <= 0.05:
                        fout.write("%s\t%.4f\t%.4f\t%d\t%d\n"%(line.strip(), levene_val, t_pval, len(ref_list), len(alt_list)))
    print("Finished")


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: python %s <in_var> <in_exp> <out_stat>"%sys.argv[0])
    else:
        in_var, in_exp, out_stat = sys.argv[1:]
        stat_var_vs_exp(in_var, in_exp, out_stat)

