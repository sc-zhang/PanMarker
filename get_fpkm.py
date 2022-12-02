#!/usr/bin/env python
import sys


def extract_fpkm(in_fpkm, gid, out_list):
    smp_db = {}
    with open(in_fpkm, 'r') as fin:
        with open(out_list, 'w') as fout:
            for line in fin:
                data = line.strip().split()
                if data[0] == 'gene_id':
                    for i in range(1, len(data)):
                        smp_db[i] = data[i]
        
                elif data[0] == gid:
                    for i in range(1, len(data)):
                        fout.write("%s\t%s\n"%(smp_db[i], data[i]))


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: python %s <in_fpkm> <gene_id> <out_list>"%sys.argv[0])
    else:
        in_fpkm, gid, out_list = sys.argv[1:]
        extract_fpkm(in_fpkm, gid, out_list)

