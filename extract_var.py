#!/usr/bin/env python
import sys


def load_aln(aln_db, in_aln):
    with open(in_aln, 'r') as fin:
        id = ""
        seq = ""
        for line in fin:
            if line[0] == '>':
                if seq != "":
                    aln_db[id] = seq
                id = line.strip().split()[0][1:].replace('_R_', '')
                seq = ""
            else:
                seq += line.strip().upper()
        aln_db[id] = seq


def get_consensus_seq(aln_db):
    consensus_seq = ""
    smp_list = sorted(aln_db.keys())
    for _ in range(len(aln_db[smp_list[0]])):
        cnt_db = {}
        for smp in smp_list:
            base = aln_db[smp][_]
            if base not in cnt_db:
                cnt_db[base] = 0
            cnt_db[base] += 1
        max_base = ""
        max_cnt = 0
        for base in cnt_db:
            if cnt_db[base] > max_cnt:
                max_cnt = cnt_db[base]
                max_base = base
        if max_cnt == len(smp_list):
            consensus_seq += max_base.upper()
        else:
            consensus_seq += max_base.lower()
    return consensus_seq


def classify_info(classified_info, full_info):
    last_info = []
    for _ in range(len(full_info)):
        cur_info = full_info[_]
        if cur_info[1] == '-':
            type = '<INS>'
        elif cur_info[2] == '-':
            type = '<DEL>'
        else:
            type = '<SNP>'
        if last_info == []:
            last_info = [type, cur_info[0], cur_info]
        else:
            if last_info[0] == '<SNP>':
                classified_info.append([last_info[2][0], last_info[0], last_info[2][1], last_info[2][2], last_info[2][3]])
                last_info = [type, cur_info[0], cur_info]
            else:
                if type != last_info[0]:
                    classified_info.append([last_info[2][0], last_info[0], last_info[2][1], last_info[2][2], last_info[2][3]])
                    last_info = [type, cur_info[0], cur_info]
                else:
                    if cur_info[0] != last_info[1]+1:
                        classified_info.append([last_info[2][0], last_info[0], last_info[2][1], last_info[2][2], last_info[2][3]])
                        last_info = [type, cur_info[0], cur_info]
                    else:
                        if last_info[2][3] != cur_info[3]:
                            classified_info.append([last_info[2][0], last_info[0], last_info[2][1], last_info[2][2], last_info[2][3]])
                            last_info = [type, cur_info[0], cur_info]
                        else:
                            last_info[1] = cur_info[0]
                            if cur_info[1] != '-':
                                last_info[2][1] += cur_info[1]
                            if cur_info[2] != '-':
                                last_info[2][2] += cur_info[2]
    if last_info:
        classified_info.append([last_info[2][0], last_info[0], last_info[2][1], last_info[2][2], last_info[2][3]])


def extract_var(in_aln, out_stat):
    print("Loading alignment")
    aln_db = {}
    load_aln(aln_db, in_aln)

    print("Getting variants")
    with open(out_stat, 'w') as fout:
        smps = sorted(aln_db.keys())
        ref_seq = get_consensus_seq(aln_db)
        seq_len = len(ref_seq)
        #fout.write("#>Consensus\n#%s\n"%ref_seq)
        fout.write("#POS\tTYPE\tREF\tALT\t%s\n"%('\t'.join(smps)))

        # Indetify each pos
        full_info = []
        for i in range(seq_len):
            info = []
            ref = ref_seq[i].upper()
            alt = ""
            alt_cnt = {}
            for smp in smps:
                base = aln_db[smp][i]
                if base != ref:
                    type = 1
                    if base not in alt_cnt:
                        alt_cnt[base] = 0
                    alt_cnt[base] += 1
                else:
                    type = 0
                info.append(type)
            max_cnt = 0
            for base in alt_cnt:
                if alt_cnt[base] > max_cnt:
                    alt = base
                    max_cnt = alt_cnt[base]
            if alt == "" or max_cnt <= 1:
                continue
            full_info.append([i, ref, alt, '\t'.join(map(str, info))])
        classified_info = []
        classify_info(classified_info, full_info)
        for info in classified_info:
            fout.write("%s\n"%('\t'.join(map(str, info))))
    
    print("Finished")
    

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python %s <in_aln> <out_stat>"%sys.argv[0])
    else:
        in_aln, out_stat = sys.argv[1:]
        extract_var(in_aln, out_stat)

