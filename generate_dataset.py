import numpy as np
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from generate_pwm import generate_pwm

d = {0:"A", 1:"T", 2:"C", 3:"G"}
dr = {"A":0, "T":1, "C":2, "G":3}

def generate_background(sl):
    return np.random.randint(4, size=sl)

def generate_motif(idx, pwm):
    r = np.random.uniform(0,1)
    # print('generate_motif')
    # print(r)

    if r < pwm[idx][0]:
        #print(0, pwm[idx][0])
        return 0
    elif r < pwm[idx][0] + pwm[idx][1]:
        #print(1, pwm[idx][1])
        return 1
    elif r < pwm[idx][0]+ pwm[idx][1] + pwm[idx][2]:
        #print(2, pwm[idx][2])
        return 2
    #print(3, pwm[idx][3])
    return 3

def generate_sequence(sl, ml, pwm):
    #nonlocal motif_count
    global d
    site = np.random.randint(sl - ml + 1)
    seq_list = []
    seq = generate_background(sl)
    for x in seq:
        seq_list.append(d[x])
    for i in range(ml):
        #generate a number with probability
        m = generate_motif(i, pwm)
        seq_list[site + i] = d[m]
    #print(seq_list[site:site + ml])
    seq_str = ''.join(seq_list)
    return seq_str, site #string, int

def generate_dataset(num, ic, ml, sl, sc):
    cwd = os.getcwd()
    os.mkdir(cwd + '/data{}'.format(num))

    seq_file = open("data{}/sequences.fa".format(num), 'w+')
    sites_file = open("data{}/sites.txt".format(num), 'a+')
    motif_file = open("data{}/motif.txt".format(num), 'a+')
    motif_len_file = open("data{}/motiflength.txt".format(num), 'w+')

    pwm = generate_pwm(ic, ml) # ml x 4 matrix
    motif_count = [[0 for i in range(4)] for j in range(ml)]
    sequence_set = []
    for i in range(sc):
        seq_str, site = generate_sequence(sl, ml, pwm)
        #update motif_count
        for j in range(ml):
            motif_count[j][dr[seq_str[site + j]]] += 1
        sequence_set.append(seq_str)
        sites_file.write(str(site) + "\n")

    print(motif_count)
    records = (SeqRecord(Seq(seq, generic_dna), str(index)) for index,seq in enumerate(sequence_set) )
    SeqIO.write(records, seq_file, "fasta")

    motif_file.write("data" + str(num) + " " + str(ml) + "\n")
    for i in range(ml):
        str_lst = [str(c) for c in motif_count[i]]
        line = ' '.join(str_lst)
        motif_file.write(line + "\n")

    motif_len_file.write(str(ml))

#generate_dataset(72, 1, 6, 500, 5)
ics = [1, 1.5, 2]
mls = [6, 7, 8]
scs = [5, 10, 20]

default_ic = 2
default_ml = 8
default_sl = 500
default_sc = 10

for i in range(7):
    base = i * 10
    if i == 0:
        ic = default_ic
        ml = default_ml
        sc = default_sc
    elif i == 1:
        ic = ics[0]
        ml = default_ml
        sc = default_sc
    elif i == 2:
        ic = ics[1]
        ml = default_ml
        sc = default_sc
    elif i == 3:
        ic = default_ic
        ml = mls[0]
        sc = default_sc
    elif i == 4:
        ic = default_ic
        ml = mls[1]
        sc = default_sc
    elif i == 5:
        ic = default_ic
        ml = default_ml
        sc = scs[0]
    else:
        ic = default_ic
        ml = default_ml
        sc = scs[2]


    for j in range(10):
        generate_dataset(base + j, ic, ml, 500, sc)
