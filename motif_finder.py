import random
from sampling import sample
import matplotlib.pyplot as plt
from Bio import SeqIO
import matplotlib.pyplot as plt
import time

nucleotides = ["A", "T", "C", "G"]
dr = {"A":0, "T":1, "C":2, "G":3}

def calculate_prob(seqs, pos, seqth, ml):
    #except i
    q = {x: [1] * ml for x in nucleotides}
    p = {x: 1 for x in nucleotides}

    for i in range(len(seqs)):
        if i == seqth:
            continue
        for j in range(len(seqs[i])):
            if j < pos[i] or j > pos[i] + ml:
                c = seqs[i][j]
                p[c] = p[c] + 1

    for i in range(len(seqs)):
        if i == seqth:
            continue
        for j in range(ml):
            start_pos = pos[i]
            c = seqs[i][start_pos + j]
            q[c][j] = q[c][j] + 1

    for c in nucleotides:
        for j in range(ml):
            q[c][j] = q[c][j] / float(len(seqs) + len(nucleotides))

    total = sum(p.values())
    for c in nucleotides:
        p[c] = p[c] / float(total)

    return p, q

def gibbs_finder(dataset_name):
    fasta_file = "{}/sequences.fa".format(dataset_name)
    motif_len_file = "{}/motiflength.txt".format(dataset_name)

    f1 = open(motif_len_file, 'r')
    line = f1.readline().strip()
    ml = int(line)

    fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')

    seqs = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        seqs.append(sequence)

    K = len(seqs)
    N = len(seqs[0])

    pos = [random.randint(0, N - ml) for x in range(K)]

    MAX_ITER = 1000
    Aj = [0] * (N - ml)
    b = [0] * MAX_ITER
    time_start = time.time()
    for it in range(MAX_ITER):
        # We pick the sequence, well, in sequence starting from index 0
        for i in range(K):
            p, q = calculate_prob(seqs, pos, i, ml)
            qx = [1] * (N - ml)
            px = [1] * (N - ml)
            for j in range(N - ml):
                for k in range(ml):
                    c = seqs[i][j + k]
                    qx[j] = qx[j] * q[c][k]
                    px[j] = px[j] * p[c]

            Aj = [x / y for (x, y) in zip(qx, px)]
            norm_c = sum(Aj)
            Aj = list(map(lambda x: x / norm_c, Aj))
            #print(Aj)

            pos[i] = sample(range(N - ml), Aj)

        b[it] = max(Aj)
    time_end = time.time()
    # plt.plot(b)
    # plt.show()
    # plt.savefig('./show.png')

    motif_count = [[0 for i in range(4)] for j in range(ml)]
    for j in range(K):
        for i in range(ml):
            start_pos = pos[j]
            motif_count[i][dr[seqs[j][start_pos + i]]] += 1

    motif_pred = open("{}/predictedmotif.txt".format(dataset_name), "w")
    sites_pred = open("{}/predictedsites.txt".format(dataset_name), "w")
    time_pred = open("{}/runningtime.txt".format(dataset_name), "w")

    motif_pred.write(dataset_name + " " + str(ml) + "\n")
    time_pred.write(str(time_end - time_start))
    for i in range(ml):
        str_lst = [str(c) for c in motif_count[i]]
        line = ' '.join(str_lst)
        motif_pred.write(line + "\n")

    for i in range(K):
        sites_pred.write(str(pos[i]) + "\n")

for i in range(70):
    data_name = "data" + str(i)
    gibbs_finder(data_name)
    print('finish ' + data_name)

