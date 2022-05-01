import numpy as np
import matplotlib.pyplot as plt

def evaluate(dataset_name):
    motif_len_file = open("{}/motiflength.txt".format(dataset_name), "r")

    motif_pred = open("{}/predictedmotif.txt".format(dataset_name), "r")
    sites_pred = open("{}/predictedsites.txt".format(dataset_name), "r")

    motif_truth = open("{}/motif.txt".format(dataset_name), "r")
    sites_truth = open("{}/sites.txt".format(dataset_name), "r")

    motif_len = int(motif_len_file.readline().strip())

    lines_pred = motif_pred.readlines()
    lines_truth = motif_truth.readlines()

    if len(lines_pred) != len(lines_truth) :
        print(len(lines_pred), len(lines_truth))
        exit(1)

    relative_entropy = 0
    alpha = 1e-3
    for i in range(motif_len + 1):
        if i == 0:
            continue
        count_pred = [int(c) for c in lines_pred[i].strip().split()]
        total = sum(count_pred)
        count_truth = [int(c) for c in lines_truth[i].strip().split()]

        for i in range(4):
            #(plnp - plnq)
            # 𝛼𝑈 | 𝑥 | +(1−𝛼)x
            p = float(count_truth[i]) / total
            q = float(count_pred[i]) / total
            p = alpha * 0.25 + (1 - alpha) * p
            q = alpha * 0.25 + (1 - alpha) * q

            relative_entropy += p * np.log(p) - p * np.log(q)

    relative_entropy = float(relative_entropy / motif_len)
    print('relative entropy ', relative_entropy)
    re_file = open("{}/relative_entropy.txt".format(dataset_name), "w")
    re_file.write(str(relative_entropy))

    sites_pred = sites_pred.readlines()
    sites_truth = sites_truth.readlines()

    overlapped_pos = 0
    for i in range(len(sites_pred)):
        d = abs(int(sites_pred[i].strip()) - int(sites_truth[i].strip()))
        if d < motif_len:
            overlapped_pos += motif_len - d

    print('overlapped_pos ', overlapped_pos)
    op_file = open("{}/overlap.txt".format(dataset_name), "w")
    op_file.write(str(overlapped_pos))

def plot_ic():
    ics = [1, 1.5, 2]
    re_means = []
    se_means = []
    runtimes = []
    se_means_2 = []
    overlaps = []
    se_means_3 = []
    for i in [1, 2, 0]:
        base = i * 10
        res = []
        runtime = []
        overlap = []
        for j in range(10):
            dataset_name = "data" + str(base + j)
            re_file = open("{}/relative_entropy.txt".format(dataset_name), "r")
            re = float(re_file.readlines()[0].strip())
            res.append(re)

            runtime_file = open("{}/runningtime.txt".format(dataset_name), "r")
            rt = float(runtime_file.readlines()[0].strip())
            runtime.append(rt)

            overlap_file = open("{}/overlap.txt".format(dataset_name), "r")
            ol = float(overlap_file.readlines()[0].strip())
            overlap.append(ol)

        data = np.array(res)
        re_mean = np.average(data)
        se_mean = np.std(data, ddof=1) / np.sqrt(np.size(data))
        re_means.append(re_mean)
        se_means.append(se_mean)

        data2 = np.array(runtime)
        runtimes.append(np.average(data2))
        se_means_2.append(np.std(data2, ddof=1) / np.sqrt(np.size(data2)))

        data3 = np.array([float(overlap[i]/10) for i in range(3)])
        overlaps.append(np.average(data3))
        print(overlaps)
        se_means_3.append(np.std(data3, ddof=1) / np.sqrt(np.size(data3)))

    plt.plot(ics, re_means, label="average relative entropy")
    plt.plot(ics, se_means, label="standard errors of mean")
    plt.legend()

    plt.savefig('./ic_entropy.png')
    plt.show()

    plt.plot(ics, runtimes, label="average running time")
    plt.plot(ics, se_means_2, label="standard errors of mean")
    plt.legend()

    plt.savefig('./ic_runtime.png')
    plt.show()

    plt.plot(ics, overlaps, label="average overlap")
    plt.plot(ics, se_means_3, label="standard errors of mean")
    plt.legend()

    plt.savefig('./ic_overlap.png')
    plt.show()

def plot_ml():
    mls = [6, 7, 8]
    re_means = []
    se_means = []
    runtimes = []
    se_means_2 = []
    overlaps = []
    se_means_3 = []
    for i in [3, 4, 0]:
        base = i * 10
        res = []
        runtime = []
        overlap = []
        for j in range(10):
            dataset_name = "data" + str(base + j)
            re_file = open("{}/relative_entropy.txt".format(dataset_name), "r")
            re = float(re_file.readlines()[0].strip())
            res.append(re)

            runtime_file = open("{}/runningtime.txt".format(dataset_name), "r")
            rt = float(runtime_file.readlines()[0].strip())
            runtime.append(rt)

            overlap_file = open("{}/overlap.txt".format(dataset_name), "r")
            ol = float(overlap_file.readlines()[0].strip())

            overlap.append(ol)

        data = np.array(res)
        re_mean = np.average(data)
        se_mean = np.std(data, ddof=1) / np.sqrt(np.size(data))
        re_means.append(re_mean)
        se_means.append(se_mean)

        data2 = np.array(runtime)
        runtimes.append(np.average(data2))
        se_means_2.append(np.std(data2, ddof=1) / np.sqrt(np.size(data2)))

        print(overlaps)
        data3 = np.array([float(overlap[i]/10) for i in range(3)])
        overlaps.append(np.average(data3))
        se_means_3.append(np.std(data3, ddof=1) / np.sqrt(np.size(data3)))

    plt.plot(mls, re_means, label="average relative entropy")
    plt.plot(mls, se_means, label="standard errors of mean")
    plt.legend()

    plt.savefig('./ml_entropy.png')
    plt.show()

    plt.plot(mls, runtimes, label="average running time")
    plt.plot(mls, se_means_2, label="standard errors of mean")
    plt.legend()

    plt.savefig('./ml_runtime.png')
    plt.show()

    plt.plot(mls, overlaps, label="average overlap")
    #plt.plot(mls, se_means_3, label="standard errors of mean")
    plt.legend()

    plt.savefig('./ml_overlap.png')
    plt.show()

def plot_sc():
    scs = [5, 10, 20]
    re_means = []
    se_means = []
    runtimes = []
    se_means_2 = []
    overlaps = []
    se_means_3 = []
    for i in [5, 0, 6]:
        base = i * 10
        res = []
        runtime = []
        overlap = []
        for j in range(10):
            dataset_name = "data" + str(base + j)
            re_file = open("{}/relative_entropy.txt".format(dataset_name), "r")
            re = float(re_file.readlines()[0].strip())
            res.append(re)

            runtime_file = open("{}/runningtime.txt".format(dataset_name), "r")
            rt = float(runtime_file.readlines()[0].strip())
            runtime.append(rt)

            overlap_file = open("{}/overlap.txt".format(dataset_name), "r")
            ol = float(overlap_file.readlines()[0].strip())
            overlap.append(ol)

        data = np.array(res)
        re_mean = np.average(data)
        se_mean = np.std(data, ddof=1) / np.sqrt(np.size(data))
        re_means.append(re_mean)
        se_means.append(se_mean)

        data2 = np.array(runtime)
        runtimes.append(np.average(data2))
        se_means_2.append(np.std(data2, ddof=1) / np.sqrt(np.size(data2)))

        data3 = np.array([float(overlap[i]/scs[i]) for i in range(3)])

        overlaps.append(np.average(data3))
        se_means_3.append(np.std(data3, ddof=1) / np.sqrt(np.size(data3)))

    plt.plot(scs, re_means, label="average relative entropy")
    plt.plot(scs, se_means, label="standard errors of mean")
    plt.legend()

    plt.savefig('./sc_entropy.png')
    plt.show()

    plt.plot(scs, runtimes, label="average running time")
    plt.plot(scs, se_means_2, label="standard errors of mean")
    plt.legend()

    plt.savefig('./sc_runtime.png')
    plt.show()

    plt.plot(scs, overlaps, label="average overlap")
    plt.plot(scs, se_means_3, label="standard errors of mean")
    plt.legend()

    plt.savefig('./sc_overlap.png')
    plt.show()

for i in range(70):
    data_name = "data" + str(i)
    evaluate(data_name)
    print('finish evaluate ' + data_name)

# plot_ic()
plot_ml()
# plot_sc()