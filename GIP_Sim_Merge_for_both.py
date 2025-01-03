import numpy as np
def integrate_disease_similarity(KD, SS):
    nd = KD.shape[0]
    SD = np.zeros((nd, nd))

    for i in range(nd):
        for j in range(nd):
            if SS[i, j] > 0:  # Assuming semantic similarity exists if SS[i, j] > 0
                SD[i, j] = (KD[i, j] + SS[i, j]) / 2
            else:
                SD[i, j] = KD[i, j]

    return SD


def integrate_snorna_similarity(KS, FS):  #KS = Kernel Similarity and FS= Functional Similarity of snoRNA
    ns = KS.shape[0]
    SM = np.zeros((ns, ns))

    for i in range(ns):
        for j in range(ns):
            if FS[i, j] > 0:  # Assuming functional similarity exists if FS[i, j] > 0
                SM[i, j] = (KS[i, j] + FS[i, j]) / 2
            else:
                SM[i, j] = KS[i, j]

    return SM

# 2 similarity metrices of disease
KD = np.loadtxt("MDRF_output/KD_matrix.csv", delimiter=",")
SS = np.loadtxt("MDRF_Data/disease_similarity.csv", delimiter=",")
# 2 similarity metrices of snoRNAs
KS = np.loadtxt("MDRF_output/KS_matrix.csv", delimiter=",")
FS = np.loadtxt("MDRF_output/snoRNA_similarity.csv", delimiter=",")


IDS = integrate_disease_similarity(KD, SS)
IRS = integrate_snorna_similarity(KS, FS)

print("Integrated disease similarity (SD):")
print(IDS)
np.savetxt("MDRF_output/IDS_matrix.csv", IDS, delimiter=",")

print("Integrated snoRNA similarity (SM):")
print(IRS)
np.savetxt("MDRF_output/IRS_matrix.csv", IRS, delimiter=",")
