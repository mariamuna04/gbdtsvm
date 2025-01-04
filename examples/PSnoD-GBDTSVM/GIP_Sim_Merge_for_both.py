import os
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

# Define input and output directories
input_data_dir = "data"
output_dir = "results"
# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Set the input file path
KD_path = os.path.join(output_dir, "KD_matrix.csv")  
SS_path = os.path.join(input_data_dir, "disease_similarity.csv")  
KS_path = os.path.join(output_dir, "KS_matrix.csv")  
FS_path = os.path.join(output_dir, "snoRNA_similarity.csv")  


# 2 similarity metrices of disease
KD = np.loadtxt(KD_path, delimiter=",")
SS = np.loadtxt(SS_path, delimiter=",")
# 2 similarity metrices of snoRNAs
KS = np.loadtxt(KS_path, delimiter=",")
FS = np.loadtxt(FS_path, delimiter=",")


IDS = integrate_disease_similarity(KD, SS)
IRS = integrate_snorna_similarity(KS, FS)

IDS_file = os.path.join(output_dir, "IDS_matrix.csv")
IRS_file = os.path.join(output_dir, "IRS_matrix.csv")

print("Integrated disease similarity (SD):")
print(IDS)
np.savetxt(IDS_file, IDS, delimiter=",")

print("Integrated snoRNA similarity (SM):")
print(IRS)
np.savetxt(IRS_file, IRS, delimiter=",")
print(f"IDS and IRD matrix successfully saved to: {output_dir}")
