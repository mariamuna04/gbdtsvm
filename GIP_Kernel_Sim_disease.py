import numpy as np

def calculate_KD(A, rd):
    # Number of diseases
    nd = A.shape[1]

    # Initializing KD matrix
    KD_matrix = np.zeros((nd, nd))

    # Calculating Kernel Disease values
    for i in range(nd):
        for j in range(nd):
            distance_squared = np.linalg.norm(A[:, i] - A[:, j]) ** 2
            KD_matrix[i, j] = np.exp(-rd * distance_squared)

    return KD_matrix


def calculate_rd(A):
    # Number of diseases
    nd = A.shape[1]

    # Calculating the average squared norm of interaction profiles
    sum_squared_norms = sum(np.linalg.norm(A[:, i]) ** 2 for i in range(nd))
    avg_squared_norm = sum_squared_norms / nd

    # r'_d is set to 1
    rd_prime = 1

    # Calculating rd-- parameter needed for calculation
    rd = rd_prime / avg_squared_norm

    return rd

# Example usage
# Assuming IP_matrix is a numpy array where each row represents IP(d(i))

ksd = "MDRF_Data/known_snoRNA_disease.csv" #ksd = kernel similarity disease

IP_matrix = np.loadtxt(ksd, delimiter=",")

rd = calculate_rd(IP_matrix)
KD_matrix = calculate_KD(IP_matrix, rd)

print("r_d:", rd)
print("KD matrix:")
print(KD_matrix)

#convert KD matrix to csv file
np.savetxt("MDRF_output/KD_matrix.csv", KD_matrix, delimiter=",")


