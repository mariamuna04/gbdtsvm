import numpy as np

def calculate_KS(A, rm):
    # Number of snoRNAs
    nm = A.shape[0]

    # Initializing Kernel SnoRNA (KS) matrix
    KS_matrix = np.zeros((nm, nm))

    # Calculate KS values
    for i in range(nm):
        for j in range(nm):
            distance_squared = np.linalg.norm(A[i] - A[j]) ** 2
            KS_matrix[i, j] = np.exp(-rm * distance_squared)

    return KS_matrix

def calculate_rm(A):
    # Number of snoRNAs
    nm = A.shape[0]

    # Calculating the average squared norm of interaction profiles
    sum_squared_norms = sum(np.linalg.norm(A[i]) ** 2 for i in range(nm))
    avg_squared_norm = sum_squared_norms / nm

    # r'_m is set to 1
    rm_prime = 1

    # Calculating rs
    rm = rm_prime / avg_squared_norm

    return rm

# Example usage
# Assuming A is the binary association matrix with snoRNAs as rows and diseases as columns
ksd = "MDRF_Data/known_snoRNA_disease.csv"
#load ksd csv file
IP_matrix = np.loadtxt(ksd, delimiter=",")

rs = calculate_rm(IP_matrix)
KS_matrix = calculate_KS(IP_matrix, rs)

print("r_m:", rs)
print("KS matrix:")
print(KS_matrix)

np.savetxt("MDRF_output/KS_matrix.csv", KS_matrix, delimiter=",")
