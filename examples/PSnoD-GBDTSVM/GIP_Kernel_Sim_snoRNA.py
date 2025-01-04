import os
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


# Assuming A is the binary association matrix with snoRNAs as rows and diseases as columns
# Define input and output directories
input_data_dir = "data"
output_dir = "results"
# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Set the input file path
ksd = os.path.join(input_data_dir, "known_snoRNA_disease.csv")  # Kernel similarity disease

#load ksd csv file
IP_matrix = np.loadtxt(ksd, delimiter=",")

rs = calculate_rm(IP_matrix)
KS_matrix = calculate_KS(IP_matrix, rs)

print("r_m:", rs)
print("KS matrix:")
print(KS_matrix)

# Save the KD matrix to a CSV file in the output directory
output_file = os.path.join(output_dir, "KS_matrix.csv")
np.savetxt(output_file, KS_matrix, delimiter=",")
print(f"KS matrix successfully saved to: {output_file}")

