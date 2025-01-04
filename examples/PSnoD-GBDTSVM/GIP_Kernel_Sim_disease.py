import os
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


# Assuming IP_matrix is a numpy array where each row represents IP(d(i))
# Define input and output directories
input_data_dir = "data"
output_dir = "results"
# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Set the input file path
ksd = os.path.join(input_data_dir, "known_snoRNA_disease.csv")  # Kernel similarity disease

IP_matrix = np.loadtxt(ksd, delimiter=",")

rd = calculate_rd(IP_matrix)
KD_matrix = calculate_KD(IP_matrix, rd)

print("r_d:", rd)
print("KD matrix:")
print(KD_matrix)

# Save the KD matrix to a CSV file in the output directory
output_file = os.path.join(output_dir, "KD_matrix.csv")
np.savetxt(output_file, KD_matrix, delimiter=",")
print(f"KD matrix successfully saved to: {output_file}")


