import os
import numpy as np
import pandas as pd

# Set the directory for the PSnoD datasets
disease_sim = pd.read_csv(r'raw_data\disease_sim_graph_filtered.csv')
rna_sim = pd.read_csv(r'raw_data\snoRNA_4mer_similarity.csv')
disease_rna_relation = pd.read_csv(r'raw_data\relationship_matrix_filtered.csv')

# output formatted data directory

output_dir = 'data'

disease_name = disease_sim.iloc[:, 0]
snoRNA_name = rna_sim.iloc[:, 0]
known_association = disease_rna_relation.iloc[:, 1:].T
disease_similarity = disease_sim.iloc[:, 1:]
SnoRNA_similarity = rna_sim.iloc[:, 1:]


# save to csv
disease_similarity_csv_path = os.path.join(output_dir, "disease_similarity.csv")
disease_similarity.to_csv(disease_similarity_csv_path, index=False, header=None)
