import os
import numpy as np
import pandas as pd

# Set the directory for the PSnoD datasets
sage_disease_sim = pd.read_csv(r'raw_data\disease_sim_graph_filtered.csv')
sage_rna_sim = pd.read_csv(r'raw_data\snoRNA_4mer_similarity.csv')
sage_dis_rna_relation = pd.read_csv(r'raw_data\relationship_matrix_filtered.csv')

disease_name = sage_disease_sim.iloc[:, 0]
snoRNA_name = sage_rna_sim.iloc[:, 0]
SnoRNA_similarity = pd.read_csv('SAGESDA_output/IRS_matrix.csv', header=None)
known_association = sage_dis_rna_relation.iloc[:, 1:].T
disease_similarity = pd.read_csv('SAGESDA_output/IDS_matrix.csv', header=None)
