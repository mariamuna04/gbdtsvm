import random
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from numpy import interp
from sklearn.cluster import KMeans
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.preprocessing import OneHotEncoder
from sklearn import svm
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline

from sklearn.metrics import precision_recall_curve, roc_curve, auc, confusion_matrix, accuracy_score, classification_report
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

from sklearn.metrics import f1_score

sage_disease_sim = pd.read_csv(r'SAGESDA_Data\disease_sim_graph_filtered.csv')
sage_rna_sim = pd.read_csv(r'SAGESDA_Data\snoRNA_4mer_similarity.csv')
sage_dis_rna_relation = pd.read_csv(r'SAGESDA_Data\relationship_matrix_filtered.csv')

disease_name = sage_disease_sim.iloc[:, 0]
snoRNA_name = sage_rna_sim.iloc[:, 0]
SnoRNA_similarity = pd.read_csv('SAGESDA_output/IRS_matrix.csv', header=None)
known_association = sage_dis_rna_relation.iloc[:, 1:].T
disease_similarity = pd.read_csv('SAGESDA_output/IDS_matrix.csv', header=None)
