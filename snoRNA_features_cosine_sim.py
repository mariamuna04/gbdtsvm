import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import cosine_similarity

df = pd.read_csv('MDRF_Data/snoRNA_features.csv')

# the first column contains IDs or, thats why we ignored it.
feature_matrix = df.iloc[:, 1:].values

# Step 2: Normalize the feature matrix
scaler = StandardScaler()
F_normalized = scaler.fit_transform(feature_matrix)

# Step 3: Calculate pairwise cosine similarity
similarity_matrix = cosine_similarity(F_normalized)

# Step 4: Save the similarity matrix to a CSV file
output_file = 'MDRF_output/snoRNA_similarity.csv'
pd.DataFrame(similarity_matrix).to_csv(output_file, index=False)