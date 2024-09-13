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


disease_name = pd.read_csv('Data/disease_name.csv')
snoRNA_name = pd.read_csv('Data/snoRNA_name.csv')
SnoRNA_similarity = pd.read_csv('Data/IRS_matrix.csv', header=None)
known_association = pd.read_csv('Data/known_snoRNA_disease.csv', header=None)
disease_similarity = pd.read_csv('Data/IDS_matrix.csv', header=None)

disease_semantic_similarity = np.zeros(disease_similarity.shape) ## 111 x 111
snoRNA_functional_similarity = np.zeros(SnoRNA_similarity.shape) ## 335 x 335
adjacency_matrix = np.zeros(known_association.shape) ## 335 x 111


# csv to array disease_semantic_similarity
for i in range(len(disease_name)):
    for j in range(len(disease_name)):
        disease_semantic_similarity[i, j] = disease_similarity.iloc[i, j]

# print(known_association.shape[0]) # snoRNAs --> 335
# print(known_association.shape[1]) # diseases --> 111
# csv to array adjacency_matrix
for i in range(known_association.shape[0]):
    for j in range(known_association.shape[1]):
        adjacency_matrix[i, j] = known_association.iloc[i, j]

# csv to array snoRNA_functional_similarity
for i in range(len(snoRNA_name)):
    for j in range(len(snoRNA_name)):
        snoRNA_functional_similarity[i, j] = SnoRNA_similarity.iloc[i, j]


unknown = []
known = []
for x in range(known_association.shape[0]):
    for y in range(known_association.shape[1]):
        if adjacency_matrix[x, y] == 0:
            unknown.append((x, y))
        else:
            known.append((x, y))

major = []
for z in range(len(unknown)):
    a = disease_semantic_similarity[unknown[z][1], :].tolist()
    b = snoRNA_functional_similarity[unknown[z][0], :].tolist()
    q = a + b
    major.append(q)


kmeans = KMeans(n_clusters=23, random_state=0).fit(major)
center = kmeans.cluster_centers_
labels = kmeans.labels_ # label is given to all the datapoints but within 1 to 20


# we have seperated x and y of each centre pair.
center_x = []
center_y = []
for j in range(len(center)):
    center_x.append(center[j][0])
    center_y.append(center[j][1])

disease_rna_tup = [[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []]
for i in range(len(labels)):
    if labels[i] == 0:
        disease_rna_tup[0].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 1:
        disease_rna_tup[1].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 2:
        disease_rna_tup[2].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 3:
        disease_rna_tup[3].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 4:
        disease_rna_tup[4].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 5:
        disease_rna_tup[5].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 6:
        disease_rna_tup[6].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 7:
        disease_rna_tup[7].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 8:
        disease_rna_tup[8].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 9:
        disease_rna_tup[9].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 10:
        disease_rna_tup[10].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 11:
        disease_rna_tup[11].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 12:
        disease_rna_tup[12].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 13:
        disease_rna_tup[13].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 14:
        disease_rna_tup[14].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 15:
        disease_rna_tup[15].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 16:
        disease_rna_tup[16].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 17:
        disease_rna_tup[17].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 18:
        disease_rna_tup[18].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 19:
        disease_rna_tup[19].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 20:
        disease_rna_tup[20].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 21:
        disease_rna_tup[21].append((unknown[i][0], unknown[i][1]))
    elif labels[i] == 22:
        disease_rna_tup[22].append((unknown[i][0], unknown[i][1]))



sampled_disease_rna_tup = [[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []]
for i in range(len(disease_rna_tup)):
    sampled_disease_rna_tup[i] = random.sample(disease_rna_tup[i], int((len(disease_rna_tup[i])/len(labels)) * len(known)))

dataset = []
for rna in range(known_association.shape[0]):
    for disease in range(known_association.shape[1]):
        for i in range(len(sampled_disease_rna_tup)):
            if (rna, disease) in sampled_disease_rna_tup[i]:
                dataset.append((rna, disease))

for rna in range(known_association.shape[0]):
    for disease in range(known_association.shape[1]):
        if (rna, disease) in known:
            dataset.append((rna, disease))

length = len(dataset)

selected_x = []
selected_y = []
#now I am just taking only the similarities of disease and rna of sampled data.
for data in dataset:
    a = disease_semantic_similarity[data[1], :].tolist()
    b = snoRNA_functional_similarity[data[0], :].tolist()
    q = a + b
    selected_x.append(q)

    if (data[0], data[1]) in known:
        selected_y.append(1)
    else:
        selected_y.append(0)

selected_data_np = np.array(selected_x)
selected_label_np = np.array(selected_y)

## Data prepation ends here before training the model--------------------------------

#### Gradient Boosting Classifier--------------------------------

GBDT=GradientBoostingClassifier(n_estimators = 12,max_depth=5,min_samples_leaf=13)
GBDT.fit(selected_data_np, selected_label_np)
OHE = OneHotEncoder()
OHE.fit(GBDT.apply(selected_data_np)[:, :, 0])


tprs = []
aucs = []
mean_fpr = np.linspace(0,1,100)

# Assuming xs is your input data-snoRNA-disease and ys is your target variable
# X_train, X_test, y_train, y_test = train_test_split(selected_data_np,selected_label_np, test_size=0.2, random_state=42)

SVM = svm.SVC(kernel='rbf', probability=True)

# Create a pipeline that scales the data-snoRNA-disease and then applies SVM
pipeline = Pipeline([('scaler', StandardScaler(with_mean=False)), ('SVM', SVM)])
# Define the parameter grid for C and gamma
param_grid = {'SVM__C': [0.1, 1, 10, 100], 'SVM__gamma': [1, 0.1, 0.01, 0.001]}

# Use GridSearchCV to find the optimal parameters
grid = GridSearchCV(pipeline, param_grid, cv=5)

# # Train the classifier
# grid.fit(OHE.transform(GBDT.apply(X_train)[:, :, 0]), y_train)
#
# predicted_probs = grid.predict_proba(OHE.transform(GBDT.apply(X_test)[:, :, 0]))[:, 1]
# # roc_curve calculates the true positive rate and false positive rate for different thresholds
# fpr, tpr, thresholds = roc_curve(y_test, predicted_probs)
# roc_auc = auc(fpr, tpr)
# print("Test Set ROC AUC:", roc_auc)


# using 5 fold cross validation:

stratified_k_fold = StratifiedKFold(n_splits=5)
# Initialize lists to store predicted and actual values
tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)
precisions = []
recalls = []
roc_auc_scores = []
scores = []
mean_recall = np.linspace(0, 1, 100)
y_true = []
y_pred = []

for train_index, test_index in stratified_k_fold.split(selected_data_np, selected_label_np):
    # Split the data into training and test sets
    X_train, X_test = selected_data_np[train_index], selected_data_np[test_index]
    y_train, y_test = selected_label_np[train_index], selected_label_np[test_index]

    grid.fit(OHE.transform(GBDT.apply(X_train)[:, :, 0]), y_train)   # Train the classifier

    predicted_probs = grid.predict_proba(OHE.transform(GBDT.apply(X_test)[:, :, 0]))[:, 1]  # Make predictions on the test set using the best parameters

    roc_auc = roc_auc_score(y_test, predicted_probs)  # Calculate the ROC AUC score and append it to the list
    print('ROC AUC:', roc_auc)
    roc_auc_scores.append(roc_auc)


    score = grid.score(OHE.transform(GBDT.apply(X_test)[:, :, 0]), y_test)  # Calculate the score and append it to the list
    print('Score:', score)
    scores.append(score)


    fpr, tpr, thresholds = roc_curve(y_test, predicted_probs)   # Compute ROC curve and area under the curve
    tprs.append(interp(mean_fpr, fpr, tpr))
    tprs[-1][0] = 0.0
    roc_auc = auc(fpr, tpr)
    aucs.append(roc_auc)


    precision, recall, _ = precision_recall_curve(y_test, predicted_probs)  # Compute Precision-Recall curve
    precisions.append(interp(mean_recall, recall[::-1], precision[::-1]))

    y_true.extend(y_test)
    predicted = grid.predict(OHE.transform(GBDT.apply(X_test)[:, :, 0]))
    y_pred.extend(predicted)


f1 = f1_score(y_test, predicted)
print("F1 Score:", f1)


# Plotting ROC curve
plt.figure(figsize=(8, 6))
mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
plt.plot(mean_fpr, mean_tpr, color='b', label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc), lw=2, alpha=.8)
std_tpr = np.std(tprs, axis=0)
tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2, label=r'$\pm$ 1 std. dev.')
plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Random')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC)')
plt.legend()
plt.show()


# Plotting PRC curve
plt.figure(figsize=(8, 6))
mean_precision = np.mean(precisions, axis=0)
mean_precision[0] = 1.0
mean_auc_pr = auc(mean_recall, mean_precision)
plt.plot(mean_recall, mean_precision, color='b', label='Mean Precision-Recall curve (AUPRC = {:.2f})'.format(mean_auc_pr), lw=2, alpha=.8)
std_precision = np.std(precisions, axis=0)
precisions_upper = np.minimum(mean_precision + std_precision, 1)
precisions_lower = np.maximum(mean_precision - std_precision, 0)
plt.fill_between(mean_recall, precisions_lower, precisions_upper, color='grey', alpha=.2, label=r'$\pm$ 1 std. dev.')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve')
plt.legend()
plt.show()

#
# # Now predicting associations for all unknown pairs
# unknown_pair = []
#
# for data in unknown:
#     a1 = disease_semantic_similarity[data[1], :].tolist()
#     b1 = snoRNA_functional_similarity[data[0], :].tolist()
#     q1 = a1 + b1
#     unknown_pair.append(q1)

# #so far we have used certain number of samples from all the clustered unknown pairs, but now will predict for all the unknown pairs whether they have association or not
# # here x1 contains the concatanated semantic similarity of diseases and functional similarity of snoRNAs of all unknown pairs
#
# predicted_probabilities_all_unknowns = grid.predict_proba(OHE.transform(GBDT.apply(unknown_pair)[:, :, 0]))
# unknown_pair_true_class = predicted_probabilities_all_unknowns[:, 1].tolist()
#
# # here I have sorted the true class of all unknown pairs in descending order.
# unknown_pair_true_class_np = np.array(unknown_pair_true_class)
# sorted_probabilities = -np.sort(-unknown_pair_true_class_np, axis=None, kind='heapsort')
#
#
# unknown_pair_true_class_matrix = np.matrix(unknown_pair_true_class)
# sorted_class_index = np.argsort(-unknown_pair_true_class_matrix).tolist()
# # print(len(sorted_class_index))
# sorted_class_index = sorted_class_index[0]
# # print(len(sorted_class_index))
#
# file_unknown_true = open("mdrf_unknown_true_svm.txt", 'w')
# file_unknown_true.writelines(['disease', '\t', 'SnoRNA', '\t', 'Score', '\n'])
# for i in range(len(sorted_class_index)):
#     file_unknown_true.writelines([str(unknown[sorted_class_index[i]][1]), '\t', str(unknown[sorted_class_index[i]][0]), '\t', str(unknown_pair_true_class[sorted_class_index[i]]), '\n'])
# file_unknown_true.close()
#
#
# # convert unknown_true.txt into csv
# df = pd.read_csv('mdrf_unknown_true_svm.txt', delimiter='\t')
# df.to_csv('mdrf_unknown_true_svm.csv', index=False)
#
#
# # created another csv file where only the row having greater or equal 0.75 value will be stored
# df = pd.read_csv('MDRF-Unknown-Results/mdrf_unknown_true_svm.csv')
# df['Score'] = df['Score'].astype(float)
# df = df[df['Score'] >= 0.75]
# df.to_csv('mdrf_unknown_true_75_svm.csv', index=False)