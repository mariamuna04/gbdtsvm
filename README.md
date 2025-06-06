# GBDTSVM Framework for snoRNA-Disease Association Prediction
## GBDTSVM: Combined  Support Vector Machine and Gradient Boosting Decision Tree Framework for efficient snoRNA-disease association prediction

## Key Points

- **Prediction Approach**: Utilizes Gradient Boosting Decision Tree (GBDT) followed by a Support Vector Machine (SVM) to predict snoRNA-disease associations.
- **Data Integration**: Leverages snoRNA functional similarity and disease semantic similarity with Gaussian kernel similarity measures.
- **Performance**: Outperforms existing state-of-the-art models in terms of prediction accuracy.
- **Validation**: Predictions validated using case studies and existing literature.


**Paper link**: [click here](https://doi.org/10.1016/j.compbiomed.2025.110219)

---

## Graphical Abstract

![Graphical Abstract](https://github.com/mariamuna04/gbdtsvm/blob/main/Figures/Graphical_Abstract.png?raw=true)

---

## Table of Contents

1. [Installation](#installation)
2. [Train and Evaluate the GBDTSVM](#train-and-evaluate-the-gbdtsvm)
   - [Dataset Descriptions](#dataset-descriptions)
   - [Code Descriptions](#code-descriptions)
3. [Example](#example)

---

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/mariamuna04/gbdtsvm.git
   cd gbdtsvm
   ```
2. Install required dependencies:
   ```bash
   pip install -r requirements.txt
   ```

---

## Train and Evaluate the GBDTSVM

### Dataset Descriptions

GBDTSVM expects 5 information about the datasets in csv format consists of the following files (**csv file name provided below according to our datasets for MDRF datasets mentioned in our paper**):

- **Diseases Name**: `disease_name.csv`
- **Diseases-to-Diseases Similarity**: `disease_similarity.csv`
- **snoRNA Name**: `snoRNA_name.csv`
- **snoRNA Features**: `snoRNA_features.csv`
- **Known Associations between snoRNA and Diseases**: `known_snoRNA_disease.csv` (**where 1 indicates there is an association between the snoRNA and disease, 0 indicates the association is unknown between that snoRNA and disease.**)

### Code Descriptions
We have 5 code files for different tasks where the prediction score is generated using `GBDTSVM.py` file and other 4 code files are used to prepare the final datasets that would be used as inputs in the `GBDTSVM.py` file.

#### `snoRNA_features_cosine_sim.py`
- **Input**: `snoRNA_features.csv`
- **Output**: `snoRNA_similarity.csv`
- **Description**: Calculates the cosine similarity between snoRNA features. **If your datasets already have similarity score between the snoRNAs, then you do not need to use this code.**

#### `GIP_Kernel_Sim_disease.py`
- **Input**: `known_snoRNA_disease.csv`
- **Output**: `KD_matrix.csv` (Gaussian kernel profile disease similarity matrix)
- **Description**: Measures disease similarities using the Gaussian kernel strategy.

#### `GIP_Kernel_Sim_snoRNA.py`
- **Input**: `known_snoRNA_disease.csv`
- **Output**: `KS_matrix.csv` (Gaussian kernel profile snoRNA similarity matrix)
- **Description**: Measures snoRNA similarities using the Gaussian kernel strategy.

#### `GIP_Sim_Merge_for_both.py`
- **Inputs**:
  - `KD_matrix.csv` (Disease Gaussian kernel similarity matrix)
  - `disease_similarity.csv` (Disease semantic similarity matrix)
  - `KS_matrix.csv` (snoRNA Gaussian kernel similarity matrix)
  - `snoRNA_similarity.csv` (snoRNA functional similarity matrix)
- **Outputs**:
  - `IDS_matrix.csv` (Integrated disease similarity matrix)
  - `IRS_matrix.csv` (Integrated snoRNA similarity matrix)
- **Working**:
  - Merges disease semantic similarity (`SS`) and Gaussian kernel disease similarity (`KD`) to produce `IDS_matrix.csv`.
  - Merges snoRNA functional similarity (`FS`) and Gaussian kernel snoRNA similarity (`KS`) to produce `IRS_matrix.csv`.

#### `GBDTSVM.py`

- **Inputs**:
  - `disease_name.csv`
  - `snoRNA_name.csv`
  - `known_snoRNA_disease.csv`
  - `IDS_matrix.csv`
  - `IRS_matrix.csv`
- **Working**:
  - **Data Preparation**:
    1. Converts input CSV files into arrays.
    2. Separates known and unknown associations.
    3. Concatenates features for unknown associations (e.g., combining disease semantic similarity and snoRNA functional similarity).
    4. Clusters unknown pairs into multiple clusters using K-Means clustering method.
    5. Samples unknown pairs proportionally from each cluster.
    6. Combines sampled unknown data with known data into a single dataset.
  - **Training**:
    1. Splits data for 5-fold cross-validation.
    2. Fits the dataset and labels into the GBDT model.
    3. Applies one-hot encoding to the outputs of the GBDT decision tree leaf nodes.
    4. Fits the transformed data into an SVM with an RBF kernel.
    5. Trains all data and tests all unknown pairs that were not used in the training or the testing set.
    6. Saves prediction scores (≥ 0.75). 
    7. Sorts prediction scores by disease.
  - **Validation**:
    - Top associations for 12 diseases are verified using PubMed, RNADisease, etc. database.

- **Outputs**:
  - `unknown_true_svm.txt` or `unknown_true_svm.csv`: List of prediction scores for unknown pairs.
  - `unknown_true_75_svm.csv`: List of prediction scores for unknown pairs having prediction score greater than 0.75.
  

---

## Example
GBDTSVM expects 5 inputs as csv for association predictions. However, what if a dataset has different structure? For such datasets, we show this examples which takes another datasets and predicts for association using our proposed method.
To train the GBDTSVM model:

For example, we perform GBDTSVM on datasets curated from MNDR V3.1 in the work of Sun et al. ([https://github.com/linDing-groups/PSnoD](https://github.com/linDing-groups/PSnoD)). 
This dataset has 3 files that can be accessed in `examples/PSnoD-GBDTSVM/raw_data`:
1. SnoRNA Similarity Matrix which represents the similarity between snoRNAs (`snoRNA_4mer_similarity.csv`).
2. Disease Similarity Matrix which represents the similarity between diseases (`disease_sim_graph_filtered.csv`).
3. Association Data which contains association information based on the known snoRNA-disease associations (`relationship_matrix_filtered.csv`).

As GBDTSVM.py expects 5 inputs, we first format these three files according to the need of GBDTSVM Framework using `data_format.py`. Please run this python file that will save
proper formatted files in the `examples/PSnoD-GBDTSVM/data` directory. Then, we can run the `GIP_Kernel_Sim_disease.py`, `GIP_Kernel_Sim_snoRNA.py`, `GIP_Sim_Merge_for_both.py` and `GBDTSVM.py` file sequentially after setting the respective input and output directory. We will get the prediction score in the `examples/PSnoD-GBDTSVM/results` folder.

Please note that, `PSnoD` dataset has directory snoRNA similarity, hence we do not need to run the `snoRNA_features_cosine_sim.py` code.


Feel free to reach out with any questions or issues! [Ummay Maria Muna](mailto:umuna201429@bscse.uiu.ac.bd), [Fahim Hafiz](mailto:fahimhafiz@cse.uiu.ac.bd)

## If you use our work, kindly cite our work!
@article{MUNA2025110219,
title = {GBDTSVM: Combined Support Vector Machine and Gradient Boosting Decision Tree Framework for efficient snoRNA-disease association prediction},
journal = {Computers in Biology and Medicine},
volume = {192},
pages = {110219},
year = {2025},
issn = {0010-4825},
doi = {https://doi.org/10.1016/j.compbiomed.2025.110219},
url = {https://www.sciencedirect.com/science/article/pii/S0010482525005700},
author = {Ummay Maria Muna and Fahim Hafiz and Shanta Biswas and Riasat Azim}
}