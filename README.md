# TCR Immune Repertoire Classification

The method utilizes [GIANA](https://github.com/s175573/GIANA) for ultrafast clustering of similar CDR3 sequences and outputs training and validation data that can be used to train [DeepCAT](https://github.com/s175573/DeepCAT), a computational method that identifies cancer associated beta chain TCRs.

This method is written in Python3 and requires the [pandas](https://pandas.pydata.org/) library to be installed.

## Standard Pipeline

Concatenate TCR immune repertoire data from cancer and non-cancer patient samples into one file
Run GIANA on the concatenated file to produce a cluster file
Input the cluster file and CSV file of labeled samples

## Input data format

**Cluster File:** The first and second columns are required to be the CDR3 amino acid sequence and cluster ID, respectively. The last column is also required to be the sample ID or filename. OTher columns in the input data may contain any information.

**CSV of Labeled Samples:** Sample ID and patient donor type is necessary to differentiate clusters. The labeled csv file should contain columns labeled “Sample” and “Benign/Malignant” where “Sample” holds the sample ID or filename corresponding to that of the cluster file. “Benign/Malignant” must have a value of “Benign” or “Malignant” to identify the patient donor types. 

## Usage
The following code performs TCR immune repertoire classification:
>```python getTrainingData.py GIANA_output.txt sampleLabels.csv cluster_size(int; optional; default=50) cluster_purity(float; optional; default=0.8)```

The method requires two input files: clustered output from GIANA and a labeled CSV file of the sample ID and donor type.

Users can optionally customize cluster size and purity thresholds by adjusting the third and fourth arguments, respectively. Leaving these arguments blank results in the method using a default cluster size of <= 50 and a cluster purity of 80%.

The method creates a new directory “TrainingData” and generates four output files for cancer and non-cancer training and validation sets. Non-cancer training data can be used as true negative samples and cancer training data can be used as true positive samples. Sequences are randomly shuffled and split such that 80% of sequences for each length are reserved for training and 20% are withheld for validation.
