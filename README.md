# AutoCAT: Automated Cancer-Associated TCR discovery

AutoCAT is a computational method to predict tumor associated TCRs from targeted TCR-seq data. The method utilizes [GIANA](https://github.com/s175573/GIANA) to quickly cluster similar CDR3 sequences. These clusters are then filtered to output training and validation data that can be used to train [DeepCAT](https://github.com/s175573/DeepCAT), a deep learning algorithm that identifies cancer associated beta chain TCRs. AutoCAT is written in Python3 and requires the [pandas](https://pandas.pydata.org/) library to be installed.

## Standard Pipeline

1. Run GIANA on a concatenated file of all TCR-seq data from cancer and non-cancer patient samples to produce a cluster file
2. Input the cluster file and CSV file of labeled samples to generate training and validation data files

## Input data format

**Cluster File:** The first and second columns are required to be the CDR3 amino acid sequence and cluster ID, respectively. The last column is also required to be the sample ID or filename. Other columns in the input data may contain any information.

**CSV of Labeled Samples:** Sample ID and patient donor type is necessary to differentiate clusters. The labeled CSV file should contain columns labeled “Sample” and “Benign/Malignant” where “Sample” holds the sample ID or filename corresponding to that of the cluster file. “Benign/Malignant” must have a value of “Benign” or “Malignant” to identify the patient donor types. 

## Usage
The following code performs TCR immune repertoire classification:

```python getTrainingData.py -f GIANA_output.txt -l sampleLabels.csv -s cluster_size -p cluster_purity```

| Command| Description |
| ----------- | ----------- |
|-h, --help| show this help message and exit |
| -f FILE, --file=FILE| Input cluster file of clustered CDR3 sequences |
| -l FILE, --labels=FILE| Labeled CSV with sample IDs and donor type |
| -s SIZE, --size=SIZE| Cluster size. Size must be an integer. If omitted, default cluster size is 50. |
| -p PURITY, --purity=PURITY| Cluster purity. Purity must be a float ≥0 and ≤1. If omitted, default cluster purity is 0.80.|
|-o OUTPUT_DIR, --outputDIR=OUTPUT_DIR| Output directory for training and validation files. If not given, a directory named 'TrainingData' will be created to hold the training and validation data files. |

The method requires two input files: the cluster file output from GIANA and a labeled CSV file of the sample ID and donor type. Users can optionally customize cluster size and purity. Leaving these arguments blank results in the method using a default cluster size of ≤50 and a cluster purity of 80%. This method creates a new directory and generates four output files for cancer and non-cancer training and validation sets. Users can optionally name this directory or use the default directory name “TrainingData.” Non-cancer training data can be used as true negative samples and cancer training data can be used as true positive samples. Sequences are randomly shuffled and split such that 80% of sequences for each length are reserved for training and 20% are withheld for validation.


