# AutoCAT: Automated Cancer-Associated TCR discovery

AutoCAT is a computational method to predict tumor-associated TCRs from targeted TCR-seq data. The method utilizes [GIANA](https://github.com/s175573/GIANA) to quickly cluster similar CDR3 sequences. These clusters are then filtered to output training and validation data that can be used to train [DeepCAT](https://github.com/s175573/DeepCAT), a deep learning algorithm that identifies cancer associated beta chain TCRs. AutoCAT acts as a bridge to connect GIANA and DeepCAT. 

AutoCAT is written in Python3 and requires the following dependencies to be installed:
* [biopython](https://biopython.org/)
* [faiss](https://github.com/facebookresearch/faiss)
* [scikit-learn](https://scikit-learn.org/stable/)
* [matplotlib](https://matplotlib.org/)
* [pandas](https://pandas.pydata.org/)

After installing these dependencies, please download the AutoCAT repository, which includes sample data as well as GIANA (version 4) and the associated TRBV allele data (Imgt_Human_TRBV.fasta).

## Standard Pipeline

1. To process input files, AutoCAT requires a directory with data files divided into two subdirectories, “Cancer” and “Control.” Files are then concatenated together and output as “trainingData.txt.”
2. GIANA is run on the concatenated file to produce a cluster file (“trainingData--RotationEncodingBL62.txt”). Clusters are then filtered via cluster size and purity cutoffs which can be optionally customized by users. Leaving these arguments blank results in the method using a default cluster size of ≤50 and a cluster purity of 80%.
3. Training and validation data files are generated and can be used as input for DeepCAT. AutoCAT creates a new directory which contains four output files for cancer and control training and validation sets. Users can optonally name this directory or use the default directory name “DeepCATInput.” Non-cancer training data are used as true negative samples and cancer training data can be used as true positive samples. Sequences are randomly shuffled and split such that 80% of sequences for each length are reserved for training and 20% are withheld for validation.

## Usage
The following code performs TCR immune repertoire classification:

```python AutoCAT.py -d inputDIR -s cluster_size -p cluster_purity```

| Command| Description |
| ----------- | ----------- |
|-h, --help| show this help message and exit |
| -d INPUTDIR, --inputDIR=INPUTDIR| Input directory containing subdirectories "Cancer" and "Control" |
| -s SIZE, --size=SIZE| Cluster size. Size must be an integer. If omitted, the default cluster size is 50. |
| -p PURITY, --purity=PURITY| Cluster purity. Purity must be a float ≥0 and ≤1. If omitted, the default cluster purity is 0.80.|
| -g True, --graph=True| If True, AutoCAT will produce a graph showing the number of available sequences for cluster purities 60% to 100% using cluster sizes ranging from 10 to 500. If omitted, the default value is False.|
|-o OUTPUT_DIR, --outputDIR=OUTPUT_DIR| Output directory for training and validation files. If not given, a directory named ‘DeepCATInput ’will be created to hold the training and validation data files. |
