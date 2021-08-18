
# AutoCAT: Automated Cancer-Associated TCR discovery

AutoCAT is a computational method to predict tumor-associated TCRs from targeted TCR-seq data. The method utilizes [GIANA](https://github.com/s175573/GIANA) to quickly cluster similar CDR3 sequences. These clusters are then filtered to output training and validation data that can be used to train [DeepCAT](https://github.com/s175573/DeepCAT), a deep learning algorithm that identifies cancer associated beta chain TCRs. AutoCAT acts as a bridge to connect GIANA and DeepCAT.

AutoCAT is written in Python3 and requires the following dependencies to be installed:

* [biopython](https://biopython.org/)

* [faiss](https://github.com/facebookresearch/faiss)

* [scikit-learn](https://scikit-learn.org/stable/)

* [matplotlib](https://matplotlib.org/)

* [numpy](https://numpy.org/)

* [pandas](https://pandas.pydata.org/)

After installing these dependencies, please download the AutoCAT repository, which includes sample data as well as ```GIANA4.py```, the associated TRBV allele data ```Imgt_Human_TRBV.fasta```, and ```query.py```. Intermediate output files from AutoCAT, including the concatenated file of all input samples and cluster file, are provided for download on Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5176884.svg)](https://doi.org/10.5281/zenodo.5176884)  

## Standard Pipeline

1. Input files are concatenated into a file and a label file is generated to store where sequences are derived

2. GIANA is run on the concatenated file

3. Sequences are classified via cluster size and purity criteria

4. Training and validation data files are generated as input for DeepCAT

#### Running the full AutoCAT Pipeline

All example code is provided using the sample directory and files provided in the repository. To run the full AutoCAT pipeline, run a Python file with the following:

```
from AutoCAT import *
runAutoCAT(“trainingData/”)
```

#### Running individual AutoCAT functions

Running GIANA is the most time consuming aspect of AutoCAT. GIANA only needs to be run once for one set of training data. Users can concatenate all input data files and run GIANA on the output concatenated file with the following:

```
from AutoCAT import *
getInputFiles(“trainingData/”)
runGIANA(“trainingData.txt”)
```

Once the cluster file and labels CSV have been generated, training and validation data can be produced using ```getTrainingandValidation()```, which calls ```getClusterComposition()``` and ```separateClusters()``` automatically to generate the files needed to run DeepCAT. The following will generate training and validation files using a cluster size threshold of 50 and a cluster purity of 80%.

```
from AutoCAT import *
getTrainingandValidation("trainingData--RotationEncodingBL62.txt", "labels.csv", 50, 0.8)
```  

To generate diagnostic plots), use the following. ```graphAvailableSeq()``` will produce a plot saved as ```graphAvailableSequences.png```, and ```graphSamplePurity()``` will produce a plot saved as ```graphPurityBarGraph.png```.

```
from AutoCAT import *
graphAvailableSeq("trainingData--RotationEncodingBL62.txt", "labels.csv")
graphSamplePurity("trainingData--RotationEncodingBL62.txt", "labels.csv", 50)

```  
Once the training and validation data files have been produced, users can follow the instructions found on the [DeepCAT GitHub](https://github.com/s175573/DeepCAT#training-deepcat-models) to train a model on the data. A model trained with the provided data has been included in the directory “DeepCAT_CHKP.”

## Selection of AutoCAT Parameters

Three measures were utilized to determine the optimal AutoCAT parameters: available TCR sequences, segregation of patient samples, and TCR classification error. We selected the cluster size threshold based on the “elbow” or steep decrease in the number of available sequences. We then decided on a cluster purity based on the segregation of sequences from cancer and non-cancer patient samples; we selected 80% as default, as we believe the more stringent criteria will lead to a greater specific classification of TCR sequences. Lastly, we evaluated our selection for cluster purity using TCR classification error, or the percentage of healthy control sequences classified as cancer. In contrast to the less stringent 60% and 70% purity cutoffs, the 80% cutoff had significantly lower TCR classification error.

## Selection of Training Data

We recommend that users utilize a higher ratio of healthy to cancer sequences. We explored AutoCAT’s performance using subsets of the data, using only lung and melanoma cohorts and all healthy patient data. Both models performed similarly to the RNA sequencing model and the TCR model trained with the full set of career samples. Based on this independent validation, DeepCAT’s performance is robust using a smaller sample size.


## Usage

|Function|Description|
|-----------|-----------|
|```runAutoCAT(inputDIR, userSize=None, userPurity=None)```| Runs the full AutoCAT pipeline. AutoCAT requires a directory with data files divided into two subdirectories, ```“Cancer”``` and ```“Control"```. View the ```“trainingData/”``` directory included in the repository as an example. This function calls the functions described below to generate training and validation files that cane be used as input for DeepCAT. AutoCAT parameters can be customized by setting values for ```userSize``` and ```userPurity```. ```userSize``` can be set to an integer value >0, and userPurity can be set to a floating point value >0 and ≤1. If omitted, the default ```userSize``` and ```userPurity``` values are 50 and 0.80, respectively.|
|```getInputFiles(inputDIR)```| Outputs a concatenated text file named ```“trainingData.txt”``` of all files in the input directory as well as a CSV files called ```labels.csv``` that labels each input file as cancer or non-cancer. The input diretory must follow the format described above.|
|```runGIANA(inputFilename)```| Runs GIANA to produce a cluster file ```“trainingData--RotationEncodingBL62.txt”``` using a 3.3 Smith-Waterman alignment score and 32 set as the number of threads. The input file should be a concatenated file produced by ```getInputFiles()```|
|```getClusterComposition(clusterFilename, labelsFilename)```| Returns dictionaries ```clusterSizeDict``` and ```clusterPurityDict``` that store information about each cluster ID and the total number of sequences within a cluster and where the number of sequences derived from cancer patients, respectively. These dictionaries are later used to classify sequences and to generate the diagnostic purity plots. Requires the cluster file generated by GIANA and the CSV of labels generated by ```getInputFiles()``` as input.|
|```separateClusters(clusterFilename, clusterPurityDict, clusterSizeDict, userSize=None, userPurity=None)```| Returns dictionaries of sequences classified as cancer or non-cancer separated into lengths from 12 to 17 as well as the total number of sequences available for training. Requires the cluster file generated by GIANA and the two dictionaries generated by ```getClusterComposition()```. Users can optionally customize ```userSize``` and ```userPurity```.|
|```getTrainingandValidation(clusterFilename, labelsFilename, userSize=None, userPurity=None)```| Training and validation data files are generated and can be used as input for DeepCAT. AutoCAT creates a new directory called “DeepCATInput” which contains four output files for cancer and control training and validation sets. Non-cancer training data are used as true negative samples and cancer training data can be used as true positive samples. Sequences are randomly shuffled and split such that 80% of sequences for each length are reserved for training and 20% are withheld for validation. Requires the cluster file generated by GIANA and the CSV of labels generated by ```getInputFiles()```. Users can optionally customize ```userSize``` and ```userPurity.```|
|```graphAvailableSeq(clusterFilename, labelsFilename)```| Produces a graph showing the number of available sequences for cluster purities 60% to 100% using cluster sizes ranging from 10 to 500 (graphAvailableSequences.png). |
|```graphSamplePurity(clusterFilename, labelsFilename, userSize=None)```| Produces a line plot of the TCR classification error (graphTCRClassificationError.png), or the percent of healthy control TCRs classified as cancer. A bar plot displaying the percentage of sequences classified as cancer derived from cancer and non-cancer samples (graphPurityBarGraph) is also generated. If a cluster size is not provided, purity is determined utilized without using a cluster size threshold.|