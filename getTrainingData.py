import os, random
from optparse import OptionParser
import pandas as pd
from decimal import Decimal

def CommandLineParser():
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="inputFile", default="", help="Input cluster file of clustered CDR3 sequences")
    parser.add_option("-l", "--labels", dest="labelFile", default="", help="Labeled CSV with sample IDs and donor type")
    parser.add_option("-s", "--size", dest="clusterSize", default="50", help="Cluster size")
    parser.add_option("-p", "--purity", dest="clusterPurity", default="0.80", help="Cluster purity")
    parser.add_option("-o", "--outputDIR", dest="outDIR", default="TrainingData", help="Output directory for training and validation sets. If not given, a directory named 'TrainingData' will be created to hold the training and validation data files.")
    return parser.parse_args()

(opt,_) = CommandLineParser()
inputFilename = opt.inputFile
inputLabelname = opt.labelFile
inputSize = int(opt.clusterSize)
inputPurity = float(opt.clusterPurity)
outDIR = opt.outDIR

if (len(inputFilename)==0 or len(inputLabelname)==0):
    print("Must provide an input cluster file and a labeled CSV with sample ID and donor type.")
    exit()

# clusterPurityDict stores value (+1 malignant, +0 benign) for all seq. in a cluster
clusterPurityDict = {}
# clusterSizeDict stores total number of sequences in a cluster
clusterSizeDict = {}

# Collects information about size and purity for each cluster
inputFile = open(inputFilename)
inputLabels = pd.read_csv(inputLabelname)

allBenign = inputLabels.index[inputLabels['Benign/Malignant'] == "Benign"].tolist()
allMalignant = inputLabels.index[inputLabels['Benign/Malignant'] == "Malignant"].tolist()

for fLine in inputFile.readlines():
    lineArr = list(fLine.split('\t'))
    sampleName = lineArr[-1][:-1]

    if (len(lineArr) > 1):
        clusterID = lineArr[1]

        if (clusterID not in clusterPurityDict):
            clusterPurityDict[clusterID] = 0
            clusterSizeDict[clusterID] = 0

        try:
            row_index = list(inputLabels["Sample"]).index(sampleName)
            # print (row_index)

            # Checks if sequence is from a cancer patient
            if (row_index in allMalignant):
                clusterPurityDict[clusterID] += 1
                clusterSizeDict[clusterID] += 1

            # Checks if sequence is from a healthy donor
            elif (row_index in allBenign):
                clusterPurityDict[clusterID] += 0
                clusterSizeDict[clusterID] += 1

            else:
                print (row_index)

        except:
            pass

# Separates sequences by length and type
cancerDict = {'cl12':[], 'cl13':[], 'cl14':[], 'cl15':[], 'cl16':[], 'cl17':[]}
nonCancerDict = {'nl12':[], 'nl13':[], 'nl14':[], 'nl15':[], 'nl16':[], 'nl17':[]}

# Filters clusters by size and purity
inputFile.seek(0)
for fLine in inputFile.readlines():
    lineArr = list(fLine.split('\t'))

    if (len(lineArr) > 1 and lineArr[1] in clusterPurityDict):
        seq = lineArr[0]
        clusterID = lineArr[1]

        if (clusterSizeDict[clusterID] <= inputSize):
            # Purity is calculated by the % of seq in a cluster that belong to healthy/cancer patients
            purity = clusterPurityDict[clusterID]/clusterSizeDict[clusterID]

            # Classify clusters as malignant
            if (purity >= inputPurity):
                if ('cl' + str(len(seq)) in cancerDict):
                    cancerDict['cl' + str(len(seq))].append(seq)

            # Classify clusters as healthy
            elif (purity <= float(Decimal('1')-Decimal(str(inputPurity)))):
                if ('nl' + str(len(seq)) in nonCancerDict):
                    nonCancerDict['nl' + str(len(seq))].append(seq)


def getTrainingandValidation(cancerDict, nonCancerDict):
    # Shuffle and distribute sequences for training and validations via 80/20 split
    MalignantTrain = []
    MalignantEval = []
    BenignTrain = []
    BenignEval = []

    for i in range(12, 18):
        random.Random(4).shuffle(cancerDict["cl" + str(i)])
        random.Random(4).shuffle(nonCancerDict["nl" + str(i)])

        c = cancerDict['cl' + str(i)]
        n = nonCancerDict["nl" + str(i)]

        train_cancer, val_cancer = c[: int(len(c) * .8)],  c[int(len(c) * .8):] 
        train_noncancer, val_noncancer = n[: int(len(n) * .8)],  n[int(len(n) * .8):] 

        MalignantTrain += train_cancer
        MalignantEval += val_cancer
        BenignTrain += train_noncancer
        BenignEval += val_noncancer

    # Create directory with training data files
    os.mkdir(outDIR)
    with open (outDIR + "/MalignantTrain.txt", "w") as f:
        for c in MalignantTrain:
            f.write(c + "\n")
    f.close() 

    with open (outDIR + "/MalignantEval.txt", "w") as f:
        for c in MalignantEval:
            f.write(c + "\n")
    f.close() 

    with open (outDIR + "/BenignTrain.txt", "w") as f:
        for n in BenignTrain:
            f.write(n + "\n")
    f.close()

    with open (outDIR + "/BenignEval.txt", "w") as f:
        for n in BenignEval:
            f.write(n + "\n")
    f.close()

getTrainingandValidation(cancerDict, nonCancerDict)
