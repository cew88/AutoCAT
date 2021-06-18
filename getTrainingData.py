import os, sys, random
import pandas as pd
from decimal import Decimal


# clusterPurityDict stores value (+1 malignant, +0 benign) for all seq. in a cluster
clusterPurityDict = {}
# clusterSizeDict stores total number of sequences in a cluster
clusterSizeDict = {}

try: 
    inputFilename = sys.argv[1]
    inputLabels = pd.read_csv(sys.argv[2])
except:
    print ("Please provide input files.")
    exit()

try:
    inputSize = int(sys.argv[3])
    inputPurity = float(sys.argv[4])
except:
    # Default: use cluster size <= 50 and cluster purity >= 80% if no arguments are provided
    inputSize = 50
    inputPurity = 0.8


allBenign = inputLabels.index[inputLabels['Benign/Malignant'] == "Benign"].tolist()
allMalignant = inputLabels.index[inputLabels['Benign/Malignant'] == "Malignant"].tolist()

# Collects information about size and purity for each cluster
inputFile = open(inputFilename)
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
            elif (purity <= float(Decimal('1')-Decimal(str(inputPurity)))
                ):
                if ('nl' + str(len(seq)) in nonCancerDict):
                    nonCancerDict['nl' + str(len(seq))].append(seq)


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
# newDir = "TrainingData"
# os.mkdir(newDir)
# with open (newDir + "/MalignantTrain.txt", "w") as f:
#     for c in MalignantTrain:
#         f.write(c + "\n")
# f.close() 

# with open (newDir + "/MalignantEval.txt", "w") as f:
#     for c in MalignantEval:
#         f.write(c + "\n")
# f.close() 

# with open (newDir + "/BenignTrain.txt", "w") as f:
#     for n in BenignTrain:
#         f.write(n + "\n")
# f.close()

# with open (newDir + "/BenignEval.txt", "w") as f:
#     for n in BenignEval:
#         f.write(n + "\n")
# f.close()
