import sys,os, random
from optparse import OptionParser
import matplotlib.pyplot as plt
import pandas as pd
from decimal import Decimal

def CommandLineParser():
    parser = OptionParser()
    parser.add_option("-d", "--inputDIR", dest="inputDIR", default="", help="Input directory containing subdirectories \"Cancer\" and \"Control\"")
    parser.add_option("-s", "--size", dest="clusterSize", default="50", help="Cluster size")
    parser.add_option("-p", "--purity", dest="clusterPurity", default="0.80", help="Cluster purity")
    parser.add_option("-g", "--graph", dest="generateGraph", default="False", help="If True, AutoCAT will produce a graph showing the number of available sequences for cluster purities 60% to 100% using cluster sizes ranging from 10 to 500.")
    parser.add_option("-o", "--outputDIR", dest="outDIR", default="DeepCATInput", help="Output directory for training and validation sets. If not given, a directory named 'DeepCATInput' will be created to hold the training and validation data files.")
    return parser.parse_args()

(opt,_) = CommandLineParser()
inputDIR = opt.inputDIR
inputSize = int(opt.clusterSize)
inputPurity = float(opt.clusterPurity)
outDIR = opt.outDIR
generateGraph = opt.generateGraph


if not os.path.isdir(inputDIR):
    print ("Must provide a valid input directory.")
    exit()

def getInputFiles(inputDIR):
    subDIRList = [f.path for f in os.scandir(inputDIR) if f.is_dir()]
    allID = {}
    label = None

    with open("trainingData.txt", "w") as p:
        for folderPath in subDIRList:
            if ("Cancer" in folderPath):
                label = 1
            elif ("Control" in folderPath):
                label = 0

            for filename in os.listdir(folderPath):
                file = open(folderPath + "/" + filename)
                for fLine in file.readlines()[1:]:
                    p.write(fLine)

                    lineArr = list(fLine.split('\t'))
                    if lineArr[-1] not in allID:
                        allID[lineArr[-1].rstrip("\n")] = label
    p.close()

    labelsDict = {"Sample": list(allID.keys()), "Benign/Malignant": list(allID.values())}
    # labels_df = pd.DataFrame.from_dict(labelsDict)
    # labels_df.to_csv("labels.csv", index=False)
    return allID

def runGIANA():
    os.system("python GIANA4.py -f " + "trainingData.txt" + " -b -S 3.3 -N 32")

def getClusterComposition(allID):
    # clusterPurityDict stores value (+1 malignant, +0 benign) for all seq. in a cluster
    clusterPurityDict = {}
    # clusterSizeDict stores total number of sequences in a cluster
    clusterSizeDict = {}

    # Collects information about size and purity for each cluster
    clusterFile = open("trainingData--RotationEncodingBL62.txt")
    
    for fLine in clusterFile.readlines()[1:]:
        lineArr = list(fLine.split('\t'))
        sampleName = lineArr[-1].rstrip("\n")

        if (len(lineArr) > 1):
            clusterID = lineArr[1]

            if (clusterID not in clusterPurityDict):
                clusterPurityDict[clusterID] = 0
                clusterSizeDict[clusterID] = 0

            try:
                clusterPurityDict[clusterID] += allID[sampleName]
                clusterSizeDict[clusterID] += 1

            except:
                pass

    return clusterPurityDict, clusterSizeDict

def separateClusters(clusterPurityDict, clusterSizeDict, userSize, userPurity):

    clusterFile = open("trainingData--RotationEncodingBL62.txt")
    availableSeq = 0

    # Separates sequences by length and type
    cancerDict = {'cl12':[], 'cl13':[], 'cl14':[], 'cl15':[], 'cl16':[], 'cl17':[]}
    nonCancerDict = {'nl12':[], 'nl13':[], 'nl14':[], 'nl15':[], 'nl16':[], 'nl17':[]}

    # Filters clusters by size and purity
    clusterFile.seek(0)
    for fLine in clusterFile.readlines():
        lineArr = list(fLine.split('\t'))

        if (len(lineArr) > 1 and lineArr[1] in clusterPurityDict):
            seq = lineArr[0]
            clusterID = lineArr[1]

            if (clusterSizeDict[clusterID] <= userSize):
                # Purity is calculated by the % of seq in a cluster that belong to healthy/cancer patients
                purity = clusterPurityDict[clusterID]/clusterSizeDict[clusterID]

                # Classify clusters as malignant
                if (purity >= userPurity):
                    if ('cl' + str(len(seq)) in cancerDict):
                        cancerDict['cl' + str(len(seq))].append(seq)
                        availableSeq += 1

                # Classify clusters as healthy
                elif (purity <= float(Decimal('1')-Decimal(str(userPurity)))):
                    if ('nl' + str(len(seq)) in nonCancerDict):
                        nonCancerDict['nl' + str(len(seq))].append(seq)
                        availableSeq += 1

    return cancerDict, nonCancerDict, availableSeq

def graphAvailableSeq(allID):
    availSeq = {0.6:[], 0.7:[], 0.8:[], 0.9:[], 1:[]}
    allSizes = [10, 50, 100, 200, 500]

    for size in allSizes:
        for purity in availSeq.keys():
            _, __, availableSeq = separateClusters(clusterPurityDict, clusterSizeDict, size, purity)
            availSeq[purity].append(availableSeq)

    fig, ax = plt.subplots()
    plt.rcParams["figure.figsize"] = (10,10)

    plt.rc('legend',fontsize=12)
    plt.axvline(x=50, linestyle='--', alpha=0.5, color='k')

    ax.plot(allSizes, availSeq[0.6], label="60% Purity")
    ax.plot(allSizes, availSeq[0.7], label="70% Purity")
    ax.plot(allSizes, availSeq[0.8], label="80% Purity")
    ax.plot(allSizes, availSeq[0.9], label="90% Purity")
    ax.plot(allSizes, availSeq[1], label="100% Purity")
    ax.legend()
    ax.set_xlim(500, 0)  # Decreasing x-axis
    ax.set_xlabel('Cluster Size', fontsize=16)
    ax.set_ylabel('Number of Available Seq', fontsize=16)
    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)
    ax.grid(True)

    plt.savefig("graphAvailableSequences.png", bbox_inches='tight')

def getTrainingandValidation(cancerDict, nonCancerDict):
    # Shuffle and distribute sequences for training and validations via 80/20 split
    cancerTrain = []
    cancerEval = []
    controlTrain = []
    controlEval = []

    for i in range(12, 18):
        random.Random(4).shuffle(cancerDict["cl" + str(i)])
        random.Random(4).shuffle(nonCancerDict["nl" + str(i)])

        c = cancerDict['cl' + str(i)]
        n = nonCancerDict["nl" + str(i)]

        train_cancer, val_cancer = c[: int(len(c) * .8)],  c[int(len(c) * .8):] 
        train_noncancer, val_noncancer = n[: int(len(n) * .8)],  n[int(len(n) * .8):] 

        cancerTrain += train_cancer
        cancerEval += val_cancer
        controlTrain += train_noncancer
        controlEval += val_noncancer

    # Create directory with training data files
    os.mkdir(outDIR)
    with open (outDIR + "/CancerTrain.txt", "w") as f:
        for c in cancerTrain:
            f.write(c + "\n")
    f.close() 

    with open (outDIR + "/CancerEval.txt", "w") as f:
        for c in cancerEval:
            f.write(c + "\n")
    f.close() 

    with open (outDIR + "/ControlTrain.txt", "w") as f:
        for n in controlTrain:
            f.write(n + "\n")
    f.close()

    with open (outDIR + "/ControlEval.txt", "w") as f:
        for n in controlEval:
            f.write(n + "\n")
    f.close()


allID = getInputFiles(inputDIR)
runGIANA()
clusterPurityDict, clusterSizeDict = getClusterComposition(allID)
cancerDict, nonCancerDict, availableSeq = separateClusters(clusterPurityDict, clusterSizeDict, inputSize, inputPurity)
getTrainingandValidation(cancerDict, nonCancerDict)

if (generateGraph):
    graphAvailableSeq(allID)