import sys,os, random
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from decimal import Decimal

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
    labels_df = pd.DataFrame.from_dict(labelsDict)
    labels_df.to_csv("labels.csv", index=False)

def runGIANA(inputFilename):
    os.system("python GIANA4.py -f " + inputFilename + " -b -S 3.3 -N 32")

def getClusterComposition(clusterFilename, labelsFilename):
    allID = pd.read_csv(labelsFilename, header=None, index_col=0, squeeze=True).to_dict()
    # clusterPurityDict stores value (+1 malignant, +0 benign) for all seq. in a cluster
    clusterPurityDict = {}
    # clusterSizeDict stores total number of sequences in a cluster
    clusterSizeDict = {}

    # Collects information about size and purity for each cluster
    clusterFile = open(clusterFilename)
    
    for fLine in clusterFile.readlines()[1:]:
        lineArr = list(fLine.split('\t'))
        sampleName = lineArr[-1].rstrip("\n")

        if (len(lineArr) > 1):
            clusterID = lineArr[1]

            if (clusterID not in clusterPurityDict):
                clusterPurityDict[clusterID] = 0
                clusterSizeDict[clusterID] = 0

            clusterPurityDict[clusterID] += int(allID[sampleName])
            clusterSizeDict[clusterID] += 1

    return clusterSizeDict, clusterPurityDict

def separateClusters(clusterFilename, clusterSizeDict, clusterPurityDict, userSize=None, userPurity=None):
    if userSize == None:
        userSize = 50
    if userPurity == None:
        userPurity = 0.8

    clusterFile = open(clusterFilename)
    cancer = 0
    noncancer = 0
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
                        cancer += 1
                        availableSeq += 1

                # Classify clusters as healthy
                elif (purity <= float(Decimal('1')-Decimal(str(userPurity)))):
                    if ('nl' + str(len(seq)) in nonCancerDict):
                        nonCancerDict['nl' + str(len(seq))].append(seq)
                        noncancer += 1
                        availableSeq += 1

    return cancerDict, nonCancerDict, availableSeq

def getTrainingandValidation(clusterFilename, labelsFilename, userSize=None, userPurity=None):
    if userSize == None:
        userSize = 50
    if userPurity == None:
        userPurity = 0.8

    clusterSizeDict, clusterPurityDict = getClusterComposition(clusterFilename, labelsFilename)
    cancerDict, nonCancerDict, availableSeq = separateClusters(clusterFilename, clusterSizeDict, clusterPurityDict, userSize, userPurity)
    outDIR = "DeepCATInput"

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

def runAutoCAT(inputDIR, userSize=None, userPurity=None):
    if userSize == None:
        userSize = 50
    if userPurity == None:
        userPurity = 0.8

    getInputFiles(inputDIR)
    runGIANA("trainingData.txt")
    getTrainingandValidation("trainingData--RotationEncodingBL62.txt", "labels.csv", userSize, userPurity)

def graphAvailableSeq(clusterFilename, labelsFilename):
    availSeq = {0.6:[], 0.7:[], 0.8:[], 0.9:[], 1:[]}
    allSizes = [10, 50, 100, 200, 500]

    clusterPurityDict, clusterSizeDict = getClusterComposition(clusterFilename, labelsFilename)

    for size in allSizes:
        for purity in availSeq.keys():
            _, __, availableSeq = separateClusters(clusterFilename, clusterPurityDict, clusterSizeDict, size, purity)
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

def graphSamplePurity(clusterFilename, labelsFilename, userSize=None):
    if userSize == None:
        userSize = 50

    clusterSizeDict, clusterPurityDict = getClusterComposition(clusterFilename, labelsFilename)

    # [derived from noncancer, derived from cancer]
    cancerPurityDict = {0.6: [0, 0], 0.7: [0, 0], 0.8: [0, 0], 0.9: [0, 0]}
    nonCancerPurityDict = {0.6: [0, 0], 0.7: [0, 0], 0.8: [0, 0], 0.9: [0, 0]}

    for purity in cancerPurityDict:
        for clusterID in clusterPurityDict:
            # if (clusterSizeDict[clusterID] <= userSize):
            clusterPurity = clusterPurityDict[clusterID]/clusterSizeDict[clusterID]

            if (clusterPurity >= purity):
                cancerPurityDict[purity][0] += clusterSizeDict[clusterID]-clusterPurityDict[clusterID]
                cancerPurityDict[purity][1] += clusterPurityDict[clusterID]

            elif (clusterPurity <= float(Decimal('1')-Decimal(str(purity)))):
                nonCancerPurityDict[purity][0] += clusterSizeDict[clusterID] - clusterPurityDict[clusterID]
                nonCancerPurityDict[purity][1] += clusterPurityDict[clusterID]
    #Bar Graph
    cacData = {"Derived from\nCancer": [], "Derived from\nNonCancer": []}
    cancData = {"Derived from\nCancer": [], "Derived from\nNonCancer": []}

    for p in cancerPurityDict:
        cacData["Derived from\nCancer"].append( int(round(cancerPurityDict[p][1] / sum(cancerPurityDict[p]), 2) * 100))
        cacData["Derived from\nNonCancer"].append( int(round(cancerPurityDict[p][0] / sum(cancerPurityDict[p]), 2) * 100))
        cancData["Derived from\nCancer"].append( int(round(nonCancerPurityDict[p][1] / sum(nonCancerPurityDict[p]), 2) * 100))
        cancData["Derived from\nNonCancer"].append( int(round(nonCancerPurityDict[p][0] / sum(nonCancerPurityDict[p]), 2) * 100))

        # print (cancerPurityDict[p][0]/(cancerPurityDict[p][0] + nonCancerPurityDict[p][0]))
    
    classifiedAsCancer = pd.DataFrame(data=cacData, index=["60%", "70%", "80%", "90%"])
    classifiedAsNonCancer = pd.DataFrame(data=cancData,index=["60%", "70%", "80%", "90%"])

    # alldf = [classifiedAsCancer, classifiedAsNonCancer]
    alldf = [classifiedAsCancer]

    numDfs = len(alldf)
    numCols = len((alldf[0].columns))
    numRows = len(alldf[0].index)

    plt.figure(figsize=(7,5))
    barplot = plt.subplot()

    # Make a bar plot for each dataframe
    for df in alldf :
        barplot = df.plot(kind="bar", linewidth=0, stacked=True, ax=barplot, legend=False, grid=False)

    h,l = barplot.get_legend_handles_labels()
    for i in range(0, numDfs * numCols, numCols):
        for j, pa in enumerate(h[i:i+numCols]):
            for rect in pa.patches:
                rect.set_x(rect.get_x() + 1 / float(numDfs + 1) * i / float(numCols))
                rect.set_hatch("/" * int(i / numCols))   
                rect.set_width(1 / float(numDfs + 1))

    barplot.set_xticks((np.arange(0, 2 * numRows, 2) + 1 / float(numDfs + 1)) / 2.)
    barplot.set_xticklabels(df.index, rotation = 0)
    barplot.set_title("TCR Classification Across Different Purities")
    barplot.set_xlabel("Purity Criteria")
    barplot.set_ylabel("Fraction of Sequences Wrongly Assigned")

    n=[]        
    for i in range(numDfs):
        n.append(barplot.bar(0, 0, color="gray", hatch="/" * i))

    l1 = barplot.legend(h[:numCols], l[:numCols], loc=[1, 0.75], prop={'size': 10})
    # l2 = plt.legend(n, ["Classified\nCancer", "Classified\nNoncancer"], loc=[1, 0.5], prop={'size': 10}) 
    barplot.add_artist(l1)

    for rec in barplot.patches:
        height = rec.get_height()
        if height != 0:
            barplot.text(rec.get_x() + rec.get_width() / 2,
                         rec.get_y() + height / 2,
                         "{}%".format(int(height)),
                         ha='center',
                         va='bottom')

    plt.tight_layout()
    plt.savefig("graphPurityBarGraph.png", bbox_extra_artists=(l1,), bbox_inches='tight')
