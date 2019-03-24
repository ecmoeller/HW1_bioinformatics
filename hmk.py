import math
import plotly
plotly.tools.set_credentials_file(username='moell162', api_key='XWeDQVjHzfvEN01WsE7p')

import plotly.plotly as py
import plotly.tools as tls
import plotly.graph_objs as go
from plotly.tools import FigureFactory as FF

import matplotlib.pyplot as plt
import pylab as plt2
import numpy as np
import scipy.stats


import random
from operator import itemgetter 

def q3pa1(array):
    fig = plt.figure()
    
    n = array
    bins = 4
    plt.hist(n, bins)
    plt.title("All probe data")
    plt.xlabel("Expression levels")
    plt.ylabel("Probe Counts")

    #plotly_fig = tls.mpl_to_plotly( fig )
    #py.plot(plotly_fig, filename='all_probe_data_before_log')


def q3pa(totalLogVals):
    #HISTOGRAM FOR QUESTION 2 PART A
    fig = plt.figure()

    x = 0.25*random.randint(1,1000)
    print("What is x", x)
    #y = 0.3*random.randint(1,1000)
    print("Honestly what's in the first row of totalLogVals", totalLogVals[0])
    n = totalLogVals
    bins = 30
    plt.hist(n, bins)
    plt.title("All probe data after log transform")
    plt.xlabel("Expression levels")
    plt.ylabel("Probe Counts")

    #plotly_fig = tls.mpl_to_plotly( fig )
    #py.plot(plotly_fig, filename='all_probe_data')

def q3pb(first, second, third, fourth):
    print("Part2b")

    #GSM36777
    fig = plt.figure()
    n = first
    bins = 30
    plt.hist(n, bins)
    plt.title("GSM36777")
    plt.xlabel("Expression levels")
    plt.ylabel("Probe Counts")

    #plotly_fig = tls.mpl_to_plotly( fig )
    #py.plot(plotly_fig, filename='log_transform_GSM36777')

    #GSM36778
    fig = plt.figure()
    n = second
    bins = 28
    plt.hist(n, bins)
    plt.title("GSM36778")
    plt.xlabel("Expression levels")
    plt.ylabel("Probe Counts")

    plotly_fig = tls.mpl_to_plotly( fig )
    py.plot(plotly_fig, filename='log_transform_GSM36778')

    #GSM36779
    fig = plt.figure()
    n = third
    bins = 30
    plt.hist(n, bins)
    plt.title("GSM36779")
    plt.xlabel("Expression levels")
    plt.ylabel("Probe Counts")

    #plotly_fig = tls.mpl_to_plotly( fig )
    #py.plot(plotly_fig, filename='log_transform_GSM36779')

    #GSM36780
    fig = plt.figure()
    n = fourth
    bins = 29
    plt.hist(n, bins)
    plt.title("GSM36780")
    plt.xlabel("Expression levels")
    plt.ylabel("Probe Counts")

    #plotly_fig = tls.mpl_to_plotly( fig )
    #py.plot(plotly_fig, filename='log_transform_GSM36780')

def transpose(m):
    return [[m[j][i] for j in range(len(m))] for i in range(len(m[0]))]

def q3pc(logVals):
    #quantile normalization
    print("In q2pc")

    
    means = computeMeanAcrossProbes(logVals)

    #Sort means with indices key = lambda test_list: test_list[1]
    means.sort()
    
    #Attach i to every probe value
    probesWI = []
    for i in range(len(logVals)):
        anotherTemp = []
        for j in range(len(logVals[i])):
            betterTemp = []
            betterTemp.append(logVals[i][j])
            betterTemp.append(i)
            anotherTemp.append(betterTemp)
        probesWI.append(anotherTemp)
    
    #Sort probes per gene with indices
    #This involves looping through each column so we must transpose the matrix
    
    tLogVals = transpose(probesWI)
    sortedTLogVals = []
    for row in tLogVals:
        tempSort = row
        tempSort.sort()
        sortedTLogVals.append(tempSort)

    #map mean to proper sorted index from above
    for i in range(len(sortedTLogVals)):
        for j in range(len(sortedTLogVals[i])):
            index = sortedTLogVals[i][j][1]
            tLogVals[i][index][0] = means[j]
        

    #transpose back 
    newLogVals = transpose(tLogVals)

    column = []
    for i in range(len(newLogVals)):
        column.append(newLogVals[i][0][0])

    #Actually plotting the quant norm data
    fig = plt.figure()
    n = column
    bins = 30
    plt.hist(n, bins)
    plt.title("Quantile Normalization")
    plt.xlabel("Expression levels")
    plt.ylabel("Probe Counts")

    #plotly_fig = tls.mpl_to_plotly( fig )
    #py.plot(plotly_fig, filename='log_transform_GSM36780')

    return newLogVals

    

def q4paT(logVals, status, probeVals): #it's not actually logVals it's quan norm
    #t-test and rank sum
    print("In q4paT")
    #print("Here are the means", means)
    #t-test: per-probe significance level of p < 0.05
    # You should divide the gene expression data into two groups (metastasis vs. nonmetastasis),
    #  using the relapse variable in the clinical data, and test each probe independently

    #print("Log values from 1851", logVals[1851])
    
    pVals = []
    for i in range(len(logVals)):
        noRelapseL = []
        relapseL = []

        for j in range(len(logVals[i])):
            if(int(status[j+2]) == 0):
                noRelapseL.append(logVals[i][j][0])
            else:
                relapseL.append(logVals[i][j][0])

        twosample_results = scipy.stats.ttest_ind(relapseL, noRelapseL)
        indexArray = []
        indexArray.append(twosample_results[1][0][0])
        indexArray.append(i)
        pVals.append(indexArray)

    pVals = sorted(pVals, key = itemgetter(0))
    count = 0
    print("Printing the pVals for t test")
    while count < 10:
        print(pVals[count])
        print("name and id?",probeVals[pVals[count][1]][0],probeVals[pVals[count][1]][1])
        count += 1

    return pVals

def q4paR(logVals, status, probeVals):
    print("In q4paR")

    pVals = []
    for i in range(len(logVals)):
        noRelapseL = []
        relapseL = []

        for j in range(len(logVals[i])):
            if(int(status[j+2]) == 0):
                noRelapseL.append(logVals[i][j][0])
            else:
                relapseL.append(logVals[i][j][0])

        relapseL.sort()
        noRelapseL.sort()
        twosample_results = scipy.stats.ranksums(relapseL, noRelapseL)
        indexArray = []
        indexArray.append(twosample_results[1])
        indexArray.append(i)
        pVals.append(indexArray)

    pVals = sorted(pVals, key = itemgetter(0))
    count = 0
    print("Printing the pVals for rank sum")
    while count < 10:
        print(pVals[count])
        print("name and id?",probeVals[pVals[count][1]][0],probeVals[pVals[count][1]][1])
        count += 1

    return pVals

def q4pb(pValsT, pValsR):

    #count all p vals less than .05
    #make array of just p values
    realPVals = []
    count = 0
    for i in range(len(pValsT)):
        if pValsT[i][0] < .05:
            count += 1
        #math.log(float(probeVals[i][j]), 10)
        realPVals.append(-1*math.log(pValsT[i][0], 10))
    
    print("This is the count of p vals < .05 for t test" , count)

    #make histogram of all pvals for t test
    fig = plt.figure()
    n = realPVals
    bins = 29
    plt.hist(n, bins)
    plt.title("All p-values for t test")
    plt.xlabel("p-values")
    plt.ylabel("p-value counts")

    #plotly_fig = tls.mpl_to_plotly( fig )
    #py.plot(plotly_fig, filename='all_pvalues_ttest')


    realPVals = []
    count = 0
    for i in range(len(pValsR)):
        if pValsR[i][0] < .05:
            count += 1
        #math.log(float(probeVals[i][j]), 10)
        realPVals.append(-1*math.log(pValsR[i][0], 10))
    
    print("This is the count of p vals < .05 for rank sum" , count)

    #make histogram of all pvals for t test
    fig = plt.figure()
    n = realPVals
    bins = 29
    plt.hist(n, bins)
    plt.title("All p-values for rank sum")
    plt.xlabel("p-values")
    plt.ylabel("p-value counts")

    plotly_fig = tls.mpl_to_plotly( fig )
    py.plot(plotly_fig, filename='all_pvalues_ranksum')



def q4pc():
    # finding the probes in common with t test and rank-sum
    print("In q4pc")

def q5pa():
    #Bonferroni correction with the rank-sum test
    print("In q5pa")

def computeMeanAcrossProbes(logVals):
    print("Computing means")
    means = []
    for i in range(len(logVals)):
        sum = 0
        for j in range(len(logVals[i])):
            sum += logVals[i][j]
        mean = sum / len(logVals[i])
        meanTemp = []
        meanTemp.append([mean, i])
        means.append(meanTemp)

    return means

def main():
    with open('wang_data.txt','r') as f: 
        data = f.readlines()

    count = 0
    sampleids = []
    status = []
    probeVals = []
    geneList = []
    for line in data:
        line = line.strip()
        if (count == 0):
            sampleids = line.split('\t')
        if(count == 1):
            status = line.split('\t')
        if(count > 1):
            vals = line.split('\t')
            probeVals.append(vals)
            if vals[1] not in geneList:
                geneList.append(vals[1])
        count += 1

    noRelapse = 0
    relapse = 0
    for i in status:
        if(i == '0'):
            noRelapse += 1
        if(i == '1'):
            relapse += 1



    #print (sampleids)
    #print(status)
    print("Num of patients ", len(status)) 
    print("Num of no relapse ", noRelapse) 
    print("Num of relapse ", relapse) 
    #print(probeVals[0])
    print("Length of geneList ", len(geneList))
    print("Length of sampleIDs ", len(sampleids))
    print("Length of first row of probe vals", len(probeVals[0]))

    print("PROBE VALS 0", probeVals[0])
    

    logVals = []
    totalLogVals = [] # this is a long list of all probe data
    first = []
    second = []
    third = []
    fourth = []
    for i in range(len(probeVals)):
        trick = 0
        temp = []
        for j in range(len(probeVals[i])):
            #ignore the first two elements that are strings
            
            if(trick > 1):
                tempLogVal = math.log(float(probeVals[i][j]), 10)

                if(tempLogVal <= 0):
                    #This is creating the arrays for the individual histograms per gene in part 2b
                    if(trick == 2):
                        first.append(1.0)
                    if(trick == 3):
                        second.append(1.0)
                    if(trick == 4):
                        third.append(1.0)
                    if(trick == 5):
                        fourth.append(1.0)
                    temp.append(1.0)
                    totalLogVals.append(1.0)
                else:
                    #This is creating the arrays for the individual histograms per gene in part 2b
                    if(trick == 2):
                        first.append(tempLogVal)
                    if(trick == 3):
                        second.append(tempLogVal)
                    if(trick == 4):
                        third.append(tempLogVal)
                    if(trick == 5):
                        fourth.append(tempLogVal)
                    temp.append(tempLogVal)
                    totalLogVals.append(tempLogVal)
            trick += 1    
        logVals.append(temp)      

    probes = []
    for i in range(len(probeVals)):
        trick = 0
        for j in range(len(probeVals[i])):
            #ignore the first two elements that are strings
            
            if(trick > 1):
                probes.append(float(probeVals[i][j]))

            trick += 1
        

    q3pa1(probes)

    q3pa(totalLogVals)

    q3pb(first, second, third, fourth)

    quanNorm = q3pc(logVals)

    pValsT = q4paT(quanNorm, status, probeVals)

    pValsR = q4paR(quanNorm, status, probeVals)

    q4pb(pValsT,pValsR)

    q5pa()


if __name__ == '__main__':
    main()