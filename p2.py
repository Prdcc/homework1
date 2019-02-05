"""M345SC Homework 1, part 2
Enrico Ancilotto
01210716
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import time
import random 

def nsearchv2(L,P,target):
    """Input:
    L: list containing *N* sub-lists of length M. Each sub-list
    contains M numbers (floats or ints), and the first P elements
    of each sub-list can be assumed to have been sorted in
    ascending order (assume that P<M). L[i][:p] contains the sorted elements in
    the i+1th sub-list of L
    P: The first P elements in each sub-list of L are assumed
    to be sorted in ascending order
    target: the number to be searched for in L

    Output:
    Lout: A list consisting of Q 2-element sub-lists where Q is the number of
    times target occurs in L. Each sub-list should contain 1) the index of
    the sublist of L where target was found and 2) the index within the sublist
    where the target was found. So, Lout = [[0,5],[0,6],[1,3]] indicates
    that the target can be found at L[0][5],L[0][6],L[1][3]. If target
    is not found in L, simply return an empty list (as in the code below)
    """
    Lout = []
    N = len(L)
    if(N == 0):
        return Lout

    M = len(L[0])        
    
    found = False
    mid = 0
    for i in range(N):                   
        Lout += [[i,j] for j in range(P,M) if L[i,j] == target]

        start = 0
        end = P-1
        if(not found):          #if the previous iteration found the target there is a good chance that it'll appear close to the previous mid otherwise reset mid
            mid = (P-1)//2

        found = False
        while start <= end:     #binary search
            if target == L[i, mid]:
                Lout += [[i,mid]]
                found = True
                break
            elif target < L[i, mid]:
                end = mid-1
            else:
                start = mid+1
            mid = (start + end)//2  
        
        if(start < mid):        #check for targets before midpoint
            for j in range(mid-1, start-1, -1):
                if(target == L[i, j]):
                    Lout += [[i,j]]
                else:
                    break
        
        if(end > mid):          #check for targets after midpoint
            for j in range(mid+1, end + 1):
                if(target == L[i, j]):
                    Lout += [[i,j]]
                else:
                    break
    return Lout

def nsearch(L,P,target):
    """Input:
    L: list containing *N* sub-lists of length M. Each sub-list
    contains M numbers (floats or ints), and the first P elements
    of each sub-list can be assumed to have been sorted in
    ascending order (assume that P<M). L[i][:p] contains the sorted elements in
    the i+1th sub-list of L
    P: The first P elements in each sub-list of L are assumed
    to be sorted in ascending order
    target: the number to be searched for in L

    Output:
    Lout: A list consisting of Q 2-element sub-lists where Q is the number of
    times target occurs in L. Each sub-list should contain 1) the index of
    the sublist of L where target was found and 2) the index within the sublist
    where the target was found. So, Lout = [[0,5],[0,6],[1,3]] indicates
    that the target can be found at L[0][5],L[0][6],L[1][3]. If target
    is not found in L, simply return an empty list (as in the code below)
    """
    Lout = []
    N = len(L)
    if(N == 0):
        return Lout

    M = len(L[0])        
    
    for i in range(N):                   
        Lout += [[i,j] for j in range(P,M) if L[i,j] == target]

        start = 0
        end = P-1
        mid = 0

        while start <= end:     #binary search
            if target == L[i, mid]:
                Lout += [[i,mid]]
                break
            elif target < L[i, mid]:
                end = mid-1
            else:
                start = mid+1
            mid = (start + end)//2  
        
        if(start < mid):        #check for targets before midpoint
            for j in range(mid-1, start-1, -1):
                if(target == L[i, j]):
                    Lout += [[i,j]]
                else:
                    break
        
        if(end > mid):          #check for targets after midpoint
            for j in range(mid+1, end + 1):
                if(target == L[i, j]):
                    Lout += [[i,j]]
                else:
                    break
    return Lout

def nsearchv3(L,P,target):
    """Input:
    L: list containing *N* sub-lists of length M. Each sub-list
    contains M numbers (floats or ints), and the first P elements
    of each sub-list can be assumed to have been sorted in
    ascending order (assume that P<M). L[i][:p] contains the sorted elements in
    the i+1th sub-list of L
    P: The first P elements in each sub-list of L are assumed
    to be sorted in ascending order
    target: the number to be searched for in L

    Output:
    Lout: A list consisting of Q 2-element sub-lists where Q is the number of
    times target occurs in L. Each sub-list should contain 1) the index of
    the sublist of L where target was found and 2) the index within the sublist
    where the target was found. So, Lout = [[0,5],[0,6],[1,3]] indicates
    that the target can be found at L[0][5],L[0][6],L[1][3]. If target
    is not found in L, simply return an empty list (as in the code below)
    """
    Lout = []
    N = len(L)
    if(N == 0):
        return Lout

    M = len(L[0])        
    
    for i in range(N):                   
        for j in range(P, M):
            if(L[i,j] == target):
                Lout += [[i,j]]

        start = 0
        end = P-1
        mid = 0

        while start <= end:     #binary search
            if target == L[i, mid]:
                Lout += [[i,mid]]
                break
            elif target < L[i, mid]:
                end = mid-1
            else:
                start = mid+1
            mid = (start + end)//2  
        
        if(start < mid):        #check for targets before midpoint
            for j in range(mid-1, start-1, -1):
                if(target == L[i, j]):
                    Lout += [[i,j]]
                else:
                    break
        
        if(end > mid):          #check for targets after midpoint
            for j in range(mid+1, end + 1):
                if(target == L[i, j]):
                    Lout += [[i,j]]
                else:
                    break
    return Lout

def generateList(N,M,P,integer = False, maxInt = 100):
    L = np.random.rand(N,M)
    L[:,:P].sort()
    if(integer):
        L = (maxInt*L).astype(int)
    return L

def analyse(Ns,MPs,maxInt=100):
    df = pd.DataFrame(columns=['N','Execution time','M','P','G','M - P'])
    for N in Ns:
        for M,P in MPs:
            L = generateList(N,M,P,True,maxInt)
            target = random.randint(0,maxInt-1)
            startTime = time.time()
            Lout = nsearch(L,P,target)
            df = df.append({'N' : int(N),
                            'M' : int(M),
                            'P' : int(P),
                            'M - P': int(M-P),
                            'G' : int(len(Lout)),
                            'Execution time' : time.time() - startTime}, ignore_index=True)
    return df

def analyseG(NMP, maxVals):
    df = pd.DataFrame(columns=['N','Execution time','M','P','G','M - P'])
    N,M,P = NMP
    for maxVal in maxVals:
        L = generateList(N,M,P,True,maxVal)
        target = random.randint(0,maxVal-1)
        startTime = time.time()
        Lout = nsearch(L,P,target)
        df = df.append({'N' : int(N),
                        'M' : int(M),
                        'P' : int(P),
                        'M - P': int(M-P),
                        'G' : int(len(Lout)),
                        'Execution time' : time.time() - startTime}, ignore_index=True)
    return df

def plotNMP(dfN,dfM,dfP,dfMP, saveNumber=-1):
    f,axes = plt.subplots(1,3)
    sns.regplot(y='Execution time', x='N', data=dfN, ax=axes[0])
    sns.regplot(y='Execution time', x='M', data=dfM, ax=axes[1])
    sns.regplot(y='Execution time', x='P', data=dfP, ax=axes[2])
    f.suptitle("Execution time for various values of the parameters")
    plt.savefig("figures/fig"+str(saveNumber)+".png", bbox_inches='tight')
    plt.show()

    sns.lmplot(y='Execution time', x='M', hue='M - P', data=dfMP).fig.suptitle("Execution time as M-P is fixed and M is varied")
    plt.savefig("figures/fig"+str(saveNumber+1)+".png", bbox_inches='tight')
    plt.show()

def plotG(dfG, saveNumber=-1):
    sns.regplot(y='Execution time', x='G', data=dfG).set_title("Execution time against various number of matching elements")
    plt.savefig("figures/fig"+str(saveNumber)+".png", bbox_inches='tight')
    plt.show()

    





def nsearch_time(saveStart=-100):
    """Analyze the running time of nsearch.
    Add input/output as needed, add a call to this function below to generate
    the figures you are submitting with your codes.

    Discussion: (add your discussion here)
    """

    Ns = range(1,1000,10)
    MPs = [(1000,200)]
    dfN = analyse(Ns,MPs)
    print("Finished N")

    Ns = [500]
    Ms = range(100,1000,10)
    Ps = [0]*len(Ms)
    MPs = zip(Ms, Ps)
    dfM = analyse(Ns,MPs)
    print("Finished M")

    Ns = [500]
    Ps = range(0,1501,30)
    Ms = [2000]*len(Ps)
    MPs = list(zip(Ms,Ps))
    dfP = analyse(Ns,MPs)
    print("Finished P")

    Ns = [500]
    Ms = np.array(range(100,2001,50))
    Ps = Ms - Ms[0]
    MPs = list(zip(Ms,Ps))
    Ms = np.array(range(200,2001,50))
    Ps = Ms - Ms[0]
    MPs += list(zip(Ms,Ps))
    Ms = np.array(range(500,2001,50))
    Ps = Ms - Ms[0]
    MPs += list(zip(Ms,Ps))
    Ms = np.array(range(1000,2001,50))
    Ps = Ms - Ms[0]
    MPs += list(zip(Ms,Ps))
    dfMP = analyse(Ns,MPs)
    print("Finished MP")


    plotNMP(dfN,dfM,dfP,dfMP,saveStart)    
    saveStart += 2


    NMP = (200,500,500)
    maxVals = range(1,100)
    dfG = analyseG(NMP, maxVals)



    

    plotG(dfG, saveStart)    
    saveStart += 1
    return #Modify as needed


if __name__ == '__main__':

    #add call(s) to nsearch here which generate the figure(s)
    #you are submitting
    sns.set()
    #nsearch_time()