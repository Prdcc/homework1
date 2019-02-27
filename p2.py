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
    
    for i in range(N):             #iterate thorugh each sublist  
        Lout += [[i,j] for j in range(P,M) if L[i,j] == target] #add matches in unsorted part of the list

        start = 0
        end = P-1
        mid = 0

        while start <= end:     #binary search
            mid = (start + end)//2  
            if target == L[i, mid]:
                Lout += [[i,mid]]
                break
            elif target < L[i, mid]:
                end = mid-1
            else:
                start = mid+1


        for j in range(mid-1, start-1, -1):     #check for targets before midpoint
            if(target == L[i, j]):
                Lout += [[i,j]]
            else:       #no more matches as list is ordered
                break
        
        for j in range(mid+1, end + 1):     #check for targets after midpoint
            if(target == L[i, j]):
                Lout += [[i,j]]
            else:
                break
    return Lout

def generateList(N,M,P,integer = False, maxInt = 100):
    """
    Generates a list of partially ordered sublists as in the problem statement
    """
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
    sns.regplot(y='Execution time', x='N', data=dfN, ax=axes[0], robust=True)
    sns.regplot(y='Execution time', x='M', data=dfM, ax=axes[1], robust=True)
    sns.regplot(y='Execution time', x='P', data=dfP, ax=axes[2], robust=True)
    f.suptitle("Enrico Ancilotto - plotNMP\nExecution time for various values of the parameters")
    plt.savefig("fig"+str(saveNumber)+".png", bbox_inches='tight')
    plt.show()

    sns.lmplot(y='Execution time', x='M', hue='M - P', data=dfMP, robust=True).fig.suptitle("Enrico Ancilotto - plotNMP\nExecution time as M-P is fixed and M is varied")
    plt.savefig("fig"+str(saveNumber+1)+".png", bbox_inches='tight')
    plt.show()

def plotG(dfG, saveNumber=-1):
    plt.clf()
    sns.regplot(y='Execution time', x='G', data=dfG, robust=True).set_title("Enrico Ancilotto - plotG\nExecution time against various number of matching elements")
    plt.savefig("fig"+str(saveNumber)+".png", bbox_inches='tight')
    plt.show()




def nsearch_time(saveStart=-100):
    """Analyze the running time of nsearch.
    Add input/output as needed, add a call to this function below to generate
    the figures you are submitting with your codes.

    Discussion: (add your discussion here)
    The algorithm is quite simple: we treat each sublist independently, first each
    unsorted part is searched with a normal linear search. Then a binary search is 
    performed on the sorted part, if no match is found the code moves on to the 
    next row, otherwise elements before and after the match are added as long
    as they are also a match.
    There are two main time consuming operations: comparing elements against the 
    target and appending elements to the end of Lout. Define N, M, P as in nsearch,
    furthermore let G be the number of matches (ie len(Lout)), r=G/NM the proportion
    of matches, s the time for a comparison and access, u the time for appending 
    to Lout. Appending values to the array takes uG=urNM time. Each search is independent
    of all the others, the ordered part of the search takes sNlog_2 P on average for finding the
    first match and a further srPN for the other matches. Finally the unsorted part
    of the search takes sN(M-P). Adding everything up we get an estimate of
    sN(M-(1-r)P+log_2 P)+arNM. Working under the assumption that r is small (if 
    not the algorithm can be made more efficient as will be described later), this
    gives us an estimate of O(N(M-P)) for our running time.
    The figures seem to confirm this analysis: figure 1 shows a strong linear trend for 
    N, M, P with the last one being inversely correlated to the execution time.
    Figure 2 shows that M is much less important than M-P for execution time (the
    small increasing tend is given by the fact that r is non-identically 0).
    The final figure shows the effect that large values of G can have on the runtime:
    N, M and P where kept constant (with P=0), while the value of r was increased
    resulting in more matches. Since appending has a high cost, this had a major
    effect on the execution time.

    Finally, a few notes on various choices I made while implementing the algorithm.
    I decided to optimise this algorithm for the case where r is small which is
    more probable in real world applications, for example this means that I've 
    implemented a linear search after the target has been found in the ordered
    part of the list (see line 58) for finding the other possible matches. A binary
    search here would only make the program more complex, and in any case accessing
    and comparing an item in an array is a lot faster than appending item to Lout,
    so the relative gains are minor.
    Working on the assumption that r is small also led me to treat every binary
    search as independent from the others: I experimented in setting the starting 
    point for the search at the coordinate of the previous match but with no increase
    in speed. If the data was known to follow some sort of distribution some adjustments
    would probably be possible but this algorithm is meant for general purposes.

    """
    Ns = range(1,1000,10)       #test running times for various values of the parameters
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

    print("Plotting first two figures")
    plotNMP(dfN,dfM,dfP,dfMP,saveStart)    
    saveStart += 2
    
    NMP = (200,500,500)
    maxVals = range(1,100)
    dfG = analyseG(NMP, maxVals)
    print("Finished G")

    plotG(dfG, saveStart)    
    saveStart += 1
    return #Modify as needed


if __name__ == '__main__':

    #add call(s) to nsearch here which generate the figure(s)
    #you are submitting
    sns.set()
    nsearch_time(1)