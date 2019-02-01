"""M345SC Homework 1, part 1
Your name and CID here
"""
from collections import namedtuple
from Bio import SeqIO       #only used to read genome, delete before handing in


class KmerOccurrence:
    def __init__(self,kmer,pos):
        self.kmer = kmer
        self.locations = [pos]
        self.frequency = 1

    def addLocation(self, pos):
        self.locations += [pos]
        self.frequency += 1
    
    def literal(self):
        return toLiteral(self.kmer)

    def __str__(self):
        return("kmer="+self.literal()+"\nFrequency="+str(self.frequency))

def _try_composite(a, d, n, s):
    """
    code from: https://rosettacode.org/wiki/Miller–Rabin_primality_test
    """
    if pow(a, d, n) == 1:
        return False
    for i in range(s):
        if pow(a, 2**i * d, n) == n-1:
            return False
    return True # n  is definitely composite
 
def isPrime(n): # accurate as long as n<2'152'302'898'747
    """
    Implements Miller-Rabin primality test
    code from: https://rosettacode.org/wiki/Miller–Rabin_primality_test
    """
    if n in (0, 1):
        return True
    d, s = n - 1, 0
    while not d % 2:
        d, s = d >> 1, s + 1
    return not any(_try_composite(a, d, n, s) for a in (2, 3, 5, 7, 11))
    

def nextPrime(n):
    """
    Find next prime after 2n+1
    """
    index = 2*n+1
    while(not isPrime(index)):
        index += 2
    return index

def toNumeral(char):
    if(char == 'A'):
        return 0
    elif(char == 'T'):
        return 1
    elif(char == 'G'):
        return 2
    elif(char == 'C'):
        return 3

def toLiteral(intTuple):
    s=""
    for x in intTuple:
        if(x == 0):
            s += 'A'
        elif(x == 1):
            s += 'T'
        elif(x == 2):
            s += 'G'
        elif(x == 3):
            s += 'C'
    return s

def getHash(intTuple, q):
    currHash = 0
    power = 1
    for s in intTuple[::-1]:
        currHash = (currHash + power * s) % q
        power = (4*power) % q
    currHash = currHash % q
    return currHash

def ksearch(S,k,f,x):
    """
    Search for frequently-occurring k-mers within a DNA sequence
    and find the number of point-x mutations for each frequently
    occurring k-mer.
    Input:
    S: A string consisting of A,C,G, and Ts
    k: The size of k-mer to search for (an integer)
    f: frequency parameter -- the search should identify k-mers which
    occur at least f times in S
    x: the location within each frequently-occurring k-mer where
    point-x mutations will differ.

    Output:
    L1: list containing the strings corresponding to the frequently-occurring k-mers
    L2: list containing the locations of the frequent k-mers
    L3: list containing the number of point-x mutations for each k-mer in L1.

    Discussion: Add analysis here
    """
    assert(x <= f)
    length=len(S)
    if(length < k):
        return [],[],[]

    q = nextPrime(length)
    occurrences = {}     #initialise empty dictionary
    numeralS = [toNumeral(c) for c in S]

    bm = pow(4,k,q)
    currHash = getHash(numeralS[:k],q)

    occurrences[currHash] = [KmerOccurrence(numeralS[:k],0)]

    for index in range(1,length-k+1):
        if(index%100000 == 0):
            print("Progress: ",(index / length)*100,"%")
        currHash = (4*currHash- numeralS[index-1]*bm + numeralS[index-1+k]) % q
        kmer = numeralS[index : index+k]

        if currHash in occurrences:     #if hash already appears check if string has already appeared
            appeared = False
            for i,other in enumerate(occurrences[currHash]):
                if(kmer == other.kmer):
                    appeared = True
                    occurrences[currHash][i].addLocation(index)     #add location to already existing kmer
                    break
                
            if(not appeared):           #hash collision: add new kmer
                occurrences[currHash] += [KmerOccurrence(kmer,index)]
        else:
            occurrences[currHash] = [KmerOccurrence(kmer,index)]

    L1,L2,L3=[],[],[]
    for kmerOccurrences in occurrences.values():
        for kmerOccurrence in kmerOccurrences:
            if(kmerOccurrence.frequency >= f):
                L1.append(kmerOccurrence.literal())
                L2 += [kmerOccurrence.locations]

                L3 += [0]
                for i in range(4):
                    if(i == kmerOccurrence.kmer[x]):
                        continue
                    kmerMutated = kmerOccurrence.kmer[:]
                    kmerMutated[x] = i
                    kHash = getHash(kmerMutated,q)
                    
                    if(kHash in occurrences):
                        for occurrence in occurrences[kHash]:
                            if(occurrence.kmer == kmerMutated):
                                L3[-1] += occurrence.frequency
                                break  
    return L1,L2,L3




if __name__=='__main__':
    #Sample input and function call. Use/modify if/as needed
    S=""
    length=1000000
    for seq_record in SeqIO.parse("genomes/sequence.fasta", "fasta"):
        S=seq_record.seq[:]
    #print(S)
    k=20
    x=0
    f=2
    L1,L2,L3=ksearch(S,k,f,x)
