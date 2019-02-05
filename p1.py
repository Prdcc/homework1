"""M345SC Homework 1, part 1
01210716
"""
from Bio import SeqIO       #only used to read genome
import sys

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
    Find next prime after n 
    """
    index = n +1 - (n%2) #make n odd if it's even
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
    raise ValueError("Nucleotide %c not recognised, should be one of A,T,G,C"%char)

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
    global powers
    for i in range(len(intTuple)):
        currHash = (currHash + powers[i]*intTuple[i]) % q
    return currHash



def ksearchv2(S,k,f,x):
    occurrences = {}     #initialise empty dictionary
    nucleotides = "ATGC"
    S = "".join([c for c in S if c in nucleotides])
    length=len(S)
    if(length < k):
        return [],[],[]
    print("starting search")
    occurrences[S[:k]] = [0]
    for index in range(1,length-k+1):
        if(index%100000 == 0):          #delete this later
            print("Progress: ",(index / length)*100,"%")

        kmer = S[index : index+k]

        if kmer in occurrences:     #if kmer has already appeared add location
            occurrences[kmer].append(index)
        else:
            occurrences[kmer] = [index]

    L1,L2,L3=[],[],[]
    for kmer in occurrences:
        frequency = len(occurrences[kmer])
        if(frequency >= f):
            L1.append(kmer)
            L2 += [occurrences[kmer]]

            L3.append(0)
            for i in "ATGC":
                if(i == kmer[x]):
                    continue
                kmerMutated = list(kmer)
                kmerMutated[x] = i
                
                if("".join(kmerMutated) in occurrences):
                    L3[-1] += len(occurrences[kmer])
    return L1,L2,L3

def getHashStr(S, q):
    return getHash([toNumeral(c) for c in S],q)

def ksearchstr(S,k,f,x):
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
    q = 472884059     #just over 2^60, a prime, chosen so that all computations are kept below 64 bit precision 
    global powers 
    powers = tuple(pow(4,i,q) for i in reversed(range(k)))
    bm = pow(4,k,q)
    occurrences = {}     #initialise empty dictionary
    nucleotides = "ATGC"
    S = [c for c in S if c in nucleotides]
    length=len(S)
    if(length < k):
        return [],[],[]
    print("Starting search")

    currHash = getHashStr(S[:k],q)
    occurrences[currHash] = [(S[:k],[0])]

    for index in range(1,length-k+1):
        if(index%100000 == 0):          #delete this later
            print("Progress: ",(index / length)*100,"%")
        currHash = (4*currHash- toNumeral(S[index-1])*bm + toNumeral(S[index-1+k])) % q
        kmer = S[index : index+k]

        if currHash in occurrences:     
            toAdd = True                                #check for collisions
            for i in range(len(occurrences[currHash])):
                if(kmer == occurrences[currHash][i][0]):
                    occurrences[currHash][i][1].append(index)
                    toAdd = False
                    break
            
            if(toAdd):
                occurrences[currHash].append((kmer,[index]))
        else:
            occurrences[currHash] = [(kmer,[index])]

    L1,L2,L3 = [],[],[]
    for results in occurrences.values():
        for kmer, locations in results:
            if(len(locations) >= f):
                kmerString = kmer
                L1.append(kmerString)
                L2.append(locations)
                L3.append(0)
                for i in "ATGC":
                    if(i == kmer[x]):
                        continue
                    kmerMutated = list(kmer[:])
                    kmerMutated[x] = i
                    kmerMutated = "".join(kmerMutated)
                    currHash = getHashStr(kmerMutated,q)
                    if(currHash in occurrences):
                        for i in range(len(occurrences[currHash])):
                            if(kmerMutated == occurrences[currHash][i][0]):
                                L3[-1] += len(occurrences[currHash][i][1])
                                break
    return L1,L2,L3




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
    q = 472884059     #just over 2^60, a prime, chosen so that all computations are kept below 64 bit precision 
    global powers 
    powers = tuple(pow(4,i,q) for i in reversed(range(k)))
    bm = pow(4,k,q)
    occurrences = {}     #initialise empty dictionary
    nucleotides = "ATGC"
    S = [toNumeral(c) for c in S if c in nucleotides]
    length=len(S)
    if(length < k):
        return [],[],[]
    print("Starting search")

    currHash = getHash(S[:k],q)
    occurrences[currHash] = [(S[:k],[0])]

    for index in range(1,length-k+1):
        if(index%100000 == 0):          #delete this later
            print("Progress: ",(index / length)*100,"%")
        currHash = (4*currHash- S[index-1]*bm + S[index-1+k]) % q
        kmer = S[index : index+k]

        if currHash in occurrences:     
            toAdd = True                                #check for collisions
            for i in range(len(occurrences[currHash])):
                if(kmer == occurrences[currHash][i][0]):
                    occurrences[currHash][i][1].append(index)
                    toAdd = False
                    break
            
            if(toAdd):
                occurrences[currHash].append((kmer,[index]))
        else:
            occurrences[currHash] = [(kmer,[index])]

    L1,L2,L3 = [],[],[]
    for results in occurrences.values():
        for kmer, locations in results:
            if(len(locations) >= f):
                kmerString = toLiteral(kmer)
                L1.append(kmerString)
                L2.append(locations)
                L3.append(0)
                for i in range(4):
                    if(i == kmer[x]):
                        continue
                    kmerMutated = kmer[:]
                    kmerMutated[x] = i
                    currHash = getHash(kmerMutated,q)
                    if(currHash in occurrences):
                        for i in range(len(occurrences[currHash])):
                            if(kmerMutated == occurrences[currHash][i][0]):
                                L3[-1] += len(occurrences[currHash][i][1])
                                break
    return L1,L2,L3



if __name__=='__main__':
    #Sample input and function call. Use/modify if/as needed
    S=""
    for seq_record in SeqIO.parse("genomes/sequence.fasta", "fasta"):
        S=seq_record.seq[:]
    S=S[:len(S)//4]    
    k=2
    f=1000
    x=0
    if(len(sys.argv) == 4):
        k = int(sys.argv[1])
        f = int(sys.argv[2])
        x = int(sys.argv[3])
    L1,L2,L3=ksearchstr(S,k,f,x)
