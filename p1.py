"""M345SC Homework 1, part 1
01210716
"""
from Bio import SeqIO       #only used to read genome
import sys          #used to access command-line arguments and byte size of objects

def toNumeral(char):
    """
    converts nucleotide to number
    """
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
    """
    converts numeric array representation of strand of DNA to string representation
    """
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
    """
    intTuple should be of the form [0,1,3,...]
    generates the hash of intTuple
    """
    currHash = 0
    global powers       #defined in ksearch{Num,Str}, each entry contains 4**j
    for i in range(len(intTuple)):
        currHash = (currHash + powers[i]*intTuple[i]) % q
    return currHash

def getHashStr(S, q):
    """
    Like getHash but S is assumed to be a string of nucleotides
    """
    currHash = 0
    global powers   #defined in ksearch{Num,Str}, each entry contains 4**j
    for i in range(len(S)):
        currHash = (currHash + powers[i]*toNumeral(S[i])) % q
    return currHash

def ksearchStr(S,k,f,x):
    """
    sea ksearchNum and discussion in ksearch
    """
    q = 472884059     #just over 2^30, a prime, chosen so that all computations are kept below 64 bit precision 
    global powers 
    powers = tuple(pow(4,i,q) for i in reversed(range(k)))
    bm = pow(4,k,q)
    occurrences = {}     #initialise empty dictionary
    nucleotides = "ATGC"
    S = "".join([c for c in S if c in nucleotides])
    length=len(S)
    if(length < k):
        return [],[],[]
    print("Starting search")

    currHash = getHashStr(S[:k],q)
    occurrences[currHash] = [[S[:k],0]]
    print(occurrences[currHash])



    for index in range(1,length-k+1):
        if(index%300000 == 0):          #delete this later
            print("------------------------------\nProgress: ",index)
            print("Size of dictionary: ", sys.getsizeof(occurrences))
            print("Size of random entry: ",realSize(occurrences[currHash]))
            print("Number of collisions in random entry: ", len(occurrences[currHash]))
            print("Number of entries: ",len(occurrences))
        currHash = (4*currHash- toNumeral(S[index-1])*bm + toNumeral(S[index-1+k])) % q
        kmer = S[index : index+k]

        position = findKmer(occurrences, kmer, currHash)

        if position == -2:
            occurrences[currHash] = [[kmer,index]]
        elif position == -1:
            occurrences[currHash].append([kmer,index])
        else:
            occurrences[currHash][position].append(index)

    L1,L2,L3 = [],[],[]
    for results in occurrences.values():
        for kmerInfo in results:
            if(len(kmerInfo) > f):
                kmerString = kmerInfo[0]
                locations = kmerInfo[1:]
                L1.append(kmerString)
                L2.append(locations)
                L3.append(0)
                for i in "ATGC":
                    if(i == kmerString[x]):
                        continue
                    kmerMutated = list(kmerString)
                    kmerMutated[x] = i
                    kmerMutated = "".join(kmerMutated)
                    currHash = getHashStr(kmerMutated,q)
                    if(currHash in occurrences):
                        for i in range(len(occurrences[currHash])):
                            if(kmerMutated == occurrences[currHash][i][0]):
                                L3[-1] += len(occurrences[currHash][i])-1
                                break
    return L1,L2,L3

def realSize(o, handlers={}, verbose=False):
    """ Code adapted from: http://code.activestate.com/recipes/577504-compute-memory-footprint-of-an-object-and-its-cont/
    Finds the total size of the object o, including the size of the objects it references
    """
    all_handlers = {tuple: iter,
                    list: iter,
                    set: iter,
                    frozenset: iter,
                   }
    all_handlers.update(handlers)     # user handlers take precedence
    seen = set()                      # track which object id's have already been seen
    default_size = sys.getsizeof(0)       # estimate sizeof object without __sizeof__

    def sizeof(o):
        if id(o) in seen:       # do not double count the same object
            return 0
        seen.add(id(o))
        s = sys.getsizeof(o, default_size)

        for typ, handler in all_handlers.items():
            if isinstance(o, typ):
                s += sum(map(sizeof, handler(o)))
                break
        return s

    return sizeof(o)

def findKmer(dict, kmer, kHash):
    """
    find kmer in the dictionary: returns the index if kmer is in dictionary
    otherwise it returns -1 in case of a hash collision and -2 if kHash is
    not already in dict
    """
    if kHash in dict:
        p = dict[kHash]
        for i in range(len(p)):
            if(p[i][0] == kmer):
                return i
        return -1
    else: return -2

def ksearchNum(S,k,f,x):
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

    Implementation of Rabin-Karp
    """
    q = 472884059     #just over 2^30, a prime, chosen so that all computations are kept below 32 bit precision 
    global powers 
    powers = tuple(pow(4,i,q) for i in reversed(range(k)))      #precalculate the powers
    bm = pow(4,k,q)
    occurrences = {}     #initialise empty dictionary
    nucleotides = "ATGC"
    S = [toNumeral(c) for c in S if c in nucleotides]   #sanitise S, can removed if S is known to be "clean"
    length=len(S)
    if(length < k):
        return [],[],[]
    print("Starting search")

    currHash = getHash(S[:k],q)
    occurrences[currHash] = [[S[:k],0]]

    for index in range(1,length-k+1):
        if(index%300000 == 0):          #delete this later
            print("------------------------------\nProgress: ",index)
            print("Size of dictionary: ", sys.getsizeof(occurrences))
            print("Size of random entry: ",realSize(occurrences[currHash]))
            print("Number of collisions in random entry: ", len(occurrences[currHash]))
            print("Number of entries: ",len(occurrences))
        currHash = (4*currHash- S[index-1]*bm + S[index-1+k]) % q   #update hash
        kmer = S[index : index+k]
        
        position = findKmer(occurrences, kmer, currHash)    #find element

        if position == -2:          #add element to dictionary
            occurrences[currHash] = [[kmer,index]]
        elif position == -1:
            occurrences[currHash].append([kmer,index])
        else:
            occurrences[currHash][position].append(index)

    L1,L2,L3 = [],[],[]
    for results in occurrences.values():    #process data, see ksearch for more details
        for kmerInfo in results:
            if(len(kmerInfo) > f):
                kmer = kmerInfo[0]
                locations = kmerInfo[1:]
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
                                L3[-1] += len(occurrences[currHash][i])-1
                                break
    return L1,L2,L3

def ksearch(S,k,f,x,debug=False,sanitise=False):
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

    Discussion: 
    The algorithm first sanitises the input string (if a nucleotide is unknown this may be written by using a letter other than "ATGC"). The data is stored in a dictionary with the following structure: {"kmer": [list of locations of kmer]}. The first kmer is added to the dictionary and then the code moves along the strand adding every location. Once this is done it starts computing L1,L2,L3. It iterates through the values in the dictionary (discarding those whose frequency is less than the cutoff), first appending the kmer and its locations to L1 and L2, then iterating through all possible point-x mutations and adding their frequency to L3.

    Running time:
    Let N be the length of the DNA strand. Sanitising the string is an O(N) operation, but it is still relatively quick (it can also be removed if the string S is known to be already sanitised), so I'm going to ignore this factor. The main loop executes N-k times and (ignoring the debugging print statements which only get executed every 300'000 iterations) each loop has a running time of O(k): accesses to the dictionary are constant time as the implementation of hash that python uses is constant as long as the string length isn't too long, ie below a few hundred characters; this means that the only non constant term is the allocation of kmer which is an O(k) process, even if the allocation isn't spelled out as in the code below Python will have to take an equivalent step, so there is no way to avoid this. Finally we have an array allocation or creation at each step. So the O time estimate for the main loop is O((N-k)k)=O(Nk) as k<<N.

    The second loop iterates about min(4^k, N-k). If the frequency is above the desired threshold, each loop takes a further O(k) iterations as it has to manipulate k long strings. So in the worst case scenario (f=1, each kmer only appears once) it will take O((N-k)k) time. However, it is far more likely to need much, much less, as it is unlikely that f=1 is ever an interesting case, and for k large we will have few matches (for example there ar only 2 sequences longer than 200 in the first 10'000'000bp of the genome of Arabidopsis thaliana (genome downloaded from https://www.ncbi.nlm.nih.gov/genome/?term=Arabidopsis+thaliana)). So in practice the final computation is much faster than the initial one. This means that the expected execution time is O(Nk) which is in line with the best known algorithms (although this implementation has much higher constants and memory usage).

    Explanation of development choices:
    The algorithm was developed with a small k kept in mind: most real world applications are limited to double digit applications (https://www.cbcb.umd.edu/software/jellyfish/ for example was limited to k<=31 until recently). In any case if larger values are needed a better approach would be to first find frequent k'mers where k'<<k and then performing the search through those possible matches. This drastically reduces the memory requirements. Better, more complex, data structures are also needed to hold the data (for example a bloom filter) and this felt out of scope. 
    
    I implemented two versions of Rabin-Karp but in practice they where slower and more memory intensive than this simple dictionary search: the overhead of having to store the kmer string (again, implementing a proper hash table with a good way to handle hash collisions felt out of scope) and the inherit slowness of purely pythonic functions when compared to wrapped C/Fortran means that the dictionary solution worked best. The implementations are in ksearchNum and ksearchStr: the first one converts the string at the very start into a list of numbers 0->3 and then does all the computations, in the second case the string is kept as such with conversion only happenning when needed.
    """
    if(x >= k):
        raise ValueError("x should be less than k, value of x: %d\tValue of k: %d"% (x,k))
    occurrences = {}     #initialise empty dictionary
    if sanitise:
        nucleotides = "ATGC"
        S = "".join([c for c in S if c in nucleotides])
    length=len(S)
    if(length < k):
        return [],[],[]

    print("starting search")
    occurrences[S[:k]] = [0]
    for index in range(1,length-k+1):
        if(index%300000 == 0 and debug):          #delete this later
            print("------------------------------\nProgress: ",index)
            print("Size of dictionary: ", sys.getsizeof(occurrences))   #only includes size of overhead, not referenced objects
            print("Size of random entry: ",realSize(occurrences[S[index-1 : index+k-1]]))
            print("Number of entries: ",len(occurrences))
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
                kmerMutated = "".join(kmerMutated)
                
                if("".join(kmerMutated) in occurrences):
                    L3[-1] += len(occurrences[kmerMutated])
    return L1,L2,L3

if __name__=='__main__':
    #Sample input and function call. Use/modify if/as needed
    S=""
    for seq_record in SeqIO.parse("genomes/sequence.fasta", "fasta"):
        S=seq_record.seq[:]

    S=S[:len(S)//10]    
    k=3
    f=1000
    x=0
    method = ksearch
    if(len(sys.argv) == 5):#usage: python p1.py k f x {d,s,n}
        k = int(sys.argv[1])
        f = int(sys.argv[2])
        x = int(sys.argv[3])
        m = sys.argv[4]
    
        if m == 'd':
            method = ksearch 
        elif m == 's':
            method = ksearchStr
        elif m == 'n':
            method = ksearchNum

    L1,L2,L3=method(S,k,f,x)
