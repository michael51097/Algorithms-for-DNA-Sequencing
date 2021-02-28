# -*- coding: utf-8 -*-
"""
Assignment:
Download and parse the read sequences from the provided Phi-X FASTQ file. 
Next, find all pairs of reads with an exact suffix/prefix match of length 
at least 30. Don't overlap a read with itself; if a read has a suffix/prefix 
match to itself, ignore that match.  Ignore reverse complements

How many distinct pairs of reads overlap?

How many reads have a suffix involved in an overlap?

Functions below provided by course
"""

def overlap(a, b, min_length=4):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
           # print('in start')
           # print(start)
           # print(a[start:])
            return len(a)-start
        start += 1  # move just past previous match
        


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


from itertools import permutations

"Given a set of reads, and overlap length, k, returns tuples of reads that overlap"
def naive_overlap_map(reads, k):
    olaps = {}
    for a,b in permutations(reads,2):
        olen = overlap(a,b, min_length = k)
        if olen > 0:
            olaps[(a,b)] = olen
    return olaps

"""
Code below created on Thu Nov 19 17:28:35 2020

@author: Michael
"""
f = open('ERR266411_1.for_asm.fastq','r')

reads=[]
x=3

"Loads lines of fasta into lines variable"
with open('ERR266411_1.for_asm.fastq') as f:
    lines=f.readlines()

"Extracts the reads from the lines variable"
for line in lines:
    if x%4 == 0:
        reads.append(line[:-1])
    x +=1



def overlap_all_pairs(reads,k):
    kmers = {}
    for read in reads: 
       x = 0
       "Creating a dictionary of kmers for all reads, assigning reads to each kmer"
       while x < (len(read) - (k-1)):
           kmer = read[x:x+k]
           if kmer not in kmers:
               kmers[kmer] = set([read])
           else:
               kmers[kmer].update([read])
           x +=1
           
    master = set()
    unique_members = set()
    
    "Using naive_overlap_map to compare reads that all have the same kmer"
    for reads in kmers: 
        ans=naive_overlap_map(kmers[reads],k)
        if len(ans) > 0:
            "Add the pairs of reads to master set, distinct pairs"
            for pairs in ans:
                master.add(pairs)
                "Add each read to unique_members"
                for each_read in pairs: 
                    unique_members.add(each_read)
    return master, unique_members


overlaps ,read_list = overlap_all_pairs(reads,30)

print("The number of distinct pairs is " + str(len(overlaps)) + " and the number of readswith a suffix involved in an overlap is " + str(len(read_list)) )




