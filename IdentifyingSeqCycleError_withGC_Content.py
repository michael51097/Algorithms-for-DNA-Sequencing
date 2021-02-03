# -*- coding: utf-8 -*-

"""
ASSIGNMENT:
This dataset has something wrong with it; 
one of the sequencing cycles is poor quality.

Report which sequencing cycle has the problem.  
Remember that a sequencing cycle corresponds to 
a particular offset in all the reads. For example,
 if the leftmost read position seems to have a 
 problem consistently across reads, report 0. If
 the fourth position from the left has the problem,
 report 3. Do whatever analysis you think is needed 
 to identify the bad cycle.
"""
"""
readFastq provided
"""
def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()
            seq = fh.readline().rstrip()
            fh.readline()
            qual = fh.readline().rstrip()
            if len(seq) ==0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

"""
Program below created on Mon Feb  1 20:04:01 2021

@author: Michael
"""

import matplotlib.pyplot as plt

#Different mix of GC based off of the position in the genome
def findGCByPos(reads):
    gc = [0]*100
    totals = [0]*100
    
    for read in reads:
        for i in range(len(read)):
            if read[i]== 'C' or read[i] == 'G':
                gc[i]+= 1
            totals[i] +=1
                
    for i in range(len(gc)):
        if totals[i] > 0:
            gc[i] /= float(totals[i])
    return gc

seqs,qual = readFastq('ERR037900_1.first1000.fastq')

gc = findGCByPos(seqs)
plt.plot(range(len(gc)),gc)
plt.show()

min_value = min(gc)
min_index = gc.index(min_value)



print("In order to identify where the poor sequencing cycle is. I calculated the GC content")
print(" accross all the reads and looked for an outlier in GC content ")
print("based off the graph displayed I identified an outlier, in the form of minimum")
print("value of "+ str(min_value) + " at sequencing cycle " + str(min_index)) 




#https://www.coursera.org/learn/dna-sequencing/lecture/YkDXI/practical-analyzing-reads-by-position
