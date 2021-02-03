# -*- coding: utf-8 -*-

"""
Assignment:
 make a new version of the naive function called naive_2mm 
 that allows up to 2 mismatches per occurrence. Unlike for
 the previous questions, do not consider the reverse complement 
 here.  We're looking for approximate matches for P itself, not
 its reverse complement.
"""

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

"""
Program below created on Thu Oct 22 22:32:23 2020

@author: Michael
"""

def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments t sequence p read
        match = True
        n=0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                n = n + 1   
        if n <=2:    
            match = True
        elif n >2:
            match = False

        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

#Example 1
p = 'CTGT'
ten_as = 'AAAAAAAAAA'
t = ten_as + 'CTGT' + ten_as + 'CTTT' + ten_as + 'CGGG' + ten_as
occurrences = naive_2mm(p, t)
print(occurrences)

#Example 2
phix_genome = readGenome('phix.fa')
occurrences = naive_2mm('GATTACA', phix_genome)
print('offset of leftmost occurrence: %d' % min(occurrences))
print('# occurrences: %d' % len(occurrences))

#Example 3
lambda_virus = readGenome('lambda_virus.fa')
occurrences = naive_2mm('TTCAAGCC', lambda_virus)
print('offset of leftmost occurrence: %d' % min(occurrences))
print('# occurrences: %d' % len(occurrences))

#Example 4
lambda_virus = readGenome('lambda_virus.fa')
occurrences = naive_2mm('AGGAGGTT', lambda_virus)
print('offset of leftmost occurrence: %d' % min(occurrences))
print('# occurrences: %d' % len(occurrences))
      







