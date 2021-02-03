# -*- coding: utf-8 -*-
#Homework #1 week 1

"""
Assignment:

First, implement a version of the naive exact
 matching algorithm that is strand-aware. That is,
 instead of looking only for occurrences of P in T, 
 additionally look for occurrences of the reverse complement 
 of P in T. If P(the read) is ACT, your function should find 
 ccurrences of both ACT and its reverse complement AGT in T(the sequence).
 
 If P and its reverse complement are identical (e.g. AACGTT),
 then a given match offset should be reported only once. So
 if your new function is called naive_with_rc, then the old 
 naivefunction and your new naive_with_rc function should 
 return the same results when P equals its reverse complement.
 
"""
"""
Functions below provided by course
"""
def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


def naive(p, t): # t is the sequence and p is the read
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences


def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t
"""
Programs below created on Wed Oct 21 15:31:13 2020

@author: Michael
"""

def naive_with_rc(p,t): 
    #Convert the pattern to its reverse complement
    occurrences = []
    pr = reverseComplement(p) 
    
    searches = [p,pr]
    
    if p != pr:
        print("made it past if")
        #print(p)
        #print(pr)
        
        for search in searches:
            print(search)
            for i in range(len(t) - len(search) + 1):  # loop over alignments t sequence p read
                match = True
                for j in range(len(search)):  # loop over characters
                    
                    if t[i+j] != search[j]:  # compare characters
                        match = False
                        break
                if match:
                    occurrences.append(i)  # all chars matched; record
                    
    else:
        for i in range(len(t) - len(p) + 1):  # loop over alignments
            match = True
            for j in range(len(p)):  # loop over characters
                if t[i+j] != p[j]:  # compare characters
                    match = False
                    break
            if match:
                occurrences.append(i)  # all chars matched; record
        return occurrences
                
    return occurrences
    
#Example 1
p = 'CCC'
ten_as = 'AAAAAAAAAA'
t = ten_as + 'CCC' + ten_as + 'GGG' + ten_as
print(t)
occurrences = naive_with_rc(p, t)
print(occurrences)

#Example 2
p = 'CGCG'
t = ten_as + 'CGCG' + ten_as + 'CGCG' + ten_as
occurrences = naive_with_rc(p, t)
print(occurrences)

#Example 3
DNA = readGenome('phix.fa')
occurrences = naive_with_rc('ATTA', DNA)

print('offset of leftmost occurrence: %d' % min(occurrences))
print('# occurrences: %d' % len(occurrences))

#Example 4
DNA = readGenome('lambda_virus.fa')
occurrences = naive_with_rc('AGGT', DNA)

print('offset of leftmost occurrence: %d' % min(occurrences))
print('# occurrences: %d' % len(occurrences))

#Example 5
DNA = readGenome('lambda_virus.fa')
occurrences = naive_with_rc('TTAA', DNA)

print('offset of leftmost occurrence: %d' % min(occurrences))
print('# occurrences: %d' % len(occurrences))

#Example 6 
DNA = readGenome('lambda_virus.fa')
occurrences = naive_with_rc('ACTAAGT', DNA)

print('offset of leftmost occurrence: %d' % min(occurrences))
print('# occurrences: %d' % len(occurrences))

#Example 7 
DNA = readGenome('lambda_virus.fa')
occurrences = naive_with_rc('AGTCGA', DNA)

print('offset of leftmost occurrence: %d' % min(occurrences))
print('# occurrences: %d' % len(occurrences))



