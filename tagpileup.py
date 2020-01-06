from __future__ import division
import csv 
import matplotlib.pyplot as plt
import math
import numpy as np

class Motif_Loc:
    def __init__ (self): 
        self.idChr=[]
        self.firstBp=[]
        self.lastBp=[]
        self.centerBp=[]
        self.strandDir=[]

class SamInfo:
    def __init__ (self): 
        self.numTagsP=[]
        self.numTagsN=[]

class TagsPileup:
    def __init__ (self):
        self.N=[]
        self.P=[]
        self.Dist=[]
        self.prob=[]
        self.cumulative_prob=[]
        self.avgDistN=None
        self.avgDistP=None
        self.stdDistN=None
        self.stdDistP=None
        self.sumN=None
        self.sumP=None




def DecimalToBinary(n): 
    return bin(n).replace("0b", "") 


def NthDigit(number,n): 
    return ((number // (10**(n-1))) % 10)


def ReadInputFile(fileName):
    with open(fileName) as data:                                                                                          
        data_reader = csv.reader(data, delimiter='\t')
        d=list (data_reader)
    
    return d


def AssessGenomeSize(MappedData):
    numberBps=[]
    for i in range (len(MappedData)):
        tmpString1=str(MappedData[i][0])
        if (tmpString1[0:1] != '@'):
            break

        tmpString2=str(MappedData[i][1])
        tmpString3=str(MappedData[i][2])
        if (tmpString2[3:6] == 'chr' and tmpString2[3:] != 'chrM'):
            numberBps.append (int(tmpString3[3:]))
    if i == 0:
         raise Exception ( 'probably your SAM file does not have a header. Please include a SAM file with a header')
    return numberBps

## to do: does SAM file starts with zero or one? Now it is assumed it starts with zero
def ParseSAMForChipexo (MappedData,numberBps,first_read_len): 
    samInfo = SamInfo() 
    for i in range (len(numberBps)):
        samInfo.numTagsP.append ([0 for j in range (numberBps[i] + first_read_len )])
        samInfo.numTagsN.append ([0 for j in range (numberBps[i] + first_read_len )])

    for i in range (len(MappedData)):
        tmpstring=str(MappedData[i][2])
        if (tmpstring[0:3] == 'chr' and tmpstring != 'chrM'):
            chrm               = int (tmpstring[3:])
            left_most_basepair = int (MappedData[i][3])
            mapq               = int (MappedData[i][4])
            bitwiseFlag        = int (MappedData[i][1])
            bitwiseFlagBinary  = int ( DecimalToBinary (bitwiseFlag) )
            if (NthDigit(bitwiseFlagBinary,7) == 1 and NthDigit(bitwiseFlagBinary,3) == 0 and mapq >= 20):     # read is first pair & read is not unmapped & mapping is accurate with 99% probability  
                if (NthDigit(bitwiseFlagBinary,5) == 1): # bitwiseflag for negative strand
                    samInfo.numTagsN [chrm - 1][left_most_basepair + first_read_len] += 1 ## first_read_len is added to left_most_basepair to reach the 5' end of negative strand
                else:
                    samInfo.numTagsP [chrm - 1][left_most_basepair] += 1                 ## left_most_basepair is already at the 5' end of positive strand
      
    
    return samInfo


def ParseGffFile (motifWindow):
    motif_locs = Motif_Loc()
    for i in range (len(motifWindow)):
        tmpstring=str(motifWindow[i][0])
        if (tmpstring[0:3] == 'chr' and tmpstring != 'chrM'):
            motif_locs.idChr.append     ( int (tmpstring[3:])   )
            motif_locs.firstBp.append   ( int (motifWindow[i][3]) )
            motif_locs.lastBp.append    ( int (motifWindow[i][4]) )
            motif_locs.centerBp.append  ( int (0.5 * ( int (motifWindow[i][3]) + int (motifWindow[i][4]) ) )) 
            motif_locs.strandDir.append ( str(motifWindow[i][6])  )
    return (motif_locs)

    
def CountPileupTags(samInfo,motif_locs,numberBps, expand_size):
    tagsPileup = TagsPileup()

    tagsPileup.P = [0 for x in range(-expand_size, expand_size + 1)]  # The range is [-250, 250] in total 501
    tagsPileup.N = [0 for x in range(-expand_size, expand_size + 1)]
    
    for i in range (len(motif_locs.centerBp)):
        chrom = motif_locs.idChr[i]
        ref   = motif_locs.centerBp[i]
        for j in range (numberBps[chrom - 1]):
            if (j >= (motif_locs.centerBp[i] - expand_size) and j <= (motif_locs.centerBp[i] + expand_size) and motif_locs.strandDir[i] == '+'):
                tagsPileup.P[j - ref + expand_size] += samInfo.numTagsP[chrom - 1][j]
            if (j >= (motif_locs.centerBp[i] - expand_size) and j <= (motif_locs.centerBp[i] + expand_size) and motif_locs.strandDir[i] == '-'):
                tagsPileup.N[j - ref + expand_size] += samInfo.numTagsN[chrom - 1][j]
    
    return tagsPileup

def Generate_tagsPileup_dict (tagsPileup, expand_size):
    tagsPileup_dict=dict()
    
    for i in range (-expand_size, expand_size + 1):    
        tagsPileup_dict [i] = [tagsPileup.P[i + expand_size], tagsPileup.N[i + expand_size]]

    
    print  (tagsPileup_dict)
    return (tagsPileup_dict)
    
def StatsTagsPileup_first(tagsPileup,expand_size):
    tagsPileup.avgDistN = 0.0
    tagsPileup.avgDistP = 0.0
    tagsPileup.stdDistN = 0.0 
    tagsPileup.stdDistP = 0.0
    tagsPileup.sumN = 0.0
    tagsPileup.sumP = 0.0


    for i in range (len (tagsPileup.N)):
        tagsPileup.sumN += tagsPileup.N[i]

    for i in range (len (tagsPileup.P)):
        tagsPileup.sumP += tagsPileup.P[i]
       
    for i in range (len (tagsPileup.N)):
        tagsPileup.avgDistN += abs(i - expand_size) * tagsPileup.N[i] /tagsPileup.sumN
    for i in range (len ( tagsPileup.P)):
        tagsPileup.avgDistP += abs(i - expand_size) * tagsPileup.P[i] /tagsPileup.sumP
    
    tmpSum = 0.0
    for i in range (len ( tagsPileup.N)):
        tmpSum += ((abs(i - expand_size) - tagsPileup.avgDistN) ** 2) * tagsPileup.N[i]
    
    tagsPileup.stdDistN = math.sqrt(tmpSum / (tagsPileup.sumN - 1))
    
    tmpSum = 0.0
    for i in range (len ( tagsPileup.P)):
        tmpSum += ((abs(i - expand_size) - tagsPileup.avgDistP) ** 2) * tagsPileup.P[i] 
        
    tagsPileup.stdDistP = math.sqrt(tmpSum / (tagsPileup.sumP - 1))
    
    print ( 'Standard deviation of tags distance from center of motif for positive strand is %s' %(str(tagsPileup.stdDistP)))
    print ( 'Standard deviation of tags distance from center of motif for negative strand is %s' %(str(tagsPileup.stdDistN)))
    print ( 'Average of tags distance from center of motif for positive strand is %s'            %(str(tagsPileup.avgDistP)))
    print ( 'Average of tags distance from center of motif for negative strand is %s'            %(str(tagsPileup.avgDistN)))


def StatsTagsPileup_second(tagsPileup,expand_size):

    tagsPileup.Dist = [0 for i in range (expand_size + 1)] #1 is added because the center might be different from left or right by one basepair
    for i in range (len (tagsPileup.N)):
        j = int (abs(i - expand_size)) 
        tagsPileup.Dist[j] += tagsPileup.N[i]

    for i in range (len (tagsPileup.P)):
        j = int (abs(i - expand_size)) 
        tagsPileup.Dist[j] += tagsPileup.P[i]
    
    
    tagsPileup.prob = [float ( x / sum(tagsPileup.Dist)) for x in tagsPileup.Dist]
    plt.figure(1)    
    plt.plot (range (len(tagsPileup.prob)), tagsPileup.prob,'r', label = 'Tags probability distribution based on distance from motif')  
    
    tagsPileup.cumulative_prob = [float ( x / sum(tagsPileup.Dist)) for x in np.cumsum(tagsPileup.Dist)]
    plt.figure(2)    
    plt.plot (range (len(tagsPileup.cumulative_prob)), tagsPileup.cumulative_prob, 'b', label = 'Tags cumulative probability distribution based on distance from motif')  

    return (tagsPileup)

def main():
    expand_size     = 250   # expansion size at the right and left side of the center of motif for plotting tag pileup
    first_read_len  = 40  # to do: this should be extracted from SAM file
    gffFile         = ReadInputFile('Reb1_396_24465_Chexmix.gff')
    motif_locs      = ParseGffFile (gffFile)
    samFile         = ReadInputFile('Reb1_396_24465_Chexmix.sam')
    numberBps       = AssessGenomeSize(samFile)
    samInfo         = ParseSAMForChipexo (samFile, numberBps, first_read_len)
    tagsPileup      = CountPileupTags (samInfo, motif_locs, numberBps, expand_size) 
    tagsPileup_dict = Generate_tagsPileup_dict(tagsPileup, expand_size)

    fileM = open('output_tagPileup.txt','w')
    fileM.write ('Distance from center of motif, Number of tags in positive strand, Number of tags in negative strand\n')
    for i in range (-expand_size, expand_size + 1):
        fileM.write(str(i))
        fileM.write('  ')
        fileM.write(str(tagsPileup_dict[i][0]))  
        fileM.write('  ')
        fileM.write(str(tagsPileup_dict[i][1]))  
        fileM.write('  \n')
       
    fileM.close() 
    
    StatsTagsPileup_first (tagsPileup, expand_size)
    StatsTagsPileup_second(tagsPileup, expand_size)
    plt.figure(3)    
    plt.plot (range (-expand_size, expand_size + 1) , tagsPileup.N,'r',label='Negative strand')  
    plt.plot (range (-expand_size, expand_size + 1) , tagsPileup.P,'b',label='Positive strand') 
    plt.show()


if __name__ == "__main__":
    main()
        
