import csv 
import matplotlib.pyplot as plt
import math

class Peak:
    def __init__ (self): 
        self.idChr=[]
        self.firstBp=[]
        self.lastBp=[]
        self.centerBp=[]
        self.strandDir=[]
        self.expandSize=None 


class SamInfo:
    def __init__ (self): 
        self.numTagsP=[]
        self.numTagsN=[]

class TagsPileup:
    def __init__ (self):
        self.N=[]
        self.P=[]
        self.Dist=[]
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

    return numberBps


def ParseSAMForChipexo (MappedData,numberBps): 
    samInfo = SamInfo() 
    for i in range (len(numberBps)):
        samInfo.numTagsP.append ( [0 for j in range (numberBps[i])])
        samInfo.numTagsN.append ( [0 for j in range (numberBps[i])])

    for i in range (len(MappedData)):
        tmpstring=str(MappedData[i][2])
        if (tmpstring[0:3] == 'chr' and tmpstring != 'chrM'):
            chrm              = int (tmpstring[3:])
            basepair          = int (MappedData[i][3])
            bitwiseFlag       = int (MappedData[i][1])
            bitwiseFlagBinary = int ( DecimalToBinary (bitwiseFlag) )
            if ( NthDigit(bitwiseFlagBinary,7) == 1 and NthDigit(bitwiseFlagBinary,3) == 0 ):     # read is first pair & read is not unmapped  
                if (NthDigit(bitwiseFlagBinary,5) == 1): # bitwiseflag for negative strand
                    samInfo.numTagsN [chrm-1][basepair-1] += 1
                else:
                    samInfo.numTagsP [chrm-1][basepair-1] += 1
    
    return samInfo


def ParseBedFile (motifWindow):
    peaks=Peak()
    for i in range (len(motifWindow)):
        tmpstring=str(motifWindow[i][0])
        if (tmpstring[0:3] == 'chr' and tmpstring != 'chrM'):
            peaks.idChr.append     ( int (tmpstring[3:])   )
            peaks.firstBp.append   ( int (motifWindow[i][1]) )
            peaks.lastBp.append    ( int (motifWindow[i][2]) )
            peaks.centerBp.append  ( int (0.5 * ( int (motifWindow[i][1]) + int (motifWindow[i][2]) ) )) 
            peaks.strandDir.append ( str(motifWindow[i][5])  )
    # Expansion size around peaks. It is costant for all the peaks in the bed file
    peaks.expandSize=int (peaks.lastBp[0] - peaks.firstBp[0] + 1)  
    return (peaks)

    
def CountPileupTags(samInfo,peaks,numberBps):
    tagsPileup = TagsPileup()

    tagsPileup.P = [0 for x in range(peaks.expandSize)]
    tagsPileup.N = [0 for x in range(peaks.expandSize)]
    
    for i in range (len(peaks.centerBp)):
        chrom = peaks.idChr[i]
        ref   = peaks.centerBp[i]
        for j in range (numberBps[chrom-1]):
            if (j >= peaks.firstBp[i] and j <= peaks.lastBp[i] and peaks.strandDir[i] == '+'):
                tagsPileup.P[j - ref + int(peaks.expandSize / 2) - 1] += samInfo.numTagsP[chrom - 1][j]
            if (j >= peaks.firstBp[i] and j <= peaks.lastBp[i] and peaks.strandDir[i] == '-'):
                tagsPileup.N[j - ref + int(peaks.expandSize / 2) - 1] += samInfo.numTagsN[chrom - 1][j]
    
    return tagsPileup


def StatsTagsPileup_first(tagsPileup,halfExpandSize):
    tagsPileup.avgDistN = 0.0
    tagsPileup.avgDistP = 0.0
    tagsPileup.stdDistN = 0.0 
    tagsPileup.stdDistP = 0.0
    tagsPileup.sumN = 0.0
    tagsPileup.sumP = 0.0


    for i in range (len (tagsPileup.N)):
        tagsPileup.sumN +=tagsPileup.N[i]

    for i in range (len (tagsPileup.P)):
        tagsPileup.sumP +=tagsPileup.P[i]
       
    for i in range (len (tagsPileup.N)):
        tagsPileup.avgDistN += abs(i - halfExpandSize) * tagsPileup.N[i] /tagsPileup.sumN
    for i in range (len ( tagsPileup.P)):
        tagsPileup.avgDistP += abs(i - halfExpandSize) * tagsPileup.P[i] /tagsPileup.sumP
    
    tmpSum = 0.0
    for i in range (len ( tagsPileup.N)):
        tmpSum += ((abs(i - halfExpandSize) - tagsPileup.avgDistN) ** 2) * tagsPileup.N[i]
    
    tagsPileup.stdDistN = math.sqrt(tmpSum / (tagsPileup.sumN - 1))
    
    tmpSum = 0.0
    for i in range (len ( tagsPileup.P)):
        tmpSum += ((abs(i - halfExpandSize) - tagsPileup.avgDistP) ** 2) * tagsPileup.P[i] 
        
    tagsPileup.stdDistP = math.sqrt(tmpSum / (tagsPileup.sumP - 1))
    
    print ( 'Standard deviation of tags distance from center of motif for positive strand is %s' %(str(tagsPileup.stdDistP)))
    print ( 'Standard deviation of tags distance from center of motif for negative strand is %s' %(str(tagsPileup.stdDistN)))
    print ( 'Average of tags distance from center of motif for positive strand is %s'            %(str(tagsPileup.avgDistP)))
    print ( 'Average of tags distance from center of motif for negative strand is %s'            %(str(tagsPileup.avgDistN)))


def StatsTagsPileup_second(tagsPileup,halfExpandSize):

    tagsPileup.Dist=[0 for i in range (halfExpandSize + 1)] #1 is added because the center might be different from left or right by one basepair
    for i in range (len (tagsPileup.N)):
        j = int (abs(i - halfExpandSize)) 
        tagsPileup.Dist[j] += tagsPileup.N[i]

    for i in range (len (tagsPileup.P)):
        j=int (abs(i - halfExpandSize)) 
        tagsPileup.Dist[j] += tagsPileup.P[i]
    
    plt.plot (range (len(tagsPileup.Dist)) , tagsPileup.Dist,'r',label='Tags distribution based on distance from motif')  

    return (tagsPileup)

def main():
    bedFile         = ReadInputFile('MotifReb1.bed')
    peaks           = ParseBedFile (bedFile)
    samFile         = ReadInputFile('out.sam')
    numberBps       = AssessGenomeSize(samFile)
    samInfo         = ParseSAMForChipexo (samFile,numberBps)
    tagsPileup      = CountPileupTags (samInfo,peaks,numberBps) 

    StatsTagsPileup_first(tagsPileup, int(peaks.expandSize / 2))
    StatsTagsPileup_second(tagsPileup, int(peaks.expandSize / 2))

    plt.plot (range (-1*int(peaks.expandSize / 2), int(peaks.expandSize / 2)) , tagsPileup.N,'r',label='Negative strand')  
    plt.plot (range (-1*int(peaks.expandSize / 2), int(peaks.expandSize / 2)) , tagsPileup.P,'b',label='Positive strand') 
    plt.show()


if __name__ == "__main__":
    main()
        
