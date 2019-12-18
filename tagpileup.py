import csv 
import matplotlib.pyplot as plt

class GenomeInfo: 
    def __int__ (self):
        self.chromosomes =[]

class Peak:
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


def DecimalToBinary(n): 
    return bin(n).replace("0b", "") 


def NthDigit(number,n): 
    tmp =number // (10**(n-1))
    tmp2=tmp %10
    return (tmp2)


def ReadInputFile(fileName):
    with open(fileName) as data:                                                                                          
        data_reader = csv.reader(data, delimiter='\t')
        d=list (data_reader)
    
    return d

def AssessGenomeSize   (MappedData):

    
    numberBps=[]
    for i in range ( len (MappedData )):
        
        if ( str(MappedData[i][0]) != '@SQ'):
            break

        tmpstring1=str(MappedData[i][1])
        tmpstring2=str(MappedData[i][2])

        if (tmpstring1[3:6]=='chr' and tmpstring1[6:] !='M'):
            chrm = int (tmpstring1[6:])
            numberBp= int (tmpstring2[3:])
            numberBps.append (numberBp)

    return numberBps



def ParseSAMForChipexo (MappedData,numberBps): 

    samInfo = SamInfo() 

    for i in range (len (numberBps)):
        samInfo.numTagsP.append ( [0 for j in range (numberBps[i])])
        samInfo.numTagsN.append ( [0 for j in range (numberBps[i])])

    for i in range ( len (MappedData )):
        tmpstring=str(MappedData[i][2])
        if (tmpstring[0:3]=='chr' and tmpstring != 'chrM'):
            chrm = int (tmpstring[3:])
            basepair = int (MappedData[i][3])
            bitwiseFlag = int (MappedData[i][1])
            bitwiseFlagBinary= int ( DecimalToBinary (bitwiseFlag) )
            if ( NthDigit (bitwiseFlagBinary,7)==1 and NthDigit (bitwiseFlagBinary,3)==0):     # read is first pair & read is not unmapped  
                if (NthDigit (bitwiseFlagBinary,5)==1 ): # bitwiseflag for negative strand
                    samInfo.numTagsN [chrm-1][basepair-1]=samInfo.numTagsN [chrm-1][basepair-1] +1
                else:
                    samInfo.numTagsP [chrm-1][basepair-1]=samInfo.numTagsP [chrm-1][basepair-1] +1
    
    return samInfo


def ParsePeakFile (motifWindow):
    peaks=Peak()
    for i in range ( len (motifWindow )):
        tmpstring=str(motifWindow[i][0])
        if (tmpstring[0:3]=='chr' and tmpstring != 'chrM'):
            peaks.idChr.append     ( int (tmpstring[3:])   )
            peaks.firstBp.append   ( int (motifWindow[i][1]) )
            peaks.lastBp.append    ( int (motifWindow[i][2]) )
            peaks.centerBp.append  ( int (0.5 * ( int (motifWindow[i][1]) + int (motifWindow[i][2]) ) )) 
            peaks.strandDir.append ( str(motifWindow[i][5])  )
    
    return (peaks)



def main():
    #### inputs
    windowSize=502
    
    motifsWindow  = ReadInputFile('MotifReb1.bed')
    peaks         = ParsePeakFile (motifsWindow )
    mappedData    = ReadInputFile('result.sam')
    numberBps     = AssessGenomeSize( mappedData)
    samInfo       = ParseSAMForChipexo (mappedData,numberBps)
    
    tagPileUpP=[ 0 for x in range (windowSize)]
    tagPileUpN=[ 0 for x in range (windowSize)]
    for i in range (len (peaks.centerBp)):
        chrom=peaks.idChr[i]
        ref=peaks.centerBp[i]
        for j in range (numberBps[chrom-1]):
            if ( j >= peaks.firstBp[i] and j <= peaks.lastBp[i] and peaks.strandDir[i]=='+'):
                tagPileUpP[j - ref + int(windowSize/2) - 1] = tagPileUpP[j - ref + int(windowSize/2) - 1] +  samInfo.numTagsP[chrom-1][j]
            if  ( j >= peaks.firstBp[i] and j <= peaks.lastBp[i] and peaks.strandDir[i]=='-'):
                tagPileUpN[j - ref + int(windowSize/2) - 1] = tagPileUpN[j - ref + int(windowSize/2) - 1] +  samInfo.numTagsN[chrom-1][j]

    plt.plot ( range(len (tagPileUpN)),tagPileUpN)  
    plt.plot ( range(len (tagPileUpP)),tagPileUpP) 
    plt.show()

if __name__ == "__main__":
     main()
        
