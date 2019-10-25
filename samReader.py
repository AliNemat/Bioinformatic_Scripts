# This script is written to examine if all the mapped reads which are different in two SAM files belongs to the reads
# that have altenrative option for mapping. this question arise when we were notcing that BAW MEM gives different SAM files for the same input FASTQ files.

# This code gets two input file. 
#  -The first input is the output of BamUtil (https://genome.sph.umich.edu/wiki/BamUtil) 
# $ ../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --out results/diffOrderSam.log
#  -The second input file is output of BWA-MEM. The code then extract the SAMs mapped reads which have the flag of "XA:", means alternative mapped to referebce genome does exist.
# $ bwa mem ref.fa read1.fq read2.fq > aln.sam
# The code output is the reads which are overlap and different between these two files.

import csv 

def ReadInputFile(fileName):
    with open(fileName) as data:                                                                                          
        data_reader = csv.reader(data, delimiter='\t')
        d=list (data_reader)
    
    
    return d

def main():
    #### inputs

    ## finding unique reads which is shwon by BamUtil software
    tmpDiffData =ReadInputFile('diffOrderSam4.log')
    tmpMapData =ReadInputFile('result.sam')


    ## finding unique reads which is shwon by BamUtil software
    diffData=[]
    for i in range (len (tmpDiffData)):
        if (len (tmpDiffData[i]) ==1) :
            diffData.append (str(tmpDiffData[i][0]))
    
    diffDataUnique=set (diffData)
    diffDataUnique=list (diffDataUnique)


    multipleMap2D=[]
    multipleMap1D=[]
    for i in range (20, len (tmpMapData)):
        tmpstring=str(tmpMapData[i][-1])
        #bitwiseFlag=int (tmpMapData[i][1])
        if (tmpstring[:2]=='XA'):
        #if (bitwiseFlag>=260):
            #multipleMap2D.append ([ str(tmpMapData[i][0]),str(tmpMapData[i][-1]) ] )
            multipleMap1D.append (  str(tmpMapData[i][0]) )

    multipleMap1DUnique=set (multipleMap1D)
    multipleMap1DUnique=list(multipleMap1DUnique)
    

    overlap=[]
    different=[]
    for i in range (len(diffDataUnique)): 
        if (diffDataUnique[i] in multipleMap1DUnique):
            overlap.append (diffDataUnique[i])
        else:
            different.append(diffDataUnique[i]) 

    print "Tasks finished" 
    f = open("overlap.txt", "w")
    g = open("different.txt", "w")

    for i in range (len(overlap)):
        f.write(overlap[i])
        f.write("\n")

    for i in range (len(different)):
        g.write(different[i])
        g.write("\n")
    



if __name__ == "__main__":
     main()
        
