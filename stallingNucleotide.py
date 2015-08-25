from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature,FeatureLocation
import csv

promoterFile=open('AllEcoCycSigma70.csv','rU')
reader=csv.reader(promoterFile)
outputFile=open('StallingPositionEcoCycSigma70.csv','w')
next(reader)
outputFile.write('Gene,'+'Lack Of U,'+'Lack Of A,'+'Lack Of C,'+'Lack Of G \n')

#for each row in input file from BioCyc.org SmartTable
for row in reader:
    print(row)
    #Transcrition start site is in column 2
    tss=row[2]
    print('tss:'+tss)

    #gene start is on column 5
    geneStart=row[5]
    print('geneStart:'+geneStart)

    #gene ID is in column 4
    geneID=row[4]
    print('gengID:'+geneID)

    #can be multiple gene IDs separated by ''//''
    geneIDs=geneID.split('//')
    print(geneIDs)

    #in which case the multiple gene starts will also be separated by a //
    geneStarts=geneStart.split('//')
    print(geneStarts)

    #direction of transcription is in column 7 either + or -
    direction=row[7]

    #for every gene start position
    for i in range(0,len(geneStarts)):
        #if we had a TSS (not ever gene does)
        if tss:
            #open the fasta file sequence.fasta 
            for seq_record in SeqIO.parse("sequence.fasta", "fasta"):
                #check from genestart to TSS is forward or tss to gene start if backwards with -1 or 1 strange appropraitely
                if direction is '-':
                    #print(tss)

                    test_feature=SeqFeature(FeatureLocation(int(geneStarts[i]),int(tss)),type="gene",strand=-1)
                else:
                    test_feature=SeqFeature(FeatureLocation(int(tss),int(geneStarts[i])),type="gene",strand=1)


                #test_feature=SeqFeature(FeatureLocation(40417,TSS),type="gene",strand=strandDirection)
                
                #example_seq now contains the sequence which we need to search through
                example_seq=test_feature.extract(seq_record)

                print(example_seq.seq)

                #need to transcribe the sequence
                sequence=str(example_seq.seq.transcribe())
                print sequence

                letters = ['U','A','C','G']
                
                #default the stalling position to x
                stallingPosition=['x','x','x','x']

                #for each letter in UACG
                for letter in letters:

                    #set stallPosition to the index where that letter occures in the transcribed sequence
                    stallPosition=sequence.index(letter) if letter in sequence else ''
                    
                    #store the stall position of each nucleotide in the stallingPosition array [UACG]
                    stallingPosition[letters.index(letter)]=str(stallPosition)
                    print('Lack of '+letter+' stalls gene '+geneIDs[i] +' at '+str(stallPosition))

                print(stallingPosition)

                #if all of the indexes were stored write it to the otput file
                if not all(pos <= 0 for pos in stallingPosition):
                    outputFile.write(geneIDs[i].lstrip()+','+stallingPosition[0]+','+stallingPosition[1]+','+stallingPosition[2]+','+stallingPosition[3]+'\n')
                
    #seq_record_seq=seq_record.seq
    #print(seq_record_seq[42037:42137])
outputFile.close()
