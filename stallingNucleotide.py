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
for row in reader:
    print(row)
    tss=row[2]
    print('tss:'+tss)
    geneStart=row[5]
    print('geneStart:'+geneStart)
    geneID=row[4]
    print('gengID:'+geneID)
    geneIDs=geneID.split('//')
    print(geneIDs)
    geneStarts=geneStart.split('//')
    print(geneStarts)
    direction=row[7]
    for i in range(0,len(geneStarts)):
        if tss:
            for seq_record in SeqIO.parse("sequence.fasta", "fasta"):

                if direction is '-':
                    print(tss)
                    test_feature=SeqFeature(FeatureLocation(int(geneStarts[i]),int(tss)),type="gene",strand=-1)
                else:
                    test_feature=SeqFeature(FeatureLocation(int(tss),int(geneStarts[i])),type="gene",strand=1)

                #test_feature=SeqFeature(FeatureLocation(40417,TSS),type="gene",strand=strandDirection)
                example_seq=test_feature.extract(seq_record)

                print(example_seq.seq)

                sequence=str(example_seq.seq.transcribe())
                print sequence
                letters = ['U','A','C','G']
                stallingPosition=['x','x','x','x']
                for letter in letters:
                    stallPosition=sequence.index(letter) if letter in sequence else ''
                    
                    stallingPosition[letters.index(letter)]=str(stallPosition)
                    print('Lack of '+letter+' stalls gene '+geneIDs[i] +' at '+str(stallPosition))

                print(stallingPosition)
                if not all(pos <= 0 for pos in stallingPosition):
                    outputFile.write(geneIDs[i].lstrip()+','+stallingPosition[0]+','+stallingPosition[1]+','+stallingPosition[2]+','+stallingPosition[3]+'\n')
                
    #seq_record_seq=seq_record.seq
    #print(seq_record_seq[42037:42137])
outputFile.close()
