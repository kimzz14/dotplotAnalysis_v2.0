fin = open('query.fa')

ref_DICT = {}
seqName_LIST = []
for line in fin:
    if line.startswith('>') == True:
        seqName = line.rstrip('\n').split('\t')[0].split(' ')[0][1:]
        seqName_LIST += [seqName]
        ref_DICT[seqName] = []
    else:
        sequence = line.rstrip('\n')
        ref_DICT[seqName] += [sequence]
fin.close()

for seqName in seqName_LIST:
    ref_DICT[seqName] = ''.join(ref_DICT[seqName])

readLength = 100
skipN = 100
quality = 'I'*readLength
fout = open('query_100.fastq', 'w')

for seqName in seqName_LIST:
    #if seqName != 'ptg000001l': continue
    sequence = ref_DICT[seqName]
    count = -1
    while True:
        count += 1
        sPos = count*skipN
        ePos = sPos + readLength
        if len(sequence) <= ePos: break
        readSeq = sequence[sPos:ePos]
        fout.write('@' + seqName + '_' + str(sPos + 1) + '\n')
        fout.write(readSeq + '\n')
        fout.write('+' + '\n')
        fout.write(quality + '\n')
fout.close()
