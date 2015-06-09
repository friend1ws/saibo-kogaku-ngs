#! /usr/local/bin/python

import sys, re, tabix, gzip

mutFile = sys.argv[1]
mut2spFile = sys.argv[2]
exonFile = sys.argv[3]
sample = sys.argv[4]

geneRe = re.compile('(\w+)\(N[MR]_\d+\)')
hIN = open(mut2spFile, 'r')

mut2sp_exist = []
for line in hIN:
    F = line.rstrip('\n').split('\t')

    candGenes = []
    genes = F[4].split(";")
    for gene in genes:

        geneMatch = geneRe.match(gene)
        if geneMatch is not None:
            candGenes.append(geneMatch.group(1))

    candGenes = list(set(candGenes))
    mut2sp_exist.append(F[10] + '\t' + F[11])

    for cGene in candGenes:
        print cGene + '\t' + sample + '\t' + "point mutation,indel" + '\t' + "splicing aberration"
   

hIN.close()

gene_exist = []
hIN = gzip.open(mutFile)
exon_tb = tabix.open(exonFile)

for line in hIN:
    F = line.rstrip('\n').split('\t')

    if F[0] + '\t' + F[1] in mut2sp_exist: continue
    
    # check exon and junction annotation for the side 1  
    tabixErrorFlag = 0
    try:
        records = exon_tb.query(F[0], int(F[1]) - 2, int(F[1]) + 2)
    except Exception as inst:
        # print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        # print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag = 1

    tGenes = [];
    if tabixErrorFlag == 0:
        for record in records:

            geneMatch = geneRe.match(record[3])
            if geneMatch is not None:
                tGenes.append(geneMatch.group(1))

    tGenes = list(set(tGenes))
    
    for gene in tGenes:
        if gene not in gene_exist:
            print gene + '\t' + sample + '\t' + "point mutation,indel" + '\t' + "---"
            gene_exist.append(gene)

hIN.close()

# chrX  489336041    48934304    exon-skip   WDR45(NM_001029896);WDR45(NM_007075)    6;7 s;s WDR45(NM_001029896);WDR45(NM_007075)    4;5 e;e chrX    48934184    .   C   T   60  PASS    SOMATIC chrX:48934183-48934191,-    splicing acceptor disruption
    
