#! /bin/sh


# make directories for the results
if [ -d ../result ];
then
    mkdir -p ./result
fi


touch  ../result/RK.mut2trans.txt
while read sample; 
do 
    # filter splicing junctinos with >= 5 support reads
    echo "python simpleFilt.py ../data/mapsplice/${sample}C.junctions.txt.gz 5 > ../result/${sample}.junctions.filt.txt"
    # python simpleFilt.py ../data/mapsplice/${sample}C.junctions.txt.gz 5 > ../result/${sample}.junctions.filt.txt

    # annotate splicing junctions (exon skip, splicing-site slip, pseudo exon inclusion and so on) 
    echo "python annotSplicing.py ../result/${sample}.junctions.filt.txt ../data/db/refGene.bed.gz ../data/db/refExon.bed.gz > ../result/${sample}.annot.txt"
    # python annotSplicing.py ../result/${sample}.junctions.filt.txt ../data/db/refGene.bed.gz ../data/db/refExon.bed.gz > ../result/${sample}.annot.txt

    # extract putative mutations causing splicing changes
    echo "python getSplicingMut.py ../result/${sample}.annot.txt ../data/vcf/${sample}.mutation.vcf.gz ../data/db/refExon.bed.gz > ../result/${sample}.mut2sp.txt"
    # python getSplicingMut.py ../result/${sample}.annot.txt ../data/vcf/${sample}.mutation.vcf.gz ../data/db/refExon.bed.gz > ../result/${sample}.mut2sp.txt

    echo "python mergeMut.py ../data/vcf/${sample}.mutation.vcf.gz ../result/${sample}.mut2sp.txt ../data/db/refExon.bed.gz ${sample} >> ../result/RK.mut2trans.txt"
    python mergeMut.py ../data/vcf/${sample}.mutation.vcf.gz ../result/${sample}.mut2sp.txt ../data/db/refExon.bed.gz ${sample} >> ../result/RK.mut2trans.txt

done < ../data/sample_list.txt 

