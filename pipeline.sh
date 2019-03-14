CTHOME=/data/groups/heewlee/MA/wazim/TimeSeries/ClonalTREE
ref=/data/groups/heewlee/MA/wazim/TimeSeries/FINAL/util/ecoli.v3/K12MG1655.fna
glist=/data/groups/heewlee/MA/wazim/TimeSeries/FINAL/util/ecoli.v3/K12MG1655_gList.txt
data=/data/groups/heewlee/MA/wazim/TimeSeries/test_pipeline
prefix=pop125
TRIMHOME=/data/groups/heewlee/MA/bin/Trimmomatic-0.33/
TRIMJAR=/data/groups/heewlee/MA/bin/Trimmomatic-0.33/trimmomatic-0.33.jar

declare -a arr1=("t1" "t2" "t3" "t4" "t5" "t6")
mkdir -p ${data}/trimmed

if [[ "$1" == *1* ]]
then
    for i in "${arr1[@]}"
    do
        java -jar ${TRIMJAR} PE -threads 1 -phred33 -trimlog $data/trimmed/${prefix}_${i}.trimlog $data/${prefix}_${i}_1.fq.gz $data/${prefix}_${i}_2.fq.gz $data/trimmed/${prefix}_${i}_1.fq.gz $data/trimmed/${prefix}_${i}_1.fq.singleton.gz $data/trimmed/${prefix}_${i}_2.fq.gz $data/trimmed/${prefix}_${i}_2.fq.singleton.gz ILLUMINACLIP:${TRIMHOME}adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:70 &> $data/trimmed/${prefix}_${i}.trimstat
    done
fi

if [[ "$1" == *2* ]]
then
    for i in "${arr1[@]}"
    do
        bwa mem -M -t 2 $ref $data/trimmed/${prefix}_${i}_1.fq.gz $data/trimmed/${prefix}_${i}_2.fq.gz > ${data}/${prefix}_PE
        samtools view -Sb -t $glist ${data}/${prefix}_PE -o ${data}/${prefix}_tmp.bam
        samtools sort ${data}/${prefix}_tmp.bam ${data}/${prefix}_tmp.sorted
        samtools mpileup -f $ref ${data}/${prefix}_tmp.sorted.bam > ${data}/${prefix}_${i}.pileup
        rm -f ${data}/${prefix}_tmp.bam ${data}/${prefix}_tmp.sorted.bam ${data}/${prefix}*PE
    done
fi

if [[ "$1" == *3* ]]
then
    for i in "${arr1[@]}"
    do
        perl ${CTHOME}/printBases.pl ${data}/${prefix}_${i}.pileup > ${data}/${prefix}_${i}.printBases
        python ${CTHOME}/bialleleFreq.py ${data}/${prefix}_${i}.printBases ${i}
    done
fi

if [[ "$1" == *4* ]]
then
    cat ${data}/${prefix}_*.freqs > ${data}/${prefix}.freqs
    python ${CTHOME}/filterFreqs.py ${data}/${prefix}.freqs ${data}/${prefix}
fi


