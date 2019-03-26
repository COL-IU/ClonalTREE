CTHOME=/data2/wazimoha/wastewater/ClonalTREE
ref=/data2/wazimoha/wastewater/ref/pan.fa
glist=
data=/data2/wazimoha/wastewater/data
prefix=wwmg
TRIMHOME=/data/groups/heewlee/MA/bin/Trimmomatic-0.33/
TRIMJAR=/data/groups/heewlee/MA/bin/Trimmomatic-0.33/trimmomatic-0.33.jar

# declare -a arr1=("t6" "t7")
# declare -a arr1=("t1" "t2" "t3" "t4" "t5" "t6")
# declare -a arr1=("t8" "t9" "t10" "t11" "t12" "t13" "t14" "t15" "t16" "t17" "t18" "t19" "t20" "t21" "t22" "t23" "t24" "t25" "t26" "t27" "t28" "t29" "t30" "t31" "t32" "t33" "t34" "t35" "t36" "t37" "t38" "t39" "t40" "t41" "t42" "t43" "t44" "t45" "t46" "t47" "t48" "t49" "t50" "t51" "t52")
declare -a arr1=("t1" "t2" "t3" "t4" "t5" "t6" "t7" "t8" "t9" "t10" "t11" "t12" "t13" "t14" "t15" "t16" "t17" "t18" "t19" "t20" "t21" "t22" "t23" "t24" "t25" "t26" "t27" "t28" "t29" "t30"   "t31" "t32" "t33" "t34" "t35" "t36" "t37" "t38" "t39" "t40" "t41" "t42" "t43" "t44" "t45" "t46" "t47" "t48" "t49" "t50" "t51" "t52")

if [[ "$1" == *1* ]]
then
    mkdir -p ${data}/trimmed
    for i in "${arr1[@]}"
    do
        java -jar ${TRIMJAR} PE -threads 1 -phred33 -trimlog $data/trimmed/${prefix}_${i}.trimlog $data/${prefix}_${i}_1.fq.gz $data/${prefix}_${i}_2.fq.gz $data/trimmed/${prefix}_${i}_1.fq.gz $data/trimmed/${prefix}_${i}_1.fq.singleton.gz $data/trimmed/${prefix}_${i}_2.fq.gz $data/trimmed/${prefix}_${i}_2.fq.singleton.gz ILLUMINACLIP:${TRIMHOME}adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:70 &> $data/trimmed/${prefix}_${i}.trimstat
    done
fi

if [[ "$1" == *2* ]]
then
    for i in "${arr1[@]}"
    do
        bwa mem -M -t 4 $ref $data/trimmed/${prefix}_${i}_1.fq.gz $data/trimmed/${prefix}_${i}_2.fq.gz > ${data}/${prefix}_PE
        samtools view -Sb ${data}/${prefix}_PE -o ${data}/${prefix}_tmp.bam
        samtools sort ${data}/${prefix}_tmp.bam ${data}/${prefix}_tmp.sorted
        samtools mpileup -f $ref ${data}/${prefix}_tmp.sorted.bam > ${data}/${prefix}_${i}.pileup
        rm -f ${data}/${prefix}_tmp.bam ${data}/${prefix}_tmp.sorted.bam ${data}/${prefix}*PE
    done
fi

if [[ "$1" == *3* ]]
then
    for i in "${arr1[@]}"
    do
        awk '$4 > 10' ${data}/${prefix}_${i}.pileup > ${data}/${prefix}_${i}.pileup.1
        rm -f ${data}/${prefix}_${i}.pileup
        perl ${CTHOME}/printBases.pl ${data}/${prefix}_${i}.pileup.1 > ${data}/${prefix}_${i}.printBases
        python ${CTHOME}/bialleleFreq.py ${data}/${prefix}_${i}.printBases ${i}
    done
fi

if [[ "$1" == *4* ]]
then
    cat ${data}/${prefix}_*.freqs > ${data}/${prefix}.freqs
    python ${CTHOME}/filterFreqs.py ${data}/${prefix}.freqs ${data}/${prefix}
fi


