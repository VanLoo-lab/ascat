#!/usr/bin/bash

SNPpos=$1
SIZES=$2
CORES=$3
REF=$4

perl createWindowBed.pl -s ${SIZES} -p ${SNPpos} -c ${CORES}

output=$( basename ${SNPpos} )

for i in $(eval echo "{0..$( expr $CORES - 1 )}")
do
    bedtools nuc -fi ${REF} -bed ${output}"_"${i}".bed" > ${output}"_"${i}".gcContent" &
done

wait

cat *.gcContent > ${output}".combined.gcContent"

R --no-save --args ${output}".combined.gcContent" "GC_"${output}".txt" <createGCcontentFile.R


