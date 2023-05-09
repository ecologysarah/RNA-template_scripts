#!/bin/bash

#This script will count the number of reads mapping to each gene. It relies on correctly setting the variables indicated below.
#Written by Sarah Christofides, 2022, based on a script from Rob Andrews. Released under Creative Commons BY-SA.

###VARIABLES TO BE SET###
#Set the path to your directory on scratch
myDir=
#Set your username
userProject=
#Set whether you want to use the version of the data with duplicates kept (mark) or removed (rm)
DUP=mark
#Indicate if the data is single-end (SE) or paired-end (PE)
ends=SE
#Set the path to your reference genome directory, FASTA file and GTF/GFF file
GENOME=${myDir}/resources/humanGRCh38
FASTA=GRCh38_latest_genomic.fna
GTF=GRCh38_latest_genomic.gff
#Set the slurm queue to use: defq for gomphus, epyc for iago, htc for hawk
queue=epyc
######

mem="2G"
cpu="2"
scriptBase="07counts"
runTime="00:20:00"

#Set the correct version of SAMtools
SAMTOOLS=$(module avail -L samtools | tail -n 1)
#Set the correct version of subread
SUBREAD=$(module avail -L subread | tail -n 1)
##Append this information to the report
echo -e "\nCounts generated with ${SUBREAD}" >> ${myDir}/AnalysisReport.txt

#Make the output directory
if [ ! -d "${myDir}"/07-featurecounts ]
        then mkdir ${myDir}/07-featurecounts
fi

#Count the reads for each feature in each sample
## Find the number of samples
sampleNo=$(cat ${myDir}/01-download/SampleFileNames.txt | wc -l)

## Set up a slurm array. Each array will process (up to) a thousand samples.
## Create a bash array of numbers incrementing by thousands up to the number of samples. This will be used to determine how many samples are processed in each slurm array, and the correct line for each slurm array to start at in SampleFileNames.txt
loopSetup=( $(echo `seq 0 1000 ${sampleNo}` ${sampleNo}) )

## Loop over each index of loopSetup except the last (as this is used to determine the end point, not to launch a slurm array).
index=0
while [[ $index < $(( ${#loopSetup[@]} -1 )) ]]

do

        ## Calculate the number of samples for the array
        nSamps=$(( ${loopSetup[(($index + 1))]} - ${loopSetup[$index]} ))

        ## write the script to the temp/ directory
        scriptName=${myDir}/temp/${scriptBase}${index}.sh
	rm -rf ${scriptName} || true
        touch ${scriptName}

        echo "#!/bin/bash" >> ${scriptName}
        echo "#SBATCH --partition=${queue}" >> ${scriptName} #epyc
        echo "#SBATCH --mem=${mem}" >> ${scriptName}
        echo "#SBATCH --nodes=1" >> ${scriptName}
        echo "#SBATCH --tasks-per-node=${cpu}" >> ${scriptName}
        echo "#SBATCH -t ${runTime}" >> ${scriptName}
        echo "#SBATCH --array=1-${nSamps}" >> ${scriptName}
        echo "#SBATCH --output ${myDir}/OUT/${scriptBase}${jobName}.%A-%a" >> ${scriptName}
        echo "#SBATCH --error ${myDir}/ERR/${scriptBase}${jobName}.%A-%a" >> ${scriptName}
        echo "#SBATCH --account=${userProject}" >> ${scriptName}

        echo "module load ${SAMTOOLS}" >> ${scriptName}
        echo "module load ${SUBREAD}" >> ${scriptName}

        ## specify the file to work on
        echo "sampleIter=\$(( \${SLURM_ARRAY_TASK_ID} + ${loopSetup[$index]} ))
        sampleID=\$(cat ${myDir}/01-download/SampleFileNames.txt | sed -n \${sampleIter}p)
        echo \${sampleID}" >> ${scriptName}

	## run featurecounts 

	echo "samtools sort -n ${myDir}/05-markduplicates/\${sampleID}_${DUP}dup.bam -o ${myDir}/07-featurecounts/\${sampleID}_${DUP}dup.sorted" >> ${scriptName}

	echo -n "cd ${myDir}/07-featurecounts/
	featureCounts " >> ${scriptName} 
        if [ "${ends}" = PE ]; then echo -n "-p " >> ${scriptName}; fi
	echo "-F GTF -t gene -g ID -a ${GENOME}/${GTF} -o ${myDir}/07-featurecounts/\${sampleID}_${DUP}dup.featurecount ${myDir}/07-featurecounts/\${sampleID}_${DUP}dup.sorted" >> ${scriptName}


        chmod u+x ${scriptName}

        sbatch ${scriptName}

	index=$(( $index + 1))

done

exit 0

