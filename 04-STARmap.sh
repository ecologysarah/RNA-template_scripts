#!/bin/bash

#This script will map reads to a reference genome using STAR. It first checks if the reference has been indexed and, if not, indexes it. It relies on correctly setting the variables indicated below.
#Written by Sarah Christofides, 2023. Released under Creative Commons BY-SA.

###VARIABLES TO BE SET###
#Set the path to your directory on scratch
myDir=
#Set your username
userProject=
#Set the path to your reference genome directory, FASTA file and GTF/GFF file
GENOME=${myDir}/resources/humanGRCh38
FASTA=GRCh38_latest_genomic.fna
GTF=GRCh38_latest_genomic.gff
#Indicate if the data is single-end (SE) or paired-end (PE)
ends=SE
#Set the slurm queue to use: defq for gomphus, epyc for iago, htc for hawk
queue=epyc
######

mem="40G"
cpu="4"
runTime="02:00:00"
scriptBase="04map"

#Set the correct version of STAR
STAR=$(module avail -L star | tail -n 1)
##Append this information to the report
echo -e "\nReads mapped to genome with ${STAR}" >> ${myDir}/AnalysisReport.txt


#Make the output directory
if [ ! -d "${myDir}"/04-star ]
        then mkdir ${myDir}/04-star
fi

#Step1: Index the reference genome

scriptName=${myDir}/temp/${scriptBase}.index.sh
touch ${scriptName}

echo "#!/bin/bash" > ${scriptName}

echo "#SBATCH --partition=${queue}" >> ${scriptName} 
echo "#SBATCH --mem-per-cpu=40G" >> ${scriptName} 
echo "#SBATCH --nodes=1" >> ${scriptName}
echo "#SBATCH --tasks-per-node=${cpu}" >> ${scriptName}
echo "#SBATCH -o ${myDir}/OUT/${scriptBase}${jobName}.%J" >> ${scriptName}
echo "#SBATCH -e ${myDir}/ERR/${scriptBase}${jobName}.%J" >> ${scriptName}
echo "#SBATCH --job-name=index" >> ${scriptName}
echo "#SBATCH --account=${userProject}" >> ${scriptName}

echo "if [ ! -f ${GENOME}/SAindex ]; then

	module load ${STAR}

	STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ${GENOME} --genomeFastaFiles ${GENOME}/${FASTA} --limitGenomeGenerateRAM 320000000000 --sjdbGTFfile ${GENOME}/${GTF} --sjdbGTFtagExonParentTranscript Parent

fi

exit 0" >> ${scriptName}

chmod u+x ${scriptName}

INDEXJOB=$(sbatch --parsable ${scriptName})

#Step 2: Map the samples against the genome

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
	touch ${scriptName}

	echo "#!/bin/bash" > ${scriptName} 
        echo "#SBATCH --partition=${queue}" >> ${scriptName} 
        echo "#SBATCH --mem=${mem}" >> ${scriptName}
        echo "#SBATCH --nodes=1" >> ${scriptName}
        echo "#SBATCH --tasks-per-node=${cpu}" >> ${scriptName}
        echo "#SBATCH -t ${runTime}" >> ${scriptName}
        echo "#SBATCH --array=1-${nSamps}" >> ${scriptName}
        echo "#SBATCH --output ${myDir}/OUT/${scriptBase}${jobName}.%A-%a" >> ${scriptName}
        echo "#SBATCH --error ${myDir}/ERR/${scriptBase}${jobName}.%A-%a" >> ${scriptName}
        echo "#SBATCH --account=${userProject}" >> ${scriptName}

        echo "module load ${STAR}" >> ${scriptName}

	## specify the file to work on
        echo "sampleIter=\$(( \${SLURM_ARRAY_TASK_ID} + ${loopSetup[$index]} ))
        sampleID=\$(cat ${myDir}/01-download/SampleFileNames.txt | sed -n \${sampleIter}p)
        echo \${sampleID}" >> ${scriptName}

	# run the star mapping command
	echo -n "STAR \
	--outSAMunmapped Within KeepPairs \
	--outMultimapperOrder Random \
	--outSAMmultNmax 1 \
	--runThreadN ${cpu} \
	--runMode alignReads \
	--quantMode GeneCounts \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix ${myDir}/04-star/\${sampleID}onemap_ \
	--readFilesCommand zcat \
	--genomeDir ${GENOME} \
	--readFilesIn ${myDir}/02-trim/trim_\${sampleID}_R1.f*q.gz" >> ${scriptName} 
	if [ "${ends}" = PE ]; then echo " ${myDir}/02-trim/trim_\${sampleID}_R2.f*q.gz" >> ${scriptName}; fi
	
	echo -e "\nexit 0" >> ${scriptName}

	chmod u+x ${scriptName}

	sbatch -d afterok:${INDEXJOB} ${scriptName}

	index=$(( $index + 1))

done

exit 0
