#!/bin/bash

#This script will identify duplicated reads in the BAM files produced by STAR. It produces two versions of the output: one with duplicates marked and retained, and one with duplicates removed. It relies on correctly setting the variables indicated below.
#Written by Sarah Christofides, 2022, based on a script from Rob Andrews. Released under Creative Commons BY-SA.

###VARIABLES TO BE SET###
#Set the path to your directory on scratch (do not end with a /)
myDir=
#Set your username
userProject=
#Set the slurm queue to use: defq for gomphus, epyc for iago, htc for hawk
queue=epyc
######

mem="20G"
cpu="2"
runTime="00:20:00"
scriptBase="05markdup"

#Set the correct version of picard
PICARD=$(module avail -L picard/ | tail -n 1 | sed -E 's/\s+//')
#Set the correct version of SAMtools
SAMTOOLS=$(module avail -L samtools/ | tail -n 1)
JAVA=$(module avail -L java/ | tail -n 1)
##Append this information to the report
echo -e "\nDuplicates assessed with ${PICARD} and ${SAMTOOLS}" >> ${myDir}/AnalysisReport.txt



#Make the output directory
if [ ! -d "${myDir}"/05-markduplicates ]
        then mkdir ${myDir}/05-markduplicates
fi

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

        scriptName=${myDir}/temp/${scriptBase}${index}.sh
        rm -rf ${scriptName} || true
        touch ${scriptName}

        echo "#!/bin/bash" >> ${scriptName}
        echo "#SBATCH --partition=${queue}" >> ${scriptName}
        echo "#SBATCH --mem=${mem}" >> ${scriptName}
        echo "#SBATCH --nodes=1" >> ${scriptName}
        echo "#SBATCH --tasks-per-node=${cpu}" >> ${scriptName}
        echo "#SBATCH -t ${runTime}" >> ${scriptName}
        echo "#SBATCH --array=1-${nSamps}" >> ${scriptName}
        echo "#SBATCH --output ${myDir}/OUT/${scriptBase}${jobName}.%A-%a" >> ${scriptName}
        echo "#SBATCH --error ${myDir}/ERR/${scriptBase}${jobName}.%A-%a" >> ${scriptName}
        echo "#SBATCH --account=${userProject}" >> ${scriptName}

        echo "module load ${PICARD}" >> ${scriptName}
        echo "module load ${JAVA}" >> ${scriptName}
        echo "module load ${SAMTOOLS}" >> ${scriptName}

        ## specify the file to work on
        echo "sampleIter=\$(( \${SLURM_ARRAY_TASK_ID} + ${loopSetup[$index]} ))
        sampleID=\$(cat ${myDir}/01-download/SampleFileNames.txt | sed -n \${sampleIter}p)
        echo \${sampleID}" >> ${scriptName}

        ## run the markduplicate and samtools sort commands

	echo "java -jar /trinity/shared/apps/site-local/${PICARD}/picard.jar MarkDuplicates \
	I=${myDir}/04-star/\${sampleID}onemap_Aligned.sortedByCoord.out.bam \
	O=${myDir}/05-markduplicates/\${sampleID}_markdup.bam \
	M=${myDir}/05-markduplicates/\${sampleID}_metrics_markdup.txt \
	REMOVE_DUPLICATES=false \
	VALIDATION_STRINGENCY=SILENT" >> ${scriptName}

	echo "java -jar /trinity/shared/apps/site-local/${PICARD}/picard.jar MarkDuplicates \
	I=${myDir}/04-star/\${sampleID}onemap_Aligned.sortedByCoord.out.bam \
	O=${myDir}/05-markduplicates/\${sampleID}_rmdup.bam \
	M=${myDir}/05-markduplicates/\${sampleID}_metrics_rmdup.txt \
	REMOVE_DUPLICATES=true \
	VALIDATION_STRINGENCY=SILENT" >> ${scriptName}

	echo "samtools index ${myDir}/05-markduplicates/\${sampleID}_markdup.bam" >> ${scriptName}
	echo "samtools index ${myDir}/05-markduplicates/\${sampleID}_rmdup.bam" >> ${scriptName}

        chmod u+x ${scriptName}

        sbatch ${scriptName}

	index=$(( $index + 1))
done

exit 0
