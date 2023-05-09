#!/bin/bash

#This script will quality check the raw and trimmed data using fastQC. It relies on correctly setting the variables indicated below.
#Written by Sarah Christofides, 2022, based on a script from Rob Andrews. Released under Creative Commons BY-SA.

###VARIABLES TO BE SET###
#Set the path to your directory on scratch
myDir=
#Set your username
userProject=
#Indicate if the data is single-end (SE) or paired-end (PE)
ends=SE
#Set the slurm queue to use: defq for gomphus, epyc for iago, htc for hawk
queue=epyc
######

mem="1G"
cpu="1"
runTime="00:05:00"
scriptBase="03QC"
slurmids=""

#Set the correct version of FastQC and multiqc
FASTQC=$(module avail -L fastqc/ | tail -n 1)
multiQC=$(module avail -L multiqc/ | tail -n 1)
##Append this information to the report
echo -e "\nQuality assessed with ${FASTQC} and ${multiQC}" >> ${myDir}/AnalysisReport.txt

#Make directories for the output
DIRLIST=("03-fastqc" "03-fastqc/raw" "03-fastqc/trimmed" "03-fastqc/multiqc")
for DIRECTORY in "${DIRLIST[@]}"
do
        if [ ! -d "${myDir}/$DIRECTORY" ]; then
          mkdir ${myDir}/${DIRECTORY}
        fi
done

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

        ## make an empty script for writing
        touch ${scriptName}

        ## write the SLURM parameters to the top of the script
        echo "#!/bin/bash" > ${scriptName}
        echo "#SBATCH --partition=${queue}" >> ${scriptName}
        echo "#SBATCH --mem-per-cpu=${mem}" >> ${scriptName}
        echo "#SBATCH --nodes=1" >> ${scriptName}
        echo "#SBATCH --ntasks=1" >> ${scriptName}
        echo "#SBATCH --cpus-per-task=${cpu}" >> ${scriptName}
        echo "#SBATCH --time=${runTime}" >> ${scriptName}
        echo "#SBATCH --array=1-${nSamps}" >> ${scriptName}
        echo "#SBATCH --output ${myDir}/OUT/${scriptBase}${jobName}.%A-%a" >> ${scriptName}
        echo "#SBATCH --error ${myDir}/ERR/${scriptBase}${jobName}.%A-%a" >> ${scriptName}

	echo "module load ${FASTQC}" >> ${scriptName}	

	## Specify the file to work on

        echo "sampleIter=\$(( \${SLURM_ARRAY_TASK_ID} + ${loopSetup[$index]} ))
        sampleID=\$(cat ${myDir}/01-download/SampleFileNames.txt | sed -n \${sampleIter}p)
        echo \${sampleID}

        ## find the raw datafile
        FILEPATH=\$(find ${myDir}/01-download/ -name \${sampleID}R1* | sed -E 's/(.+)\/.+/\1/')" >> ${scriptName}


	## run fastqc on the raw fastq

	echo -n "fastqc -o ${myDir}/03-fastqc/raw/ \${FILEPATH}/\${sampleID}R1*" >> ${scriptName}
        if [ "${ends}" = PE ]; then echo " \${FILEPATH}/\${sampleID}R2*" >> ${scriptName}; fi
	echo -e "\n" >> ${scriptName}

	## run fastqc on the trimmed fastq

	echo -n "fastqc -o ${myDir}/03-fastqc/trimmed/ ${myDir}/02-trim/trim_\${sampleID}_R1*" >> ${scriptName}
        if [ "${ends}" = PE ]; then echo " ${myDir}/02-trim/trim_\${sampleID}_R2*" >> ${scriptName}; fi

	echo -e "\nexit 0" >> ${scriptName}

	chmod u+x ${scriptName}

        slurmids="${slurmids}:$(sbatch --parsable ${scriptName})"

        index=$(( $index + 1))

done


#Create and run a slurm script that will do multiQC on both sets of FastQCs
scriptName=${myDir}/temp/${scriptBase}.sh
touch ${scriptName}

echo "#!/bin/bash" > ${scriptName}
echo "#SBATCH --partition=defq" >> ${scriptName}
echo "#SBATCH --mem-per-cpu=${mem}" >> ${scriptName}
echo "#SBATCH --nodes=1" >> ${scriptName}
echo "#SBATCH --tasks-per-node=1" >> ${scriptName}
echo "#SBATCH --output ${myDir}/OUT/${scriptBase}${jobName}.%J" >> ${scriptName}
echo "#SBATCH --error ${myDir}/ERR/${scriptBase}${jobName}.%J" >> ${scriptName}

echo "module load ${multiQC}" >> ${scriptName}

#Run multiQC on the raw fastQC
echo "multiqc ${myDir}/03-fastqc/raw/ -o ${myDir}/03-fastqc/multiqc/ -i Raw" >> ${scriptName}

#Run multiQC on the trimmed fastQC
echo "multiqc ${myDir}/03-fastqc/trimmed/ -o ${myDir}/03-fastqc/multiqc -i Trimmed" >> ${scriptName}

chmod u+x ${scriptName}

sbatch -d afterok${slurmids} ${scriptName}

exit 0
