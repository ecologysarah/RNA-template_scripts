#!/bin/bash

#This script will quality trim and remove adapters from sequencing data using fastp (Chen et al 2018). It relies on correctly setting the variables indicated below.
#Written by Sarah Christofides, 2022. Released under Creative Commons BY-SA.

###VARIABLES TO BE SET###
#Set the path to your directory on scratch
myDir=
#Set your username
userProject=
#Indicate if the data is single-end (SE) or paired-end (PE)
ends=SE
#Indicate if you want to specify an adapter for SE reads (set as "" to omit)
adapt=""
#Set the slurm queue to use: defq for gomphus, epyc for iago, htc for hawk
queue=epyc
######

sampleIDs=$(cat ${myDir}/01-download/SampleFileNames.txt)

mem="2G"
cpu="3"
runTime="00:05:00"
scriptBase="02trim"

#Set the correct version of fastp
FASTP=$(module avail -L fastp/ | tail -n 1)
##Append this information to the report
echo -e "\nSequences quality trimmed with ${FASTP}" >> ${myDir}/AnalysisReport.txt

#Make the output directory
if [ ! -d "${myDir}"/02-trim ]
        then mkdir ${myDir}/02-trim
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

        ## load the module
        echo "module load $FASTP" >> ${scriptName}

        ## Specify the file to work on

        echo "sampleIter=\$(( \${SLURM_ARRAY_TASK_ID} + ${loopSetup[$index]} ))
        sampleID=\$(cat ${myDir}/01-download/SampleFileNames.txt | sed -n \${sampleIter}p)
        echo \${sampleID}

        ## find the file
        FILEPATH=\$(find ${myDir}/01-download/ -name \${sampleID}R1* | sed -E 's/(.+)\/.+/\1/')

	## ascertain the format of the end of the file name
	FILEND=\$(ls \${FILEPATH}/\${sampleID}R1* | sed -E 's/.+_R1(.+)/\1/')" >> ${scriptName}

        ## run the command
        echo "fastp -h ${myDir}/02-trim/trim_\${sampleID}.html \\
	-j ${myDir}/02-trim/trim_\${sampleID}.json \\
        -i \${FILEPATH}/\${sampleID}R1\${FILEND} \\" >> ${scriptName}
        if [ ! "${adapt}" = "" ]; then echo "-a ${adapt} \\" >> ${scriptName}; fi
        if [ "${ends}" = PE ]; then echo "-I \${FILEPATH}/\${sampleID}R2\${FILEND} \\" >> ${scriptName}; fi
        echo "-o ${myDir}/02-trim/trim_\${sampleID}_R1\${FILEND} \\" >> ${scriptName}
        if [ "${ends}" = PE ]; then echo "-O ${myDir}/02-trim/trim_\${sampleID}_R2\${FILEND} \\" >> ${scriptName}; fi
	##Set the minimum phred score to 20, with no more than 10% bases allowed to drop below that
	echo "-q 20 \\
	-u 10 \\" >> ${scriptName}
	## Enable sliding window trim
	echo "--cut_right" >> ${scriptName}

        ## make the script into an 'executable'
        chmod u+x ${scriptName}

        ## submit the script to the compute queue
        sbatch ${scriptName}

        index=$(( $index + 1))

done

exit 0

