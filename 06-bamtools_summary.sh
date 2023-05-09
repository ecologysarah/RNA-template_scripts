#!/bin/bash

#This script will produced summary statistics for the mapping. It relies on correctly setting the variables indicated below.
#Written by Sarah Christofides, 2022, based on a script from Rob Andrews. Released under Creative Commons BY-SA.

###VARIABLES TO BE SET###
#Set the path to your directory on scratch
myDir=
#Set your username
userProject=
######

mem="1G"
cpu="1"
runTime="00:05:00"
scriptBase="06bamstats"
slurmids=""

#Set the correct version of BAMtools
BAMTOOLS=$(module avail -L bamtools/ | tail -n 1)
#Set the correct version of multiqc
multiQC=$(module avail -L multiqc/ | tail -n 1)
##Append this information to the report
echo -e "\nData summaries produced with ${BAMTOOLS}" >> ${myDir}/AnalysisReport.txt

#Make the output directory
if [ ! -d "${myDir}"/06-bamtools ]
        then mkdir ${myDir}/06-bamtools
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
        touch ${scriptName}

        echo "#!/bin/bash" > ${scriptName}
        echo "#SBATCH --partition=epyc" >> ${scriptName}
        echo "#SBATCH --mem=${mem}" >> ${scriptName}
        echo "#SBATCH --nodes=1" >> ${scriptName}
        echo "#SBATCH --tasks-per-node=${cpu}" >> ${scriptName}
        echo "#SBATCH -t ${runTime}" >> ${scriptName}
        echo "#SBATCH --array=1-${nSamps}" >> ${scriptName}
        echo "#SBATCH --output ${myDir}/OUT/${scriptBase}${jobName}.%A-%a" >> ${scriptName}
        echo "#SBATCH --error ${myDir}/ERR/${scriptBase}${jobName}.%A-%a" >> ${scriptName}
        echo "#SBATCH --account=${userProject}" >> ${scriptName}

        echo "module load ${BAMTOOLS}" >> ${scriptName}

        ## specify the file to work on

	echo "sampleIter=\$(( \${SLURM_ARRAY_TASK_ID} + ${loopSetup[$index]} ))
        sampleID=\$(cat ${myDir}/01-download/SampleFileNames.txt | sed -n \${sampleIter}p)
        echo \${sampleID}" >> ${scriptName}

        ## run bamtools 

        echo "bamtools stats -in ${myDir}/05-markduplicates/\${sampleID}_markdup.bam > ${myDir}/06-bamtools/\${sampleID}_markdup_dupstats.txt" >> ${scriptName}
        echo "bamtools stats -in ${myDir}/05-markduplicates/\${sampleID}_rmdup.bam > ${myDir}/06-bamtools/\${sampleID}_rmdup_dupstats.txt" >> ${scriptName}

        chmod u+x ${scriptName}

        slurmids="${slurmids}:$(sbatch --parsable ${scriptName})"

	index=$(( $index + 1))
done

#Create and run a slurm script that will do multiQC on the output
scriptName=${myDir}/temp/${jobName}.sh
touch ${scriptName}

echo "#!/bin/bash" > ${scriptName}
echo "#SBATCH --partition=epyc" >> ${scriptName}
echo "#SBATCH --mem-per-cpu=2" >> ${scriptName}
echo "#SBATCH --nodes=1" >> ${scriptName}
echo "#SBATCH --tasks-per-node=${cpu}" >> ${scriptName}
echo "#SBATCH --output ${myDir}/OUT/${scriptBase}${jobName}.%J" >> ${scriptName}
echo "#SBATCH --error ${myDir}/ERR/${scriptBase}${jobName}.%J" >> ${scriptName}

echo "module load ${multiQC}" >> ${scriptName}

#Run multiQC
echo "multiqc ${myDir}/06-bamtools -o ${myDir}/03-fastqc/multiqc/ -i MarkDuplicates" >> ${scriptName}

chmod u+x ${scriptName}

sbatch -d afterok${slurmids} ${scriptName}

exit 0
