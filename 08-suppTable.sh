#!/bin/bash

#This script will tar compress and archive your data. The OUT, ERR and temp directories are not included. It relies on correctly setting the variables indicated below.
#Written by Sarah Christofides, 2022. Released under Creative Commons BY-SA.

###VARIABLES TO BE SET###
#Set the path to your directory on scratch - do not include a trailing /
myDir=
#Set your username
userProject=
#Were duplicates removed (rm) or not (mark)?
DUP=mark
#Set the name and path for the tar file (do not include the file extension)
TAR=
#Set the path to the archive location. N.B. If this is a remote server, it will need a ssh key pair in place.
ARCHIVE=
#Set the slurm queue to use: defq for gomphus, epyc for iago, htc for hawk
queue=mammoth
######

mem="40G"
cpu="2"
runTime="05:00:00"
scriptBase="08-summaryArchive"

##Create the slurm script
scriptName=${myDir}/temp/${scriptBase}.sh
rm -rf ${scriptName} || true
touch ${scriptName}

echo "#!/bin/bash" >> ${scriptName}
echo "#SBATCH --partition=${queue}" >> ${scriptName}
echo "#SBATCH --mem-per-cpu=${mem}" >> ${scriptName}
echo "#SBATCH --nodes=1" >> ${scriptName}
echo "#SBATCH --cpus-per-task=${cpu}" >> ${scriptName}
echo "#SBATCH --time=${runTime}" >> ${scriptName}
echo "#SBATCH --output ${myDir}/OUT/${scriptBase}.%J" >> ${scriptName}
echo "#SBATCH --error ${myDir}/ERR/${scriptBase}.%J" >> ${scriptName}

##Create the summary table
echo "sampleIDs=\$(cat ${myDir}/01-download/SampleFileNames.txt)" >> ${scriptName}

echo -e "echo -e SampleID\\tRawReads\\tReadsPassingQC\\tNumberUniquelyMapped\\tPercentUniquelyMapped\\tNumberMultimapped\\tPercentMultimapped\\tTotalMapped\\tPercentMapped\\tDeduplicatedAndMapped\\tPercentDupliction\\tEstLibSize\\tAssignedToGene\\tPercentAssignedToGene > ${myDir}/SuppTable.txt" >> ${scriptName}

echo -e "
for sampleID in \$sampleIDs
do
nraw=\$(cat ${myDir}/03-fastqc/multiqc/Raw_multiqc_report_data/multiqc_general_stats.txt | grep \${sampleID}R1 | cut -f 6)
nstart=\$(cat ${myDir}/04-star/\${sampleID}onemap_Log.final.out | grep 'Number of input reads' | cut -f 2)
nuniq=\$(cat ${myDir}/04-star/\${sampleID}onemap_Log.final.out | grep 'Uniquely mapped reads number' |  cut -f 2)
puniq=\$(cat ${myDir}/04-star/\${sampleID}onemap_Log.final.out | grep 'Uniquely mapped reads %' | cut -f 2)
nmulti=\$(cat ${myDir}/04-star/\${sampleID}onemap_Log.final.out | grep 'Number of reads mapped to multiple loci' | cut -f 2)
pmulti=\$(cat ${myDir}/04-star/\${sampleID}onemap_Log.final.out | grep '% of reads mapped to multiple loci' | cut -f 2)
ntotal=\$((\$nuniq+\$nmulti))
ptotal=\$((\${ntotal}*100))
ptotal=\$(echo \$((\${ptotal} / \${nstart}))%)
ndedup=\$(cat ${myDir}/06-bamtools/\${sampleID}_rmdup_dupstats.txt | grep 'Total reads:' | awk '{print \$3}')
pdup=\$(cat ${myDir}/06-bamtools/\${sampleID}_markdup_dupstats.txt  | grep 'Duplicates:' | awk '{print \$3}' | sed -E 's/[\(\)]//g')
estlib=\$(cat ${myDir}/05-markduplicates/\${sampleID}_metrics_markdup.txt | grep -A 1 'ESTIMATED_LIBRARY_SIZE' |  awk '{print \$10}' | tail -n 1)
nassigned=\$(cat ${myDir}/07-featurecounts/\${sampleID}_${DUP}dup.featurecount.summary | grep 'Assigned' | awk '{print \$2}')
assigTot=\$(cat ${myDir}/07-featurecounts/\${sampleID}_${DUP}dup.featurecount.summary | awk '{sum+=\$2} END {print sum}')
passigned=\$(echo -e \$((\$nassigned*100/\$assigTot))%)

echo -e \$sampleID\\t\$nraw\\t\$nstart\\t\$nuniq\\t\$puniq\\t\$nmulti\\t\$pmulti\\t\$ntotal\\t\$ptotal\\t\$ndedup\\t\$pdup\\t\$estlib\\t\$nassigned\\t\$passigned >> ${myDir}/SuppTable.txt
done
" >> ${scriptName}

##Create the tar archive
echo "tar -czvf ${myDir}/${TAR}.tar.gz --exclude={"${myDir}/ERR","${myDir}/OUT","${myDir}/temp"} ${myDir}" >> ${scriptName}

##Transfer it to the backup location
echo "rsync ${myDir}/${TAR}.tar.gz ${ARCHIVE}" >> ${scriptName}

echo "exit 0" >> ${scriptName}

##Make the script executable
chmod u+x ${scriptName}

##Submit the script to the compute queue
sbatch ${scriptName}

exit 0

