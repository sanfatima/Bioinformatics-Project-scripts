
#!/bin/sh
#$ -V
#$ -cwd
#$ -l h_rt=40:00:0
#$ -pe smp 2
#$ -t 2-35
#$ -l h_vmem=10G
#$ -m bea
#$ -M bt20810@qmul.ac.uk



source /data/home/hfx365/.bashrc
source activate lohhla_env
# Define the folder for bam files and all other data
tableFile=/data/BCI-EvoCa2/syeda/HLA/samples.txt
bamFolder=/data/BCI-EvoCa2/syeda/HLA/BAM
baseFolder=/data/BCI-EvoCa2/syeda/Sequenza
# Read in file infos from tab-separated table
namePrefix=$(awk -v var="$SGE_TASK_ID" 'NR ==var { OFS="\t";print $1}' $tableFile)
# Create output folder to store results in
outputFolder=$baseFolder/LOHHLAfinalresults/${namePrefix}
echo '----Creating output folder: '$outputFolder'------';
mkdir -p $outputFolder
# Run LOHHLA
Rscript /data/home/hfx365/Software/lohhla/lohhla/LOHHLAscript_verbose.R --patientId ${namePrefix} --outputDir ${outputFolder} --normalBAMfile ${bamFolder}/${namePrefix}/${namePrefix}_normal.mkdup.bam --BAMDir ${bamFolder}/${namePrefix} --CopyNumLoc ${baseFolder}/tumour_purity.txt --hlaPath /data/BCI-EvoCa2/syeda/HLA/HLA_types/${namePrefix}/winners.hla.txt --HLAfastaLoc /data/home/hfx365/Software/hla-polysolver/data/abc_complete.fasta --mappingStep TRUE --minCoverageFilter 10 --fishingStep TRUE --cleanUp FALSE --gatkDir /data/home/hfx365/anaconda2/envs/polysolver_env/jar --novoDir /data/home/hfx365/anaconda2/envs/lohhla_env/bin/ --HLAexonLoc /data/home/hfx365/Software/lohhla/lohhla/data/hla.dat --genomeBuild hg19 --BAMpattern *.mkdup.bam


