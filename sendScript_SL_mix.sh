#!/bin/bash

# Submission script for GridEngine (GE). Each job will 
# be executed via the jobScript.sh
# This jobScript supports up to 7 parameters. Edit 
# the user specific part of the script according to 
# your program.
#
# Input to the script is a filelist with 1 file per line.
# For each file a job is started. With the parameter 
# nFilesPerJob a comma separated filelist will be 
# generated and handed to the job script. This feature
# is usefull when running many small jobs. Each
# job has its own logfile. All needed directories for the
# logfiles will be created if non existing.
#
# IMPORTANT: the lustre/nyx/prometheus cluster jobs will only
# see the /lustre/nyx file system. All needed scripts, programs
# and parameters have to be located on /lustre/nyx or the job
# will crash. This script syncs your working dir to the submission
# dir on /lustre/nyx . Make sure your scripts use the submission dir!
# Software should be taken from /cvmfs/hades.gsi.de/install/
#                                                      
# job log files will be named like inputfiles. If nFilesPerJob > 1
# the log files will contain the partnumber.
#
######################################################################
#   CONFIGURATION

user=$(whoami)
currentDir=$(pwd)
#day=Sim

# 074 and 075 already submitted
day=${1}
mix=yes
submmissionbase=/lustre/nyx/hades/user/${user}/sub/GE/${day}_mix/
submissiondir=${submmissionbase}/gen2/
 nFilesPerJob=100 # for SIM or EXP extract
    jobscript=${submissiondir}/jobScript_SL.sh                        # exec script (full path, call without dot, set it executable!)
    outputdir=/lustre/nyx/hades/user/${user}/mar19/gee/day${day}_gen2_mix   # outputdir for files AND logFiles
#    outputdir=/lustre/nyx/hades/user/${user}/mar19/analysis_exp/sim_proposal_gen2_enCorr_MiX/   # outputdir for files AND logFiles
pathoutputlog=${outputdir}/out                                  # protocol from batch farm for each file
     filename=march19_${day}                             # filename of log file if nFilesPerJob > 1 (partnumber will be appended)
par1=/cvmfs/hades.gsi.de/install/5.34.34/hydra2-5.2/defall.sh     # optional par1 : environment script
#par1=/lustre/nyx/hades/user/rlalik/emc/profile.sh   # optional par1 : environment script
par2=${submissiondir}/emc_diphoton_ana_mix                               # optional par2 : executable
par3=""                                                         # optional par3 : input file list
par4=${outputdir}                                               # optional par4 : outputfile (part number will be appended (_num.root))
par5=1000000000                                                 # optional par5 : number of events
par6=${mix}                                                       # optional par6
par7=                                                       # optional par7
resources="--mem=6000 --time=0-6:00:00"                        # runtime < 10h, mem < 2GB
jobarrayFile="jobarray_${1}.dat"
#filelist=${currentDir}/day_${day}.list
#filelist=${currentDir}/day_no_enhancement_gcalor.list
#filelist=${currentDir}/ag165ag_gen2.list
filelist=${currentDir}/lists/day_${day}_ag158ag_3200A.list # file list in local dir! not in submissiondir!!!
#filelist=${currentDir}/ag165ag_gen2.list # file list in local dir! not in submissiondir!!!
######################################################################




createJobarray=yes # create jobarry if yes
#---------------------------------------------------------------------
# create a file list for submission (simulation, testing etc.)
# for real data you will have a filelist with real filenames
#---------------------------------------------------------------------
nFiles=$( cat $filelist | wc -l)

#---------------------------------------------------------------------
# create needed dirs
if [ ! -d $submmissionbase ]
then
    echo "===> CREATE SUBMISSIONBASEDIR : $submmissionbase"
    mkdir -p $submmissionbase
else
    echo "===> USE SUBMISSIONBASEDIR : $submmissionbase"
fi

#---------------------------------------------------------------------
# output dirs

if [ ! -d $pathoutputlog ]
then
   echo "===> CREATE LOGDIR : $pathoutputlog"
   mkdir -p $pathoutputlog
else
   echo "===> USE LOGDIR : $pathoutputlog"
fi

if [ ! -d $outputdir ]
then
   echo "===> CREATE OUTPUTDIR : $outputdir"
   mkdir -p $outputdir
else
   echo "===> USE OUTPUTDIR : $outputdir"
fi

#if [ ! -d $outputdir/crash ]
#then
#   echo "===> CREATE CRASHDIR : $outputdir/crash"
#   mkdir -p $outputdir/crash
#fi
#---------------------------------------------------------------------


ctF=0          # counter for file number
ctJ=0          # counter for job number
partNumber=0   # counter for part number


if [ -f $jobarrayFile ]
then
  rm -f $jobarrayFile
fi

echo "===> CREATING JOB ARRAY FILE"
#---------------------------------------------------------------------
# read the files list into an job array
declare -a jobarray
ct1=0
for file in $(cat $filelist)
do
   jobarray[$ct1]=$file
   ((ct1+=1))
done
#---------------------------------------------------------------------


#---------------------------------------------------------------------
# loop over the job array and submit parts with
# nFilesPerJob to GE

while ((ctF<$nFiles))
do
     #---------------------------------------------------------------------
     # build comma separated file list
     # per job
     if [ $nFilesPerJob -gt 1 ]
     then
        infileList=${jobarray[${ctF}]}
        ((ctF+=1))
        for (( ctList=1;ctList<$nFilesPerJob; ctList++ ))
        do   	
            if [ $ctF -lt ${nFiles} ]
            then
               infileList="${infileList},${jobarray[${ctF}]}"
               ((ctF+=1))
            fi
        done
     else 
        infileList=${jobarray[${ctF}]}
        ((ctF+=1))
     fi
     #---------------------------------------------------------------------
     
     ((partNumber+=1))

     logfile="${pathoutputlog}/${filename}_${partNumber}.log"

     if [ $nFilesPerJob -eq 1 ]
     then
        file=$(basename ${infileList})
        logfile="${pathoutputlog}/${file}.log"
        par4=${outputdir}/${filename}_${file}.root
     else
        par4=${outputdir}/${filename}_${partNumber}.root
     fi
     
     if [ -f ${logfile} ]
     then
        rm -f ${logfile}
     fi

     ######################################################################
     #  SEND NEW JOB (USER SPECIFIC)
     
     par3=${infileList}

     #     defall.sh prog  filelist outfile  nev
     echo "${par1} ${par2} ${par3} ${par4} ${par5} ${par6} ${par7}" >>  $jobarrayFile
     

     ######################################################################
     
done
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# sync the local modified stuff 
# to the submission dir
echo "===> SYNC CURENTDIR TO SUBMISSIONDIR : rsync  -vHaz $currentDir ${submmissionbase}"
rsync  -vHaz $currentDir ${submmissionbase}/

syncStat=$?

if [ ! $syncStat -eq 0 ]
then
     echo "===> ERROR : SYNCHRONIZATION ENCOUNTERED PROBLEMS"
else

  echo "-------------------------------------------------"
  
  nFiles=$( cat $jobarrayFile | wc -l)
  ctsend=0
  block=1000
  while ((${ctsend} * ${block} < ${nFiles}))
  do
     ((start=${ctsend}*${block}+1))
     ((stop= ${start}+${block}-1))
     ((rest=${nFiles}-${start}))
     if [ $rest -le $block ]
     then
        ((stop=$start+$rest))
     fi

     command="--array=${start}-${stop} ${resources} -D ${submissiondir}  --output=${pathoutputlog}/slurm-%A_%a.out ${jobscript} ${submissiondir}/${jobarrayFile} ${pathoutputlog}"
     #echo $command
     sbatch $command

     ((ctsend+=1))
  done

  echo "${nFiles} jobs for day ${day} submitted"
fi

