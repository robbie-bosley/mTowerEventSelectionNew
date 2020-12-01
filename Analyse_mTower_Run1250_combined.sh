#!/bin/bash

#run with run number as argument, e.g. ./Analyse_mTower.sh 1250

#to submit this script in the quark queue
#        qsub -cwd -V ./Analyse_mTower.sh 1250
# -cwd uses the current working directory
# -V keeps the environment variables

echo "====================================="
echo "starting job Analyse_mTower($1)"
echo "====================================="

root -l -b <<SHELL
.!pwd
.L /afs/cern.ch/user/r/rbosley/public/EventSelectionFinal/classes/mTowerHit.cxx+
.L /afs/cern.ch/user/r/rbosley/public/EventSelectionFinal/classes/mTowerEvent.cxx+
.L /afs/cern.ch/user/r/rbosley/public/EventSelectionFinal/classes/mTowerCluster.cxx+
.L /afs/cern.ch/user/r/rbosley/public/EventSelectionFinal/classes/mTowerChip.cxx+
.L /afs/cern.ch/user/r/rbosley/public/EventSelectionFinal/Analyse_mTower_combine.cxx+
Analyse_mTower(1250)
.q
SHELL
