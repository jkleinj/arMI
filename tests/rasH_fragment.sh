#! /bin/sh
#===============================================================================
# arMI computation on human Ras fragment 
#===============================================================================
../src/armi --mali rasH_fragment.fasta --nsubset 100 --prefix rasH_fragment_ || exit 1
