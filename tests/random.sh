#! /bin/sh
#===============================================================================
# arMI computation on random sequences (created internally)
#===============================================================================
../src/armiali --random --nseq 100 --lseq 10 --prefix random_ || exit 1
