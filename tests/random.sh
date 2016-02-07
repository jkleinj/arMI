#! /bin/sh
#===============================================================================
# arMI computation on random sequences (created internally)
#===============================================================================
../src/armi --random --nseq 100 --lseq 10 --prefix random_ || exit 1
