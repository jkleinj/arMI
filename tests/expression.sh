#! /bin/sh
#===============================================================================
# arMI computation on gene expression samples 
#===============================================================================
valgrind ../src/armiexp --expr expr_tf.dat --row rownames.dat --col colnames.dat --prefix expr_ || exit 1
