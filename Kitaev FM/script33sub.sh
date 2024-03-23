#!/bin/bash
#
#$ -l hostname=micro4.local 
#$ -N sub33_S1
#$ -l h_rt=180:00:00                   #estimate max run time
#$ -q all.q
#$ -m ea
#$ -M lucasrdf@fisica.ufrn.br
#$ -cwd
#$ -o /home/rodrigues/ED-QM-Impurity/out
#$ -e /home/rodrigues/ED-QM-Impurity/err
#
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/sge/lib/lx-amd64
export CLASSPATH=/opt/sge/lib/drmaa.jar:/opt/sge/lib/juti.jar:/opt/sge/lib/jgdi.jar
#
math -script substitution33.m
