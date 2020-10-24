#!/bin/bash

prefix='VC101'

# input base directory
inpbasedir=`$pwd`

# input directory for this particular data
inpdir=$inpbasedir$prefix'/rawdata/'

# configuration file containing the execution parameters of HiCPro
configfile='configfile'

# directory which will contain the output results
outdir=$inpbasedir$prefix'/HiCPro'

HiC-Pro -i $inpdir -o $outdir -c $configfile


