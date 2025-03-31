#!/usr/bin/bash

export PATH=/apps/statistics2/R-3.6.0/bin/:$PATH
MYINPUT=$@
Rscript ${MYINPUT[*]}
