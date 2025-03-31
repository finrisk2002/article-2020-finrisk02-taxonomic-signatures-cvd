#!/bin/bash

for f in ~/data/tmp/BAM/*.bam
do
  echo $f;
  samtools depth -a $f > $f.coverage
done	 

