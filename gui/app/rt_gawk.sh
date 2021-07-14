#!/bin/bash

gawk -F\\\\t '!/^@/ { print }' $1/master.sam | gawk -F\\\\t '!and($2,256) && !and($2,2048) { print } END { print NR }' | gawk -F\\\\t '!and($2,4) { print }'| grep -v "|Hc" | awk 'NF>=2' > $1/master.gawk.sam
