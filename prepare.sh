#!/bin/bash

# split by chromosome
tail -n +2 NAM_Map_20090730.txt |egrep -o '^[[:digit:]]+'|uniq | xargs -I{} bash -c "egrep '^{}[[:space:]]' NAM_Map_20090730.txt | cut -f1,3- | sed -re 's/m(0)*//' > map.{}.txt"
#tail -n +2 NAM_Map_20090730.txt | cut -f1,3-|sed -re 's/m(0)*//' > allmap.txt

ls fastphase_chr*.txt | xargs -I{} bash -c "tail -n +2 {} | cut -f3,4,12- > {}.filtered"

# format: MO###; the rest is the same as the GWAS markers
ls imputedIBMMarkers*.txt | xargs -I{} bash -c "tail -n +2 {} | sed -re 's/^MO(0)*/17\t/' > {}.filtered"

ls imputedMarkersGWAS.chr*.082809.txt | xargs -I{} bash -c "tail -n +2 {}|sed -re 's/Z([[:digit:]]+)E([[:digit:]]+)/\1\t\2/' \
|sed -re 's/([[:space:]]|^)(0)*([1-9.]+)/\1\3/g' > {}.filtered"

# merge GWAS and IBM markers
ls imputed*chr*.txt.filtered |egrep -o 'chr[[:digit:]]+'|sort | uniq | xargs -I{} bash -c "cat imputed*{}.*.filtered > imputedMarkers.{}.merged"

# filter residuals
tail -n +2 residuals_distancebased_for_asi.txt|sed -re 's/Z([[:digit:]]+)E([[:digit:]]+)/\1\t\2/' | sed -re 's/([[:space:]]|^)(0)*([1-9.]+)/\1\3/g'| sed -re 's/^MO(0)*/17\t/' > residuals_distancebased_for_asi.txt.filtered
