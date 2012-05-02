#!/bin/bash

# split by chromosome
tail -n +2 NAM_Map_20090730.txt |egrep -o '^[[:digit:]]+'|uniq | xargs -I{} bash -c "egrep '^{}[[:space:]]' NAM_Map_20090730.txt | cut -f1,3- | sed -re 's/m(0)*//' > map.{}.txt"

cat map.* | cut -f2 | awk 'NR-1{if($1!=(_+1))print _}{_=$1}' |sort -n  > missing.txt

#tail -n +2 NAM_Map_20090730.txt | cut -f1,3-|sed -re 's/m(0)*//' > allmap.txt

ls fastphase_chr*.txt | xargs -I{} bash -c "tail -n +2 {} | cut -f3,4,12- > {}.filtered"

# format: MO###; the rest is the same as the GWAS markers
ls imputedIBMMarkers*.txt | xargs -I{} bash -c "tail -n +2 {} | sed -re 's/^MO/017\t/' > {}.filtered"

ls imputedMarkersGWAS.chr*.082809.txt | xargs -I{} bash -c "tail -n +2 {}|sed -re 's/Z([[:digit:]]+)E([[:digit:]]+)/\1\t\2/' > {}.filtered"

# merge GWAS and IBM markers
# the GWAS and IBM markers are independently sorted, and need to be merged.
# sort -m
ls imputed*chr*.txt.filtered |egrep -o 'chr[[:digit:]]+'|sort|uniq | xargs -I{} bash -c "sort -m imputed*{}.*.filtered > imputedMarkers.{}.merged"
#| sed -re 's/([[:space:]]|^)(0)*([1-9.]+)/\1\3/g'

# filter residuals
# the sorting could be faster; MO and Z sections are already sorted
tail -n +2 residuals_distancebased_for_asi.txt | grep -n Z018|cut -d':' -f1|head -1 > .first18
(head -n `cat .first18 ` residuals_distancebased_for_asi.txt; \
grep MO residuals_distancebased_for_asi.txt|sed -re 's/^MO(0)*/17\t/';  \
tail -n +2 residuals_distancebased_for_asi.txt | tail -n +`cat .first18` | grep -v MO ) | \
sed -re 's/Z([[:digit:]]+)E([[:digit:]]+)/\1\t\2/' | sed -re 's/([[:space:]]|^)(0)*([1-9.]+)/\1\3/g' > residuals_distancebased_for_asi.txt.filtered

#tail -n +2 residuals_distancebased_for_asi.txt|sed -re 's/Z([[:digit:]]+)E([[:digit:]]+)/\1\t\2/' | sed -re 's/^MO(0)*/17\t/' | sort -g +0 -2 | sed -re 's/([[:space:]]|^)(0)*([1-9.]+)/\1\3/g' > residuals_distancebased_for_asi.txt.filtered
