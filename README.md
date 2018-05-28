# INTERSECTparser
calculate percentage intersected using BEDTools intersect reports

## obtain gene cds bed (target)<br>
cat UCSC1.final.noSEQ.gff | awk -v OFS='\t' -v FS='\t' '$2=="maker"&&$3=="CDS" {print $1, $4-1, $5, $1 ":" $4 "-" $5 ";" $9, "1", $7}' >maker.bed
<br>

## obtain feature bed (queries)<br>
for feature in \`echo  repeatmasker protein2genome genemark ori_snap ori_augustus\`;
do 
cat UCSC1.final.noSEQ.gff | awk -v OFS='\t' -v FS='\t' '$2~"'$feature'"&&$3=="match_part" {print $1, $4-1, $5, $1 ":" $4 "-" $5 ";" $9, "1", $7}' >$feature.bed;
done
<br>

## bedtool intersect target v.s each queries<br>
* repeatmasker (both strand)
for feature in \`echo repeatmasker\`; do ./bedtools intersect -wo -a maker.bed -b $feature.bed >maker.$feature.intersect; done<br>
* other features (strickly stranded)
for feature in \`echo protein2genome genemark ori_snap ori_augustus\`; do ./bedtools intersect -s -wo -a maker.bed -b $feature.bed >maker.$feature.intersect; done<br>

## put bed files to INTERSECTparser/bed/<br>

## put intersect files to INTERSECTparser/intersect/<br>

## run main python script to get percentage of coverage of CDS region for each gene<br>
cd INTERSECTparser/<br>
python3 compare.py<br>
