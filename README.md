# INTERSECTparser
calculate percentage intersected using BEDTools intersect reports

## U.list
UCSC1<br>
UMSG1<br>
UMSG2<br>
UMSG3<br>

## obtain gene cds bed (target; A)<br>
for i in \`cat U.list\`; do cat ${i}_5th.gff | awk -v OFS='\t' -v FS='\t' '$2=="maker"&&$3=="CDS" {print $1, $4-1, $5, $1 ":" $4 "-" $5 ";" $9, "1", $7}' >$i.A.maker.bed; done<br>

## obtain gene cds bed (control; B)<br>
for i in \`cat U.list\`; do cat ${i}_6th.gff | awk -v OFS='\t' -v FS='\t' '$2=="maker"&&$3=="CDS" {print $1, $4-1, $5, $1 ":" $4 "-" $5 ";" $9, "1", $7}' >$i.B.maker.bed; done<br>

## obtain gene cds bed (control; C)<br>
for i in \`cat U.list\`; do cat ${i}_maker_final_rename_only_rmhost_addorf_secretome_mod.gff | awk -v OFS='\t' -v FS='\t' '$2=="maker"&&$3=="CDS" {print $1, $4-1, $5, $1 ":" $4 "-" $5 ";" $9, "1", $7}' >$i.C.maker.bed; done<br>

## obtain repeatmasker bed (query)<br>
for i in \`cat U.list\`; do for feature in \`echo repeatmasker\`; do cat ${i}_5th.gff | grep -v Low_complexity | grep -v Simple_repeat | awk -v OFS='\t' -v FS='\t' '$2~"'$feature'"&&$3=="match_part" {print $1, $4-1, $5, $1 ":" $4 "-" $5 ";" $9, "1", $7}' >$i.A.$feature.bed; done; done
<br>


## obtain feature bed (queries)<br>
for i in \`cat U.list\`; do
for feature in \`echo  repeatmasker protein2genome genemark ori_snap ori_augustus est2genome\`;
do 
cat ${i}_5th.gff | awk -v OFS='\t' -v FS='\t' '$2~"'$feature'"&&$3=="match_part" {print $1, $4-1, $5, $1 ":" $4 "-" $5 ";" $9, "1", $7}' >$i.A.$feature.bed;
done;
done
<br>

## bedtool intersect target v.s each queries<br>
* repeatmasker (both strand)<br>
for i in \`cat U.list\`; do
for feature in \`echo repeatmasker\`; do ./bedtools intersect -wo -a $i.A.maker.bed -b $i.A.$feature.bed >$i.A.maker.$feature.intersect; done;
done<br>
* other features (strickly stranded)<br>
for i in \`cat U.list\`; do
for feature in \`echo protein2genome genemark ori_snap ori_augustus est2genome\`; do ./bedtools intersect -s -wo -a $i.A.maker.bed -b $i.A.$feature.bed >$i.A.maker.$feature.intersect; done;
done<br>

## put bed files to INTERSECTparser/bed/<br>

## put intersect files to INTERSECTparser/intersect/<br>

## run main python script to get percentage of coverage of CDS region for each feature<br>
cd INTERSECTparser/<br>
python3 compare.py<br>
