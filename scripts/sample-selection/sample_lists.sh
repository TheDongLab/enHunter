
# listing out the minus strand big wigs
gsutil -u amp-pd-dawg-analyses ls gs://fc-secure-4300ea7f-8e77-4f20-a20b-38417286eaaf/submissions/ae05f833-cac2-42e4-9cd9-24fb0cd797ed/bam2bigwig_workflow/'*'/call-bam2bigwig/attempt-2/'*'.minus.normalized.bw >> fixed-samplesN200.bigwig.minus.list.txt 
gsutil -u amp-pd-dawg-analyses ls gs://fc-secure-4300ea7f-8e77-4f20-a20b-38417286eaaf/submissions/ae05f833-cac2-42e4-9cd9-24fb0cd797ed/bam2bigwig_workflow/'*'/call-bam2bigwig/'*'.minus.normalized.bw >> fixed-samplesN200.bigwig.minus.list.txt

# listing out the plus strand big wigs 
gsutil -u amp-pd-dawg-analyses ls gs://fc-secure-4300ea7f-8e77-4f20-a20b-38417286eaaf/submissions/ae05f833-cac2-42e4-9cd9-24fb0cd797ed/bam2bigwig_workflow/'*'/call-bam2bigwig/attempt-2/'*'.plus.normalized.bw >> fixed-samplesN200.bigwig.plus.list.txt
gsutil -u amp-pd-dawg-analyses ls gs://fc-secure-4300ea7f-8e77-4f20-a20b-38417286eaaf/submissions/ae05f833-cac2-42e4-9cd9-24fb0cd797ed/bam2bigwig_workflow/'*'/call-bam2bigwig/'*'.plus.normalized.bw >> fixed-samplesN200.bigwig.plus.list.txt

# creating table of 8356 - 200 
grep -v -f n200-new-sample-names.tsv new-sample.tsv > n8356-200_new-sample.tsv 

lines=200
# selecting random n200 samples from CONTROL 

# selecting random n200 samples from PD 

sort -R all_PD_cases.txt | head -n $lines
# https://stackoverflow.com/questions/9245638/select-random-lines-from-a-file 


#### creating all bigwig minus list #####
gsutil ls gs://fc-secure-4300ea7f-8e77-4f20-a20b-38417286eaaf/submissions/72e10f97-bfdd-4d9e-afc2-9ad9a1b19325/bam2bigwig_workflow/'*'/call-bam2bigwig/attempt-2/'*'.minus.normalized.bw > all_minus_bigwigs.txt
gsutil ls gs://fc-secure-4300ea7f-8e77-4f20-a20b-38417286eaaf/submissions/72e10f97-bfdd-4d9e-afc2-9ad9a1b19325/bam2bigwig_workflow/'*'/call-bam2bigwig/'*'.minus.normalized.bw >> all_minus_bigwigs.txt

# for some reason there are two extra lines ... 



test=$(head -1 subset.wide_hu.xls| tr '\t' '\n' | nl -b a | grep -F -f baseline.tsv|cut -f1|tr -d "[:blank:]" |awk 1 ORS="," |sed 's/,$//')
cut -f $test,1 subset.wide_hu.xls > testing.xls
# TODO find a way to select the exonic RNA

#cut -f1 temp.txt| tr -d "[:blank:]" |awk 1 ORS="," |sed 's/,$//'| head




awk -F'/' '{split($NF, a, "."); print a[1]}' all
#sed 's//' head -1 subset.wide_hu.xls| tr '\t' '\n' | nl | //g' PPMI_baseline_samples.tsv | awk 1 ORS="," 

sample=$(head PPMI_baseline_samples.tsv | tr -d '^M' | awk 1 ORS="" | sed 's/,$//')
awk -v s=$sample 'NR==1{split(s,a,","); for(i=1; i<=NF; i++) if ($i in a){values[i] = i}; next} END{for (key in values) {print values[key]} }' subset.wide_hu.xlsm