
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
awk -F'/' '{split($NF, a, "."); print a[1]}' all