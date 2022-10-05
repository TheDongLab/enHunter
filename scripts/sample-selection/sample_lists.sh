
# listing out the minus strand big wigs
gsutil -u amp-pd-dawg-analyses ls gs://fc-secure-4300ea7f-8e77-4f20-a20b-38417286eaaf/submissions/ae05f833-cac2-42e4-9cd9-24fb0cd797ed/bam2bigwig_workflow/'*'/call-bam2bigwig/attempt-2/'*'.minus.normalized.bw >> fixed-samplesN200.bigwig.minus.list.txt 
gsutil -u amp-pd-dawg-analyses ls gs://fc-secure-4300ea7f-8e77-4f20-a20b-38417286eaaf/submissions/ae05f833-cac2-42e4-9cd9-24fb0cd797ed/bam2bigwig_workflow/*/call-bam2bigwig/*.minus.normalized.bw >> fixed-samplesN200.bigwig.minus.list.txt

# listing out the plus strand big wigs 
gsutil ls gs://fc-secure-4300ea7f-8e77-4f20-a20b-38417286eaaf/submissions/ae05f833-cac2-42e4-9cd9-24fb0cd797ed/bam2bigwig_workflow/*/call-bam2bigwig/attempt-2/*.plus.normalized.bw >> fixed-samplesN200.bigwig.plus.list.txt
gsutil ls gs://fc-secure-4300ea7f-8e77-4f20-a20b-38417286eaaf/submissions/ae05f833-cac2-42e4-9cd9-24fb0cd797ed/bam2bigwig_workflow/*/call-bam2bigwig/*.plus.normalized.bw >> fixed-samplesN200.bigwig.plus.list.txt

lines=200
# selecting random n200 samples from CONTROL 

# selecting random n200 samples from PD 

sort -R all_PD_cases.txt | head -n $lines
# https://stackoverflow.com/questions/9245638/select-random-lines-from-a-file 