#!/usr/bin/env bash
study_id=$1
if [[ ! -f $2 ]]
then
   touch null.txt
   exclude="null.txt"
else
   exclude="$2"
fi 

rm -f samples.tsv

printf 'sample\tobject\tsample_supplier_name\tid_run\tis_paired_read\tstudy_id\tstudy\n' > samples.tsv

jq --arg study_id $study_id -n '{avus: [
       {attribute: "type", value: "cram", o: "="}, 
       {attribute: "manual_qc", value: "1", o: "="}, 
      {attribute: "target", value: "1", o: "="},
      {attribute: "study_id", value: $study_id, o: "="}]}' |\
/software/sciops/pkgg/baton/2.0.1+1da6bc5bd75b49a2f27d449afeb659cf6ec1b513/bin/baton-metaquery \
		--zone seq --obj --avu |\
jq '.[] as $a| 
"\($a.avus | .[] | select(.attribute == "sample") | .value)____\($a.collection)/\($a.data_object)____\($a.avus | .[] | select(.attribute == "sample_supplier_name") | .value)____\($a.avus | .[] | select(.attribute == "id_run") | .value)____\($a.avus | .[] | select(.attribute == "is_paired_read") | .value)____\($a.avus | .[] | select(.attribute == "study_id") | .value)____\($a.avus | .[] | select(.attribute == "study") | .value)"' |\
    sed s"/$(printf '\t')//"g |\
    sed s"/\"//"g |\
    sed s"/____/$(printf '\t')/"g |\
sort | uniq | grep -v -f ${exclude} >> samples.tsv

echo jq search study id done
