#!/bin/bash    

# This script creates database files (one per organism in Atlas) to be used in the on-the-fly gene set enrichment against differentially expressed genes in each comparisons in a given Atlas organism

scriptDir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source ${scriptDir}/../../bash_util/generic_routines.sh

# quit if not prod user
check_prod_user
if [ $? -ne 0 ]; then
    exit 1
fi

pushd $ATLAS_EXPS
wDir=~/tmp/gse
mkdir -p $wDir
rm -rf $wDir/*
gseFTPDir=/ebi/ftp/pub/databases/arrayexpress/data/atlas/gse/

echo "About to assemble organism-experiments tuples"
grep organism contrastdetails.tsv | grep -P '\torganism\t' | awk -F"\t" '{print $6"\t"$1}' | sort | uniq | sed 's| |_|g' > ${wDir}/exps.aux
awk '$1 = tolower($1)' ${wDir}/exps.aux > ${wDir}/exps.aux.lc && mv ${wDir}/exps.aux.lc ${wDir}/exps.aux
perl -pi -e 's|^oryza_sativa |oryza_sativa_japonica_group |g' ${wDir}/exps.aux
perl -pi -e 's|^oryza_sativa_japonica |oryza_sativa_japonica_group |g' ${wDir}/exps.aux

echo "About to assemble tsv input lists for gsa_prepare_data.R" 
for l in $(cat ${wDir}/exps.aux); do
    o=`echo $l | awk '{print $1}'`
    e=`echo $l | awk '{print $2}'`
    ls $ATLAS_EXPS/$e/${e}*-analytics.tsv > /dev/null
    if [ $? -eq 0 ]; then 
       ls $ATLAS_EXPS/$e/$e-configuration.xml > /dev/null
       if [ $? -eq 0 ]; then 
       	  for f in $(ls $ATLAS_EXPS/$e/${e}*-analytics.tsv); do
    	      echo "$f" >> ${wDir}/${o}.tsvlist.aux
    	  done
    	  echo "$ATLAS_EXPS/$e/$e-configuration.xml" >> ${wDir}/${o}.xmllist.aux
       else
	   echo "WARNING: No $ATLAS_EXPS/$e/$e-configuration.xml present"
       fi
    else
	   echo "WARNING: No $ATLAS_EXPS/$e/${e}*-analytics.tsv present"
    fi
done

echo "About to prepare database files against which enrichment of the user-provided gene sets will be tested"
for o in $(awk '{print $1}' ${wDir}/exps.aux | sort | uniq); do
    echo $o
    $ATLAS_PROD/sw/atlasinstall_prod/atlasprod/irap/atlas_gse/scripts/gsa_prepare_data.R -c 4 -i ${wDir}/${o}.tsvlist.aux -o ${wDir}/${o}.po 1> ${wDir}/${o}.out 2> ${wDir}/${o}.err 
     if [ $? -ne 0 ]; then 
     	echo "Command: '$ATLAS_PROD/sw/atlasinstall_prod/atlasprod/irap/atlas_gse/scripts/gsa_prepare_data.R -c 4 -i ${wDir}/${o}.tsvlist.aux -o ${wDir}/${o}.po' failed"
	break
     fi
done

echo "About to prepare mapping between expAcc-contrastID and contrastTitle"
for l in $(cat $ATLAS_EXPS/contrastdetails.tsv | awk -F"\t" '{print $1"\t"$2}' | sort | uniq); do
    expAcc=`echo $l | awk -F"\t" '{print $1}'`
    contrastID=`echo $l | awk -F"\t" '{print $2}'`
    contrastTitle=`curl -H 'Content-type:application/json' -s "http://wwwdev.ebi.ac.uk/gxa/rest/contrast-summary?experimentAccession=$expAcc&contrastId=$contrastID" | perl -p -e 's|.*\"contrastDescription\":\"(.*?)\".*|$1|g'`
    echo -e "${l}\t$contrastTitle"
done |  perl -p -e "s|\\\\u0027|'|g" >  ${wDir}/contrastTitles.tsv

echo "About to move the files to the Atlas ftp sever"
mv ${wDir}/*.po $gseFTPDir
mv ${wDir}/contrastTitles.tsv $gseFTPDir

popd