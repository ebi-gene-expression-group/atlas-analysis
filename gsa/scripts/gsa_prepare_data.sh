#!/bin/bash    

[ -z ${ATLAS_EXPS+x} ] && echo "Env var ATLAS_EXPS for the path to atlas experiments needs to be defined." && exit 1
[ -z ${ATLAS_GSA_FTP_DIR+x} ] && echo "Env var ATLAS_GSA_FTP_DIR for the path to the atlas ftp gsa folder needs to be defined." && exit 1

WORKDIR="${ATLAS_GSA_WORKDIR:-${TMPDIR}/gsa}"

# This script creates database files (one per organism in Atlas) to be used in the on-the-fly gene set enrichment against differentially expressed genes in each comparisons in a given Atlas organism

scriptDir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# quit if not prod user
atlas-bash-util check_prod_user
if [ $? -ne 0 ]; then
    exit 1
fi

pushd $ATLAS_EXPS > /dev/null
mkdir -p $WORKDIR
rm -rf $WORKDIR/*

echo "About to assemble organism-experiments tuples"
grep organism contrastdetails.tsv | grep -P '\torganism\t' | awk -F"\t" '{print $6"\t"$1}' | sort | uniq | sed 's| |_|g' > ${WORKDIR}/exps.aux
awk '$1 = tolower($1)' ${WORKDIR}/exps.aux > ${WORKDIR}/exps.aux.lc && mv ${WORKDIR}/exps.aux.lc ${WORKDIR}/exps.aux
perl -pi -e 's|^oryza_sativa |oryza_sativa_japonica_group |g' ${WORKDIR}/exps.aux
perl -pi -e 's|^oryza_sativa_japonica |oryza_sativa_japonica_group |g' ${WORKDIR}/exps.aux

echo "About to assemble tsv input lists for gsa_prepare_data.R" 

cat ${WORKDIR}/exps.aux | while read -r l; do
    o=`echo $l | awk '{print $1}'`
    e=`echo $l | awk '{print $2}'`
    ls $ATLAS_EXPS/$e/${e}*-analytics.tsv > /dev/null
    if [ $? -eq 0 ]; then 
       ls $ATLAS_EXPS/$e/$e-configuration.xml > /dev/null
       if [ $? -eq 0 ]; then 
       	  for f in $(ls $ATLAS_EXPS/$e/${e}*-analytics.tsv); do
    	      echo "$f" >> ${WORKDIR}/${o}.tsvlist.aux
    	  done
    	  echo "$ATLAS_EXPS/$e/$e-configuration.xml" >> ${WORKDIR}/${o}.xmllist.aux
       else
	   echo "WARNING: No $ATLAS_EXPS/$e/$e-configuration.xml present"
       fi
    else
	   echo "WARNING: No $ATLAS_EXPS/$e/${e}*-analytics.tsv present"
    fi
done

echo "About to prepare database files against which enrichment of the user-provided gene sets will be tested"
for o in $(awk '{print $1}' ${WORKDIR}/exps.aux | sort | uniq); do
    echo $o
    $scriptDir/gsa_prepare_data.R -c 4 -i ${WORKDIR}/${o}.tsvlist.aux -o ${WORKDIR}/${o}.po 1> ${WORKDIR}/${o}.out 2> ${WORKDIR}/${o}.err
     if [ $? -ne 0 ]; then
        echo "Command: '$scriptDir/gsa_prepare_data.R -c 4 -i ${WORKDIR}/${o}.tsvlist.aux -o ${WORKDIR}/${o}.po' failed"
	      break
     fi
done

echo "About to prepare mapping between expAcc-contrastID and contrastTitle"
for l in $(cat $ATLAS_EXPS/contrastdetails.tsv | awk -F"\t" '{print $1"\t"$2}' | sort | uniq); do
    expAcc=`echo $l | awk -F"\t" '{print $1}'`
    contrastID=`echo $l | awk -F"\t" '{print $2}'`
    contrastTitle=`curl -H 'Content-type:application/json' -s "http://wwwdev.ebi.ac.uk/gxa/rest/contrast-summary?experimentAccession=$expAcc&contrastId=$contrastID" | perl -p -e 's|.*\"contrastDescription\":\"(.*?)\".*|$1|g'`
    echo -e "${l}\t$contrastTitle"
done |  perl -p -e "s|\\\\u0027|'|g" >  ${WORKDIR}/contrastTitles.tsv

echo "About to move the files to the Atlas ftp sever"
mkdir -p $ATLAS_GSA_FTP_DIR
mv ${WORKDIR}/*.po $ATLAS_GSA_FTP_DIR
mv ${WORKDIR}/contrastTitles.tsv $ATLAS_GSA_FTP_DIR

# Clean up

rm -rf $WORKDIR

popd > /dev/null
