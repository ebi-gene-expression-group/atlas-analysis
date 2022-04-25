#!/usr/bin/env bats

# Extract the test data

setup() {

    test_dir="${BATS_TEST_DIRNAME}"
    data_dir="${test_dir}/data"
    output_dir="/tmp"

    array_diff_script="${test_dir}/../differential/diffAtlas_DE_limma.R"
    rnaseq_diff_script="${test_dir}/../differential/diffAtlas_DE_deseq.R"
    gsa_prepare_script="${test_dir}/../gsa/scripts/gsa_prepare_data.sh"   
 
    array_diff_exp1="E-GEOD-11166"
    array_diff_exp_dir1="${data_dir}/differential/microarray/experiments/$array_diff_exp1/"
    array_diff_exp_result1="${output_dir}/tmp/${array_diff_exp1}.g1_g2.analytics.tsv"
    
    array_diff_exp2="E-GEOD-10211"
    array_diff_exp_dir2="${data_dir}/differential/microarray/experiments/$array_diff_exp2/"
    array_diff_exp_result2="${output_dir}/tmp/${array_diff_exp2}.g1_g2.analytics.tsv"

    array_diff_exp3="E-GEOD-20333"
    array_diff_exp_dir3="${data_dir}/differential/microarray/experiments/$array_diff_exp3/"
    array_diff_exp_result3="${output_dir}/tmp/${array_diff_exp3}.g2_g1.analytics.tsv"
    echo "$array_diff_exp_result3"

    array_diff_exp4="E-GEOD-44048"
    array_diff_exp_dir4="${data_dir}/differential/microarray/experiments/$array_diff_exp4/"
    array_diff_exp_result4="${output_dir}/tmp/${array_diff_exp4}.g1_g2.analytics.tsv"

    rnaseq_diff_exp1="E-MTAB-7471"
    rnaseq_diff_exp_dir1="${data_dir}/differential/rnaseq/experiments/$rnaseq_diff_exp1/"
    rnaseq_diff_exp_result1="${output_dir}/tmp/${rnaseq_diff_exp1}.g1_g2.analytics.tsv"

    rnaseq_diff_exp2="E-MTAB-8510"
    rnaseq_diff_exp_dir2="${data_dir}/differential/rnaseq/experiments/$rnaseq_diff_exp2/"
    rnaseq_diff_exp_result2="${output_dir}/tmp/${rnaseq_diff_exp2}.g2_g1.analytics.tsv"
    rnaseq_diff_matrix2="${data_dir}/differential/rnaseq/experiments/${rnaseq_diff_exp2}/${rnaseq_diff_exp2}-raw-counts.tsv.undecorated "
    
    rnaseq_diff_exp3="E-GEOD-52687"
    rnaseq_diff_exp_dir3="${data_dir}/differential/rnaseq/experiments/$rnaseq_diff_exp3/"
    rnaseq_diff_exp_result3="${output_dir}/tmp/${rnaseq_diff_exp3}.g3_g2.analytics.tsv"
    rnaseq_diff_matrix3="${data_dir}/differential/rnaseq/experiments/${rnaseq_diff_exp3}/${rnaseq_diff_exp3}-raw-counts.tsv.undecorated "

    export ATLAS_GSA_WORKDIR="${output_dir}/gsa"
    export ATLAS_GSA_FTP_DIR="${ATLAS_GSA_WORKDIR}_ftp"
    gsa_prep_contrast_titles="$ATLAS_GSA_FTP_DIR/contrastTitles.tsv"
    gsa_ss_index="${ATLAS_GSA_FTP_DIR}/saccharomyces_cerevisiae.po"
    gsa_hs_index="${ATLAS_GSA_FTP_DIR}/homo_sapiens.po"
	
    if [ ! -d "${output_dir}/tmp" ]; then
        mkdir -p ${output_dir}/tmp
    fi
    
    export HOME=$output_dir
    export ATLAS_EXPS="${data_dir}/differential/rnaseq/experiments"
    export ATLAS_PROD_USER="$(whoami)"
}

@test "Run a differential gene expression analysis (1 color microarray, multiple contrasts (E-GEOD-11166))" {
    if  [ "$resume" = 'true' ] && [ -f "$array_diff_exp_result1" ]; then
        skip "$array_diff_exp_result1 exists"
    fi

    run rm -rf $array_diff_exp_result1 && $array_diff_script $array_diff_exp1 $array_diff_exp_dir1

    [ "$status" -eq 0 ]
    [ -f "$array_diff_exp_result1" ]
}

@test "Run a differential gene expression analysis (1 color microarray, one contrast (E-GEOD-10211))" {
    if  [ "$resume" = 'true' ] && [ -f "$array_diff_exp_result2" ]; then
        skip "$array_diff_exp_result2 exists"
    fi

    run rm -rf $array_diff_exp_result2 && $array_diff_script $array_diff_exp2 $array_diff_exp_dir2

    [ "$status" -eq 0 ]
    [ -f "$array_diff_exp_result2" ]
}

@test "Run a differential gene expression analysis (1 color microarray, one contrast, one batch effect (E-GEOD-20333))" {
    if  [ "$resume" = 'true' ] && [ -f "$array_diff_exp_result3" ]; then
        skip "$array_diff_exp_result3 exists"
    fi

    run rm -rf $array_diff_exp_result3 && $array_diff_script $array_diff_exp3 $array_diff_exp_dir3

    [ "$status" -eq 0 ]
    [ -f "$array_diff_exp_result3" ]
}

@test "Run a differential gene expression analysis (2 color microarray, one contrast (E-GEOD-44048))" {
    if  [ "$resume" = 'true' ] && [ -f "$array_diff_exp_result4" ]; then
        skip "$array_diff_exp_result4 exists"
    fi

    run rm -rf $array_diff_exp_result4 && $array_diff_script $array_diff_exp4 $array_diff_exp_dir4

    [ "$status" -eq 0 ]
    [ -f "$array_diff_exp_result4" ]
}

@test "Run a differential gene expression analysis (RNA-seq, one contrast (E-MTAB-7471))" {
    if  [ "$resume" = 'true' ] && [ -f "$rnaseq_diff_exp_result1" ]; then
        skip "$rnaseq_diff_exp_result1 exists"
    fi

    run rm -rf $rnaseq_diff_exp_result1 && $rnaseq_diff_script $rnaseq_diff_exp1 $rnaseq_diff_exp_dir1

    [ "$status" -eq 0 ]
    [ -f "$rnaseq_diff_exp_result1" ]
}

@test "Run a differential gene expression analysis (RNA-seq, one contrast, one batch effect (E-MTAB-8510))" {
    if  [ "$resume" = 'true' ] && [ -f "$rnaseq_diff_exp_result2" ]; then
        skip "$rnaseq_diff_exp_result2 exists"
    fi

    run rm -rf $rnaseq_diff_exp_result2 && $rnaseq_diff_script $rnaseq_diff_exp2 $rnaseq_diff_exp_dir2

    [ "$status" -eq 0 ]
    [ -f "$rnaseq_diff_exp_result2" ]
}

@test "Run a differential gene expression analysis (RNA-seq, multiple contrasts (E-GEOD-52687))" {
    if  [ "$resume" = 'true' ] && [ -f "$rnaseq_diff_exp_result3" ]; then
        skip "$rnaseq_diff_exp_result3 exists"
    fi

    run rm -rf $rnaseq_diff_exp_result3 && $rnaseq_diff_script $rnaseq_diff_exp3 $rnaseq_diff_exp_dir3

    [ "$status" -eq 0 ]
    [ -f "$rnaseq_diff_exp_result3" ]
}

@test "Run gene set analysis preparation" {
    if  [ "$resume" = 'true' ] && [ -f "$gsa_prep_contrast_titles" ] && [ -f "$gsa_ss_index" ] && [ -f "$gsa_hs_index" ]; then
        skip "$gsa_prep_contrast_titles and $gsa_ss_index and $gsa_hs_index exist"
    fi

    run rm -rf ${ATLAS_GSA_WORKDIR}_ftp && $gsa_prepare_script

    [ "$status" -eq 0 ]
    [ -f "$gsa_prep_contrast_titles" ]	
    [ -f "$gsa_ss_index" ]	
    [ -f "$gsa_hs_index" ]	
}

