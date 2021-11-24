#!/usr/bin/env bats

# Extract the test data

setup() {

    test_dir="${BATS_TEST_DIRNAME}"
    data_dir="${test_dir}/data"
    output_dir="${TMPDIR}"

    array_diff_script="${test_dir}/../differential/diffAtlas_DE_limma.R"
    
    array_diff_exp1="E-GEOD-11166"
    array_diff_exp_dir1="${data_dir}/differential/microarray/experiments/$array_diff_exp1/"
    array_diff_exp_result1="${output_dir}/tmp/${array_diff_exp1}.g1_g2.analytics.tsv"
    
    array_diff_exp2="E-GEOD-10211"
    array_diff_exp_dir2="${data_dir}/differential/microarray/experiments/$array_diff_exp2/"
    array_diff_exp_result2="${output_dir}/tmp/${array_diff_exp2}.g1_g2.analytics.tsv"

    if [ ! -d "${output_dir}/tmp" ]; then
        mkdir -p ${output_dir}/tmp
    fi
    
    export HOME=$output_dir
}

@test "Run a differential gene expression analysis (1 color microarray, multiple contrasts (E-GEOD-11166))" {
    if  [ "$resume" = 'true' ] && [ -f "$array_diff_exp_result1" ]; then
        skip "$array_diff_exp_result exists"
    fi

    run rm -rf $array_diff_exp_result && $array_diff_script $array_diff_exp1 $array_diff_exp_dir1

    [ "$status" -eq 0 ]
    [ -f "$array_diff_exp_result1" ]
}

@test "Run a differential gene expression analysis (1 color microarray, one contrast (E-GEOD-10211))" {
    if  [ "$resume" = 'true' ] && [ -f "$array_diff_exp_result2" ]; then
        skip "$array_diff_exp_result exists"
    fi

    run rm -rf $array_diff_exp_result && $array_diff_script $array_diff_exp2 $array_diff_exp_dir2

    [ "$status" -eq 0 ]
    [ -f "$array_diff_exp_result2" ]
}
