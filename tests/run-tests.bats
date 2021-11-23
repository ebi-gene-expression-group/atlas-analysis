#!/usr/bin/env bats

# Extract the test data

setup() {
    test_dir="${BATS_TEST_DIRNAME}"
    data_dir="${test_dir}/data"
    output_dir="${test_dir}/outputs"

    array_diff_script="${test_dir}/../differential/diffAtlas_DE_limma.R"
    array_diff_exp="E-GEOD-11166"
    array_diff_exp_dir="${data_dir}/differential/microarray/experiments/E-GEOD-11166/"
    array_diff_exp_result="${output_dir}/tmp/E-GEOD-11166.g1_g2.analytics.tsv"

    if [ ! -d "${output_dir}/tmp" ]; then
        mkdir -p ${output_dir}/tmp
    fi
    
    export HOME=$output_dir
}

@test "Run a differential gene expression analysis (1 color microarray)" {
    if  [ "$resume" = 'true' ] && [ -f "$array_diff_exp_result" ]; then
        skip "$array_diff_exp_result exists"
    fi

    run rm -rf $array_diff_exp_result && $array_diff_script $array_diff_exp $array_diff_exp_dir

    [ "$status" -eq 0 ]
    [ -f "$array_diff_exp_result" ]
}
