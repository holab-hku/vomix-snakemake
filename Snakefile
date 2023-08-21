rule virsorter2:
    input:
        "data/testdata/test.fa"
    output:
        "output/intermediate/virsorter2"
    conda:
        "envs/"