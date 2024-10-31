datadir=config["datadir"]

rule validate_samples:
  name: "preprocessing.py validating samples"
  localrule: True
  input:
    config["samplelist"]
  output:
    samplejson=relpath(".vomix/samples.json"),
    assemblyjson=relpath(".vomix/assembly.json")
  params:
    script="workflow/scripts/utility/parse_sample_list.py",
    outdir=relpath(".vomix")
  log: os.path.join(datadir, "log/validate.log")
  conda: "../envs/preprocessing.yml"
  threads: 1
  shell:
    """
    python {params.script} \
        --sample_list {input} \
        --datadir {params.outdir} &> {log}
    """
