rule phabox_db:
  output:
    "workflow/database/vs2/phabox/done.log"
  params:
    output_dir = "workflow/database/", 
    fileid = "1hjACPsIOqqcS5emGaduYvYrCzrIpt2_9", 
    html = 'curl -c ./cookie -s -L "https://drive.google.com/uc?export=download&id={params.fileid}"',
    tmp_dir = "$TMPDIR/phabox/"
  log: "logs/db_phabox.log"
  shell:
    """
    mkdir -p {params.tmp_dir} {params.output_dir}

    curl -Lb ./cookie "https://drive.google.com/uc?export=download&`echo {params.html}|grep -Po '(confirm=[a-zA-Z0-9\-_]+)'`&id={params.fileid}" -o {params.tmp_dir}/tmp.zip 2> {log}

    unzip {params.tmp_dir}/phabox.zip > {output} 2> {log}
    touch {output} 2> {log}

    rm -r {params.tmp_dir}
    """
 
