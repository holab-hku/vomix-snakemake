import subprocess
import os
import sys
from subprocess import Popen, PIPE, CalledProcessError
import datetime
import json
import yaml
import shutil


class vomix_actions:
    def __init__(self):
        self.name = "vomix"
        self.version = "1.0.0"
        self.description = "vomix is a tool for viral metagenomics analysis."

    def __repr__(self):
        return f"vomix(name={self.name}, version={self.version}, description={self.description})"

    def __str__(self):
        return f"{self.name} v{self.version}: {self.description}"
    
    def get_snakefile(filename):
        # sf = os.path.join(os.path.dirname(os.path.realpath("vomix/workflow/rules/")), filename)
        sf = os.path.realpath("vomix/workflow/rules/" + filename)
        if not os.path.exists(sf):
            sys.exit("Unable to locate the Snakemake file; tried %s" % sf)
        return sf
    
    def env_setup_script() -> str:

        script_path = os.path.realpath("vomix/env_setup.sh")
        print(f"Running script: {script_path}")

        cmd = ['bash', script_path]
        try:
            result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            return result.stdout
        except subprocess.CalledProcessError as e:
            return f"Error: {e.stderr}"
        
    def createScript(self, module, module_obj):
        script = ""

        if module == "preprocess":
            script += 'snakemake --config module="preprocess" decontam-host=' + str(module_obj.decontamHost) + ' outdir="' + module_obj.outdir + '" datadir="' + module_obj.datadir + '" samplelist="' + module_obj.samplelist + '" --use-conda -j 4 --latency-wait 20'
            if module_obj.hasOptions:
                script += " --config"
                if module_obj.dwnldparams:
                    script += ' dwnldparams="' + module_obj.dwnldparams + '"'
                if module_obj.pigzparams:
                    script += ' pigzparams="' + module_obj.pigzparams + '"'
                if module_obj.fastpparams:
                    script += ' fastpparams="' + module_obj.fastpparams + '"'
                if module_obj.hostileparams:
                    script += ' hostileparams="' + module_obj.hostileparams + '"'
                if module_obj.hostilealigner:   
                    script += ' hostilealigner="' + module_obj.hostilealigner + '"'
                if module_obj.alignerparams:    
                    script += ' alignerparams="' + module_obj.alignerparams + '"'
                if module_obj.indexpath:    
                    script += ' indexpath="' + module_obj.indexpath + '"'
            
        elif module == "assembly":
            script += 'snakemake --config module="assembly" assembler=' + module_obj.assembler + ' outdir="' + module_obj.outdir + '" datadir="' + module_obj.datadir + '" samplelist="' + module_obj.samplelist + '" --use-conda -j 4 --latency-wait 20'
            if module_obj.hasOptions:
                script += " --config"
                if module_obj.megahit_minlen:
                    script += ' megahit_minlen=' + str(module_obj.megahit_minlen)
                if module_obj.megahit_params:   
                    script += ' megahit_params="' + module_obj.megahit_params + '"'
                if module_obj.spades_params:
                    script += ' spades_params="' + module_obj.spades_params + '"'
                if module_obj.spades_memory:
                    script += ' spades_memory=' + str(module_obj.spades_memory)

        elif module == "viral-identify":
            script += 'snakemake --config module="viral-identify" ' + 'outdir="' + str(module_obj.outdir) + '" datadir="' + str(module_obj.datadir) + '" samplelist="' + str(module_obj.samplelist)  + '" splits="' + str(module_obj.splits) + '" --use-conda -j 4 --latency-wait 20'
            if module_obj.hasOptions:
                script += " --config"
                if module_obj.contig_minlen:
                    script += ' contig_minlen=' + str(module_obj.contig_minlen)
                if module_obj.genomad_db:
                    script += ' genomad_db="' + module_obj.genomad_db + '"'
                if module_obj.genomad_minlen:
                    script += ' genomad_minlen=' + str(module_obj.genomad_minlen)
                if module_obj.genomad_params:
                    script += ' genomad_params="' + module_obj.genomad_params + '"'
                if module_obj.genomad_cutoff:
                    script += ' genomad_cutoff=' + str(module_obj.genomad_cutoff)
                if module_obj.checkv_original:
                    script += ' checkv_original=' + str(module_obj.checkv_original)
                if module_obj.checkv_params:
                    script += ' checkv_params="' + module_obj.checkv_params + '"'
                if module_obj.checkv_database:
                    script += ' checkv_database="' + module_obj.checkv_database + '"'
                if module_obj.clustering_fast:
                    script += ' clustering_fast=' + str(module_obj.clustering_fast)
                if module_obj.cdhit_params:
                    script += ' cdhit_params="' + module_obj.cdhit_params + '"'
                if module_obj.vOTU_ani:
                    script += ' vOTU_ani=' + str(module_obj.vOTU_ani)
                if module_obj.vOTU_targetcov:
                    script += ' vOTU_targetcov=' + str(module_obj.vOTU_targetcov)
                if module_obj.vOTU_querycov:
                    script += ' vOTU_querycov=' + str(module_obj.vOTU_querycov)
                        
        elif module == "viral-taxonomy":
            script += 'snakemake --config module="viral-taxonomy" ' + 'fasta="' + module_obj.fasta + '" outdir="' + module_obj.outdir + '" --use-conda -j 4 --latency-wait 20'
            if module_obj.hasOptions:
                script += " --config"
                if module_obj.viphogs_hmmeval:
                    script += ' viphogs_hmmeval=' + str(module_obj.viphogs_hmmeval)
                if module_obj.viphogs_prop:
                    script += ' viphogs_prop=' + str(module_obj.viphogs_prop)
                if module_obj.PhaBox2_db:
                    script += ' PhaBox2_db="' + module_obj.PhaBox2_db + '"'
                if module_obj.phagcn_minlen:
                    script += ' phagcn_minlen=' + str(module_obj.phagcn_minlen)
                if module_obj.phagcn_params:
                    script += ' phagcn_params="' + module_obj.phagcn_params + '"'
                if module_obj.diamond_params:
                    script += ' diamond_params="' + module_obj.diamond_params + '"'
                if module_obj.genomad_db:
                    script += ' genomad_db="' + module_obj.genomad_db + '"'
                if module_obj.genomad_params:
                    script += ' genomad_params="' + module_obj.genomad_params + '"'

        elif module == "viral-host":
            script += 'snakemake --config module="viral-host" ' + 'fasta="' + module_obj.fasta + '" outdir="' + module_obj.outdir + '" --use-conda -j 4 --latency-wait 20'
            if module_obj.hasOptions:
                script += " --config"
                if module_obj.CHERRY_params:
                    script += ' CHERRY_params="' + module_obj.CHERRY_params + '"'
                if module_obj.PhaTYP_params:
                    script += ' PhaTYP_params="' + module_obj.PhaTYP_params + '"'
                if module_obj.iphop_cutoff:
                    script += ' iphop_cutoff=' + str(module_obj.iphop_cutoff)
                if module_obj.iphop_params:
                    script += ' iphop_params="' + module_obj.iphop_params + '"'

        elif module == "viral-community":
            script += 'snakemake --config module="viral-community" ' + 'outdir="' + module_obj.outdir + '" datadir="' + module_obj.datadir + '" samplelist="' + module_obj.samplelist + '" --use-conda -j 4 -c 4 --latency-wait 20'
            if module_obj.hasOptions:   
                script += " --config"
                if module_obj.mpa_indexv:
                    script += ' mpa_indexv="' + module_obj.mpa_indexv + '"'
                if module_obj.mpa_params:
                    script += ' mpa_params="' + module_obj.mpa_params + '"'
        
        elif module == "viral-annotate":
            script += 'snakemake --config module="viral-annotate" ' + 'outdir="' + module_obj.outdir + '" datadir="' + module_obj.datadir + '" samplelist="' + module_obj.samplelist + '" --use-conda -j 4 --latency-wait 20'
            if module_obj.hasOptions:
                script += " --config"
                if module_obj.eggNOG_params:
                    script += ' eggNOG_params="' + module_obj.eggNOG_params + '"'
                if module_obj.PhaVIP_params:
                    script += ' PhaVIP_params="' + module_obj.PhaVIP_params + '"'
            
        elif module == "prok-community":
            script += 'snakemake --config module="prok-community" ' + 'outdir="' + module_obj.outdir + '" datadir="' + module_obj.datadir + '" samplelist="' + module_obj.samplelist + '" --use-conda -j 4 -c 4 --latency-wait 20'
            if module_obj.hasOptions:
                script += " --config"
                if module_obj.mpa_indexv:
                    script += ' mpa_indexv="' + module_obj.mpa_indexv + '"'
                if module_obj.mpa_params:
                    script += ' mpa_params="' + module_obj.mpa_params + '"'
        
        elif module == "prok-binning":
            script += f"echo 'Running prok-binning with options: {module_obj}'\n"

        elif module == "prok-annotate":
            script += 'snakemake --config module="prok-annotate" ' + 'outdir="' + module_obj.outdir + '" datadir="' + module_obj.datadir + '" samplelist="' + module_obj.samplelist + '" --use-conda -j 4 --latency-wait 20'

        elif module == "end-to-end":
           script += 'snakemake --config module="end-to-end" ' + 'outdir="' + module_obj.outdir + '" datadir="' + module_obj.datadir + '" samplelist="' + module_obj.samplelist + '" --use-conda -j 4 -c 4 --latency-wait 20'
        
        elif module == "cluster-fast":
            script += 'snakemake --config module="cluster-fast" ' + 'fasta="' + module_obj.fasta + '" outdir="' + module_obj.outdir + '" --use-conda -j 4 --latency-wait 20'
            if module_obj.hasOptions:
                script += " --config"
                if module_obj.clustering_fast:
                    script += ' clustering_fast=' + str(module_obj.clustering_fast)
                if module_obj.cdhit_params:
                    script += ' cdhit_params="' + module_obj.cdhit_params + '"'
                if module_obj.vOTU_ani:
                    script += ' vOTU_ani=' + str(module_obj.vOTU_ani)
                if module_obj.vOTU_targetcov:
                    script += ' vOTU_targetcov=' + str(module_obj.vOTU_targetcov)
                if module_obj.vOTU_querycov:
                    script += ' vOTU_querycov=' + str(module_obj.vOTU_querycov)
        
        elif module == "checkv-pyhmmer":
            script += 'snakemake --config module="checkv-pyhmmer" ' + 'fasta="' + module_obj.fasta + '" outdir="' + module_obj.outdir + '" --use-conda -j 4 --latency-wait 20'
            if module_obj.hasOptions:
                script += " --config"
                if module_obj.checkv_original:
                    script += ' checkv_original=' + str(module_obj.checkv_original)
                if module_obj.checkv_params:
                    script += ' checkv_params="' + module_obj.checkv_params + '"'
                if module_obj.checkv_database:
                    script += ' checkv_database="' + module_obj.checkv_database + '"'
        
        elif module == "setup-database":
            script += 'snakemake --config module="setup-database" ' + 'fasta="' + module_obj.fasta + '" outdir="' + module_obj.outdir + '" --use-conda -j 4 --latency-wait 20'
            if module_obj.hasOptions:
                script += " --config"
                if module_obj.PhaBox2_db:
                    script += ' PhaBox2_db="' + module_obj.PhaBox2_db + '"'
                if module_obj.genomad_db:
                    script += ' genomad_db="' + module_obj.genomad_db + '"'
                if module_obj.checkv_db:
                    script += ' checkv_db="' + module_obj.checkv_db + '"'
                if module_obj.eggNOG_db:
                    script += ' eggNOG_db="' + module_obj.eggNOG_db + '"'
                if module_obj.eggNOG_db_params:
                    script += ' eggNOG_db_params="' + module_obj.eggNOG_db_params + '"'
                if module_obj.virsorter2_db:
                    script += ' virsorter2_db="' + module_obj.virsorter2_db + '"'
                if module_obj.iphop_db:
                    script += ' iphop_db="' + module_obj.iphop_db + '"'
                if module_obj.humann_db:
                    script += ' humann_db="' + module_obj.humann_db + '"'
        
        else:
            script += f"echo 'Unknown module: {module}'\n"

        return script

    def run_module(self, module, module_obj):
        # save_script_path = os.path.realpath("vomix/runModules/" + module +".sh")

        outdir = module_obj.outdir
        datadir = module_obj.datadir

        if not (os.path.exists(outdir) and os.path.exists(os.path.join(outdir, ".vomix"))):
            os.makedirs(os.path.join(outdir, ".vomix"), exist_ok=True)

        ### Save configuration file
        now = datetime.datetime.now()
        latest_run = now.strftime("%Y%m%d_%H%M%S")
        logdir = os.path.join(outdir, ".vomix/log/vomix" + latest_run)

        os.makedirs(logdir, exist_ok=True)
        # TODO update the config file to include latest_run
        # with open(os.path.join(logdir,  "config.json"), "w") as configf:
        #     json.dump(config, configf)

        # Create a new config file with latest_run
        shutil.copy(os.path.realpath("config/config.yml"), logdir)
        with open(logdir + "/config.yml") as f:
            list_doc = yaml.safe_load(f)
            list_doc["latest_run"] = latest_run

        with open(logdir + "/config.yml", "w") as f:
            yaml.dump(list_doc, f)

        script_path = os.path.realpath(logdir + "/snakemake" +".sh")
        print(f"Running script: {script_path}")
        script = self.createScript(module, module_obj)

        # save the script to another location for history
        # with open(save_script_path, "w") as f:   
        #     f.write(script)
        
        # save the script to the main location to run
        with open(script_path, "w") as f:
            f.write(script)

        cmd = ['bash', script_path]

        try:
            # stdout=open('out.log', 'w'), stderr=open('error.log', 'a')
            # result = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
            
            with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
                for line in p.stdout:
                    print(line, end='') 
            if p.returncode != 0:
                raise CalledProcessError(p.returncode, p.args)
            # for stdout_line in iter(result.stdout.readline, ""):
                # yield stdout_line 
            # out = result.stdout.close()
            # return_code = result.wait()
            # if return_code:
            #     raise subprocess.CalledProcessError(return_code, cmd)            
            # return out
        except subprocess.CalledProcessError as e:
            return f"Error: {e.stderr}"
        
    def run_last_module(module) -> str:
        try:
            script_path = os.path.realpath("vomix/snakemake.sh")

            with open(script_path, "r") as file:
                lines = file.readlines()

            if not lines:
                return f"Error: The script file is empty."

            with open(script_path, "w") as file:
                for line in lines:
                    if 'snakemake --config module=' in line:
                        # Replace the module value in the line
                        parts = line.split('module="')
                        if len(parts) > 1:
                            rest = parts[1].split('"', 1)
                            if len(rest) > 1:
                                line = parts[0] + f'module="{module}' + '"' + rest[1][rest[1].find('"'):]
                    file.write(line)

            cmd = ['bash', script_path]

            # result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            result = subprocess.run(cmd, stdout=open('out.log', 'w'), stderr=open('error.log', 'a'),)
            return result.stdout 
        except subprocess.CalledProcessError as e:
            return f"Error: {e.stderr}"
