import subprocess
import os
import sys

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
            script += 'nohup snakemake --config module="preprocess" decontam-host=' + str(module_obj.decontamHost) + ' outdir="' + module_obj.outdir + '" datadir="' + module_obj.datadir + '" samplelist="' + module_obj.samplelist + '" --use-conda -j 4 --latency-wait 20'
            
        elif module == "assembly":
            script += 'nohup snakemake --config module="assembly" assembler=' + module_obj.assembler + ' outdir="' + module_obj.outdir + '" datadir="' + module_obj.datadir + '" samplelist="' + module_obj.samplelist + '" --use-conda -j 4 --latency-wait 20'

        elif module == "viral-identify":
            script += 'nohup snakemake --config module="viral-identify" ' + 'outdir="' + module_obj.outdir + '" datadir="' + module_obj.datadir + '" samplelist="' + module_obj.samplelist + '" --use-conda -j 4 --latency-wait 20'
                        
        elif module == "viral-taxonomy":
            script += 'nohup snakemake --config module="viral-taxonomy" ' + 'fasta="' + module_obj.fasta + '" outdir="' + module_obj.outdir + '" --use-conda -j 4 --latency-wait 20'

        elif module == "viral-host":
            script += 'nohup snakemake --config module="viral-host" ' + 'fasta="' + module_obj.fasta + '" outdir="' + module_obj.outdir + '" --use-conda -j 4 --latency-wait 20'

        elif module == "viral-community":
            script += 'nohup snakemake --config module="viral-community" ' + 'outdir="' + module_obj.outdir + '" datadir="' + module_obj.datadir + '" samplelist="' + module_obj.samplelist + '" --use-conda -j 4 -c 4 --latency-wait 20'
        
        elif module == "viral-annotate":
            script += 'nohup snakemake --config module="viral-annotate" ' + 'outdir="' + module_obj.outdir + '" datadir="' + module_obj.datadir + '" samplelist="' + module_obj.samplelist + '" --use-conda -j 4 --latency-wait 20'
            
        elif module == "prok-community":
            script += 'nohup snakemake --config module="prok-community" ' + 'outdir="' + module_obj.outdir + '" datadir="' + module_obj.datadir + '" samplelist="' + module_obj.samplelist + '" --use-conda -j 4 -c 4 --latency-wait 20'
        
        elif module == "prok-binning":
            script += f"echo 'Running prok-binning with options: {module_obj}'\n"

        elif module == "prok-annotate":
            script += 'nohup snakemake --config module="prok-annotate" ' + 'outdir="' + module_obj.outdir + '" datadir="' + module_obj.datadir + '" samplelist="' + module_obj.samplelist + '" --use-conda -j 4 --latency-wait 20'
            
        elif module == "end-to-end":
           script += 'nohup snakemake --config module="end-to-end" ' + 'outdir="' + module_obj.outdir + '" datadir="' + module_obj.datadir + '" samplelist="' + module_obj.samplelist + '" --use-conda -j 4 -c 4 --latency-wait 20'
        
        elif module == "cluster-fast":
            script += 'nohup snakemake --config module="cluster-fast" ' + 'fasta="' + module_obj.fasta + '" outdir="' + module_obj.outdir + '" --use-conda -j 4 --latency-wait 20'
        
        elif module == "checkv-pyhmmer":
            script += 'nohup snakemake --config module="checkv-pyhmmer" ' + 'fasta="' + module_obj.fasta + '" outdir="' + module_obj.outdir + '" --use-conda -j 4 --latency-wait 20'
        
        elif module == "setup-database":
            script += 'nohup snakemake --config module="setup-database" ' + 'fasta="' + module_obj.fasta + '" outdir="' + module_obj.outdir + '" --use-conda -j 4 --latency-wait 20'
        
        else:
            script += f"echo 'Unknown module: {module}'\n"

        return script

    def run_module(self, module, module_obj) -> str:

        script_path = os.path.realpath("vomix/" + "snakemake" +".sh")
        save_script_path = os.path.realpath("vomix/runModules/" + module +".sh")

        script = self.createScript(module, module_obj)

        # save the script to another location for history
        with open(save_script_path, "w") as f:   
            f.write(script)
        
        # save the script to the main location to run
        with open(script_path, "w") as f:

            f.write(script)

            # f.write('nohup snakemake --config module="' + module 
            #         + '" outdir="' + module_options["outdir"] 
            #         + '" datadir="' + module_options['datadir'] 
            #         + '" assembler="' + module_options['assembler'] 
            #         + '" samplelist="' + module_options['samplelist'] 
            #         + '" decontam-host=' + str(module_options['decontam_host'])
            #         + ' binning-consensus=' + str(module_options['binning_consensus'])
            #         + ' --use-conda --rerun-incomplete --rerun-triggers mtime --latency-wait 20 --retries 0 --sdm conda -j 88 --executor  cluster-generic --cluster-generic-submit-cmd "qsub -N {log} -l nodes=1:ppn={threads} -l mem={resources.mem_mb}m -l walltime=120:00:00 '
            #         + '-M ' + module_options['email']
            #         + ' -q cgsd -o qsub.log -e qsub.log -m a"  ' )
            #         # + ' > nohup.' + module + '.out &') 

        cmd = ['bash', script_path]

        try:
            result = subprocess.run(cmd, stdout=open('out.log', 'w'), stderr=open('error.log', 'a'),)
            return result.stdout
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


        
