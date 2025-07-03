class PreProcessingModule:
    # snakemake --config module="preprocess" decontam-host=False outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "preprocess"
    def __init__(self, decontamHost=False, outdir="", datadir="", samplelist=""):
        self.decontamHost = decontamHost
        self.outdir = outdir
        self.datadir = datadir
        self.samplelist = samplelist

class AssemblyCoAssemblyModule:
    # snakemake --config module="assembly" assembler="megahit" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "assembly"
    def __init__(self, assembler=False, outdir="", datadir="", samplelist=""):
        self.assembler = assembler
        self.outdir = outdir
        self.datadir = datadir
        self.samplelist = samplelist
    
class ViralIdentifyModule:
    # snakemake --config module="viral-identify" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "viral-identify"
    def __init__(self, outdir="", datadir="", samplelist=""):
        self.outdir = outdir
        self.datadir = datadir
        self.samplelist = samplelist
    
class ViralTaxonomyModule:
    # snakemake --config module="viral-taxonomy" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
    name = "viral-taxonomy"
    def __init__(self, fasta="", outdir=""):
        self.fasta = fasta
        self.outdir = outdir

class ViralHostModule:
    # snakemake --config module="viral-host" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
    name = "viral-host"
    def __init__(self, fasta="", outdir=""):
        self.fasta = fasta
        self.outdir = outdir

class ViralCommunityModule:
    # snakemake --config module="viral-community" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 -c 4 --latency-wait 20
    name = "viral-community"
    def __init__(self, outdir="", datadir="", samplelist=""):
        self.outdir = outdir
        self.datadir = datadir
        self.samplelist = samplelist

class ViralAnnotateModule:
    # snakemake --config module="viral-annotate" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "viral-annotate"
    def __init__(self, outdir="", datadir="", samplelist=""):
        self.outdir = outdir
        self.datadir = datadir
        self.samplelist = samplelist

class ProkaryoticCommunityModule:
    # snakemake --config module="prok-community" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 -c 4 --latency-wait 20
    name = "prok-community"
    def __init__(self, outdir="", datadir="", samplelist=""):
        self.outdir = outdir
        self.datadir = datadir
        self.samplelist = samplelist

# TBD module
class ProkaryoticBinningModule:
    name = "prok-binning"

class ProkaryoticAnnotateModule:
    # snakemake --config module="prok-annotate" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "prok-annotate"
    def __init__(self, outdir="", datadir="", samplelist=""):
        self.outdir = outdir
        self.datadir = datadir
        self.samplelist = samplelist

class EndToEndModule:
    # snakemake --config module="end-to-end" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 -c 4
    name = "end-to-end"
    def __init__(self, outdir="", datadir="", samplelist=""):
        self.outdir = outdir
        self.datadir = datadir
        self.samplelist = samplelist

class ClusterFastModule:
    # snakemake --config module="cluster-fast" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
    name = "cluster-fast"
    def __init__(self, fasta="", outdir=""):
        self.fasta = fasta
        self.outdir = outdir

class CheckVPyHMMERModule:
    # snakemake --config module="checkv-pyhmmer" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
    name = "checkv-pyhmmer"
    def __init__(self, fasta="", outdir=""):
        self.fasta = fasta
        self.outdir = outdir

class SetupDatabaseModule:
    # snakemake --config module="setup-database" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
    name = "setup-database"
    def __init__(self, fasta="", outdir=""):
        self.fasta = fasta
        self.outdir = outdir
