class PreProcessingModule:
    # snakemake --config module="preprocess" decontam-host=False outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "preprocess"
    def __init__(self, decontamHost=False, outdir="", datadir="", samplelist="", hasOptions=False, dwnldparams="", pigzparams="", fastpparams="", hostileparams="", hostilealigner="minimap2", alignerparams="-x sr", indexpath="./workflow/database/hostile/human-t2t-hla.fa.gz"):
        self.decontamHost = decontamHost
        self.outdir = outdir
        self.datadir = datadir
        self.samplelist = samplelist
        self.hasOptions = hasOptions
        self.dwnldparams = dwnldparams
        self.pigzparams = pigzparams
        self.fastpparams = fastpparams
        self.hostileparams = hostileparams 
        self.hostilealigner = hostilealigner
        self.alignerparams = alignerparams
        self.indexpath = indexpath

class AssemblyCoAssemblyModule:
    # snakemake --config module="assembly" assembler="megahit" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "assembly"
    def __init__(self, assembler=False, outdir="", datadir="", samplelist="", hasOptions=False, megahit_minlen=300, megahit_params="--prune-level 3", spades_params="--meta", spades_memory=250):
        self.assembler = assembler
        self.outdir = outdir
        self.datadir = datadir
        self.samplelist = samplelist
        self.hasOptions = hasOptions
        self.megahit_minlen = megahit_minlen
        self.megahit_params = megahit_params
        self.spades_params = spades_params
        self.spades_memory = spades_memory
    
class ViralIdentifyModule:
    # snakemake --config module="viral-identify" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "viral-identify"
    def __init__(self, outdir="", datadir="", samplelist="", splits=0, hasOptions=False, contig_minlen=0, genomad_db="workflow/database/genomad", genomad_minlen=1500, genomad_params="", genomad_cutoff=0.7, checkv_original=False, checkv_params="", checkv_database="workflow/database/checkv", clustering_fast=True, cdhit_params="-c 0.95 -aS 0.85 -d 400 -M 0 -n 5", vOTU_ani=95, vOTU_targetcov=85, vOTU_querycov=0):
        self.outdir = outdir
        self.datadir = datadir
        self.samplelist = samplelist
        self.splits = splits
        self.hasOptions = hasOptions
        self.contig_minlen = contig_minlen
        self.genomad_db = genomad_db
        self.genomad_minlen = genomad_minlen
        self.genomad_params = genomad_params
        self.genomad_cutoff = genomad_cutoff
        self.checkv_original = checkv_original
        self.checkv_params = checkv_params
        self.checkv_database = checkv_database
        self.clustering_fast = clustering_fast
        self.cdhit_params = cdhit_params
        self.vOTU_ani = vOTU_ani
        self.vOTU_targetcov = vOTU_targetcov
        self.vOTU_querycov = vOTU_querycov
    
class ViralTaxonomyModule:
    # snakemake --config module="viral-taxonomy" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
    name = "viral-taxonomy"
    def __init__(self, fasta="", outdir="", hasOptions=False, viphogs_hmmeval=0.01, viphogs_prop=0.06, PhaBox2_db="workflow/database/phabox_db_v2", phagcn_minlen=1500, phagcn_params="", diamond_params="--query-cover 50 --subject-cover 50 --evalue 1e-5 --max-target-seqs 1000", genomad_db="workflow/database/genomad", genomad_params="--enable-score-calibration --relaxed"):
        self.fasta = fasta
        self.outdir = outdir
        self.hasOptions = hasOptions
        self.viphogs_hmmeval = viphogs_hmmeval
        self.viphogs_prop = viphogs_prop
        self.PhaBox2_db = PhaBox2_db
        self.phagcn_minlen = phagcn_minlen
        self.phagcn_params = phagcn_params
        self.diamond_params = diamond_params
        self.genomad_db = genomad_db
        self.genomad_params = genomad_params

class ViralHostModule:
    # snakemake --config module="viral-host" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
    name = "viral-host"
    def __init__(self, fasta="", outdir="", hasOptions=False, CHERRY_params="", PhaTYP_params="", iphop_cutoff=90, iphop_params=""):
        self.fasta = fasta
        self.outdir = outdir
        self.hasOptions = hasOptions
        self.CHERRY_params = CHERRY_params
        self.PhaTYP_params = PhaTYP_params
        self.iphop_cutoff = iphop_cutoff
        self.iphop_params = iphop_params

class ViralCommunityModule:
    # snakemake --config module="viral-community" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 -c 4 --latency-wait 20
    name = "viral-community"
    def __init__(self, outdir="", datadir="", samplelist="", hasOptions=False, mpa_indexv="mpa_vOct22_CHOCOPhlAnSGB_202212", mpa_params="--ignore_eukaryotes"):
        self.outdir = outdir
        self.datadir = datadir
        self.samplelist = samplelist
        self.hasOptions = hasOptions
        self.mpa_indexv = mpa_indexv
        self.mpa_params = mpa_params    


class ViralAnnotateModule:
    # snakemake --config module="viral-annotate" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "viral-annotate"
    def __init__(self, outdir="", datadir="", samplelist="", hasOptions=False, eggNOG_params="", PhaVIP_params=""):
        self.outdir = outdir
        self.datadir = datadir
        self.samplelist = samplelist
        self.hasOptions = hasOptions
        self.eggNOG_params = eggNOG_params
        self.PhaVIP_params = PhaVIP_params

class ProkaryoticCommunityModule:
    # snakemake --config module="prok-community" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 -c 4 --latency-wait 20
    name = "prok-community"
    def __init__(self, outdir="", datadir="", samplelist="", hasOptions=False, mpa_params="--ignore_eukaryotes", mpa_indexv="mpa_vOct22_CHOCOPhlAnSGB_202212"):
        self.outdir = outdir
        self.datadir = datadir
        self.samplelist = samplelist
        self.hasOptions = hasOptions
        self.mpa_params = mpa_params
        self.mpa_indexv = mpa_indexv

# TBD module
class ProkaryoticBinningModule:
    name = "prok-binning"

class ProkaryoticAnnotateModule:
    # snakemake --config module="prok-annotate" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "prok-annotate"
    def __init__(self, outdir="", datadir="", samplelist="", hasOptions=False):
        self.outdir = outdir
        self.datadir = datadir
        self.samplelist = samplelist
        self.hasOptions = hasOptions

class EndToEndModule:
    # snakemake --config module="end-to-end" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 -c 4
    name = "end-to-end"
    def __init__(self, outdir="", datadir="", samplelist="", hasOptions=False):
        self.outdir = outdir
        self.datadir = datadir
        self.samplelist = samplelist
        self.hasOptions = hasOptions

class ClusterFastModule:
    # snakemake --config module="cluster-fast" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
    name = "cluster-fast"
    def __init__(self, fasta="", outdir="", hasOptions=False, clustering_fast=True, cdhit_params="", vOTU_ani=95, vOTU_targetcov=85, vOTU_querycov=0):
        self.fasta = fasta
        self.outdir = outdir
        self.hasOptions = hasOptions
        self.clustering_fast = clustering_fast
        self.cdhit_params = cdhit_params
        self.vOTU_ani = vOTU_ani
        self.vOTU_targetcov = vOTU_targetcov
        self.vOTU_querycov = vOTU_querycov

class CheckVPyHMMERModule:
    # snakemake --config module="checkv-pyhmmer" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
    name = "checkv-pyhmmer"
    def __init__(self, fasta="", outdir="", hasOptions=False, checkv_original=False, checkv_params="", checkv_database=""):
        self.fasta = fasta
        self.outdir = outdir
        self.hasOptions = hasOptions
        self.checkv_original = checkv_original
        self.checkv_params = checkv_params
        self.checkv_database = checkv_database

class SetupDatabaseModule:
    # snakemake --config module="setup-database" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
    name = "setup-database"
    def __init__(self, fasta="", outdir="", hasOptions=False, PhaBox2_db="", genomad_db="", checkv_db="", eggNOG_db="", eggNOG_db_params="", virsorter2_db="", iphop_db="", humann_db=""):
        self.fasta = fasta
        self.outdir = outdir
        self.hasOptions = hasOptions
        self.PhaBox2_db = PhaBox2_db
        self.genomad_db = genomad_db
        self.checkv_db = checkv_db
        self.eggNOG_db = eggNOG_db
        self.eggNOG_db_params = eggNOG_db_params
        self.virsorter2_db = virsorter2_db
        self.iphop_db = iphop_db
        self.humann_db = humann_db
