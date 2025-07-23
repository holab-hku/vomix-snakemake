import click
import sys
import logging
import os
from vomix.vomix_actions import vomix_actions
from vomix.modules import PreProcessingModule, AssemblyCoAssemblyModule, ViralIdentifyModule, ViralTaxonomyModule, ViralHostModule, ViralCommunityModule, ViralAnnotateModule, ProkaryoticCommunityModule, ProkaryoticBinningModule, ProkaryoticAnnotateModule, EndToEndModule, ClusterFastModule, CheckVPyHMMERModule, SetupDatabaseModule

logging.basicConfig(level=logging.INFO)

modules_list = ['assembly','checkv-pyhmmer','checkv'
,'clustering-fast','clustering-sensitive','host-cherry','host','preprocess-decontam','preprocess','prok-annotate','prok-binning','prok-community','refilter-genomad'
,'setup-database','symlink','viral-annotate','viral-benchmark','viral-binning','viral-community','viral-host','viral-identify','viral-refilter','viral-taxonomy']

def useLastOptionsCheck(ctx, param, value):
    i = 2
    totalParams = len(ctx.command.params)
    if value == True:
        while i < totalParams:
            ctx.command.params[i].prompt_required = False
            ctx.command.params[i].required = False
            i += 1

@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
def cli():
    """
    vOMIX-MEGA - a reproducible, scalable, and fast viral metagenomic pipeline
    """

@cli.command(
        'activate',
        context_settings=dict(ignore_unknown_options=True),
        short_help='Set up vomix conda environment'
    )
def activate():
    output = vomix_actions.env_setup_script()
    logging.info(f"Output: {output}")
    message = "Completed vomix environment creation"
    logging.info(message)

@cli.command(
    'preprocess',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Pre-processing module'
)
@click.option('--decontam-host', default=True, required=True)
@click.option('--outdir', required=True, default=None)
@click.option('--datadir', required=False, default=None)
@click.option('--samplelist', required=False, default=None)
@click.option('--dwnldparams', required=False, default=None)
@click.option('--pigzparams', required=False, default=None)
@click.option('--fastpparams', required=False, default=None)
@click.option('--hostileparams', required=False, default=None)
@click.option('--hostilealigner', required=False, default=None)
@click.option('--alignerparams', required=False, default=None)
@click.option('--indexpath', required=False, default=None)
def run_preprocess(decontam_host, outdir, datadir, samplelist, dwnldparams, pigzparams, fastpparams, hostileparams, hostilealigner, alignerparams, indexpath):
        logging.info(f"Running module: preprocess")
        logging.info(f"decontamHost: {decontam_host}, outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = PreProcessingModule()
        module_obj.name = "preprocess"
        # Set the attributes of the module object
        module_obj.decontamHost = decontam_host
        module_obj.outdir = outdir 
        module_obj.datadir = datadir
        module_obj.samplelist = samplelist

        # optional params
        if dwnldparams:
            module_obj.dwnldparams = dwnldparams
            module_obj.hasOptions = True
        if pigzparams:
            module_obj.pigzparams = pigzparams  
            module_obj.hasOptions = True
        if fastpparams:
            module_obj.fastpparams = fastpparams 
            module_obj.hasOptions = True
        if hostileparams:
            module_obj.hostileparams = hostileparams
            module_obj.hasOptions = True
        if hostilealigner:
            module_obj.hostilealigner = hostilealigner
            module_obj.hasOptions = True
        if alignerparams:
            module_obj.alignerparams = alignerparams
            module_obj.hasOptions = True
        if indexpath:
            module_obj.indexpath = indexpath
            module_obj.hasOptions = True

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("preprocess", module_obj)
        logging.info(f"End module run")

@cli.command(
    'assembly',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Assembly & Co-assembly module'
)
@click.option('--assembler', default="megahit", 
required=True)
@click.option('--outdir', required=True, default=None)
@click.option('--datadir', required=False, default=None)
@click.option('--samplelist', required=False, default=None)
@click.option('--megahit-minlen', required=False, default=None)
@click.option('--megahit-params', required=False, default=None)
@click.option('--spades-params', required=False, default=None)
@click.option('--spades-memory', required=False, default=None)
def run_assembly(assembler, outdir, datadir, samplelist, megahit_minlen, megahit_params, spades_params, spades_memory):
        logging.info(f"Running module: assembly")
        logging.info(f"assembler: {assembler}, outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = AssemblyCoAssemblyModule()
        module_obj.name = "assembly"
        # Set the attributes of the module object
        module_obj.assembler = assembler
        module_obj.outdir = outdir 
        module_obj.datadir = datadir
        module_obj.samplelist = samplelist

        if megahit_minlen:
            module_obj.megahit_minlen = megahit_minlen
            module_obj.hasOptions = True
        if megahit_params: 
            module_obj.megahit_params = megahit_params
            module_obj.hasOptions = True
        if spades_params:
            module_obj.spades_params = spades_params
            module_obj.hasOptions = True
        if spades_memory:
            module_obj.spades_memory = spades_memory
            module_obj.hasOptions = True

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("assembler", module_obj)
        logging.info(f"End module run")

@cli.command(
    'viral-identify',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Viral Identify module'
)
@click.option('--outdir', required=True, default=None)
@click.option('--datadir', required=False, default=None)
@click.option('--fasta', required=False, default=None)
@click.option('--samplelist', required=False, default=None)
@click.option('--splits', required=False, default=0)
@click.option('--contig-minlen', required=False, default=None)
@click.option('--genomad-db', required=False, default=None)
@click.option('--genomad-minlen', required=False, default=None)
@click.option('--genomad-params', required=False, default=None)
@click.option('--genomad-cutoff', required=False, default=None)
@click.option('--checkv-original', required=False, default=None)
@click.option('--checkv-params', required=False, default=None)
@click.option('--checkv-database', required=False, default=None)
@click.option('--clustering-fast', required=False, default=None)
@click.option('--cdhit-params', required=False, default=None)
@click.option('--vOTU-ani', required=False, default=None)
@click.option('--vOTU-targetcov', required=False, default=None)
@click.option('--vOTU-querycov', required=False, default=None)
def run_viral_identify(outdir, datadir, fasta, samplelist, splits, contig_minlen, genomad_db, genomad_minlen, genomad_params, genomad_cutoff, checkv_original, checkv_params, checkv_database, clustering_fast, cdhit_params, votu_ani, votu_targetcov, votu_querycov):
        logging.info(f"Running module: viral-identify")
        logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = ViralIdentifyModule()
        module_obj.name = "viral-identify"
        # Set the attributes of the module object
        module_obj.outdir = outdir 
        module_obj.datadir = datadir
        module_obj.samplelist = samplelist
        module_obj.splits = splits
        module_obj.fasta = fasta

        if contig_minlen:
            module_obj.contig_minlen = contig_minlen
            module_obj.hasOptions = True
        if genomad_db:  
            module_obj.genomad_db = genomad_db
            module_obj.hasOptions = True
        if genomad_minlen:
            module_obj.genomad_minlen = genomad_minlen
            module_obj.hasOptions = True
        if genomad_params:  
            module_obj.genomad_params = genomad_params
            module_obj.hasOptions = True    
        if genomad_cutoff:
            module_obj.genomad_cutoff = genomad_cutoff
            module_obj.hasOptions = True
        if checkv_original:
            module_obj.checkv_original = checkv_original
            module_obj.hasOptions = True
        if checkv_params:
            module_obj.checkv_params = checkv_params
            module_obj.hasOptions = True
        if checkv_database: 
            module_obj.checkv_database = checkv_database
            module_obj.hasOptions = True
        if clustering_fast:
            module_obj.clustering_fast = clustering_fast
            module_obj.hasOptions = True
        if cdhit_params:
            module_obj.cdhit_params = cdhit_params
            module_obj.hasOptions = True
        if votu_ani:
            module_obj.vOTU_ani = votu_ani
            module_obj.hasOptions = True
        if votu_targetcov:
            module_obj.vOTU_targetcov = votu_targetcov
            module_obj.hasOptions = True
        if votu_querycov:
            module_obj.vOTU_querycov = votu_querycov
            module_obj.hasOptions = True

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("viral-identify", module_obj)
        logging.info(f"End module run")


@cli.command(
    'viral-taxonomy',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Viral Taxonomy module'
)
@click.option('--fasta', required=True, default=None)
@click.option('--outdir', required=True, default=None)
@click.option('--viphogs-hmmeval', required=False, default=None)
@click.option('--viphogs-prop', required=False, default=None)
@click.option('--PhaBox2-db', required=False, default=None)
@click.option('--phagcn-minlen', required=False, default=None)
@click.option('--phagcn-params', required=False, default=None)
@click.option('--diamond-params', required=False, default=None)
@click.option('--genomad-db', required=False, default=None)
@click.option('--genomad-params', required=False, default=None)
def run_viral_taxonomy(fasta, outdir, viphogs_hmmeval, viphogs_prop, PhaBox2_db, phagcn_minlen, phagcn_params, diamond_params, genomad_db, genomad_params):
        logging.info(f"Running module: viral-taxonomy")
        logging.info(f"fasta: {fasta}, outdir: {outdir}")
        
        module_obj = ViralTaxonomyModule()
        module_obj.name = "viral-taxonomy"
        # Set the attributes of the module object
        module_obj.fasta = fasta 
        module_obj.outdir = outdir 

        if viphogs_hmmeval:
            module_obj.viphogs_hmmeval = viphogs_hmmeval
            module_obj.hasOptions = True
        if viphogs_prop:
            module_obj.viphogs_prop = viphogs_prop
            module_obj.hasOptions = True
        if PhaBox2_db: 
            module_obj.PhaBox2_db = PhaBox2_db
            module_obj.hasOptions = True
        if phagcn_minlen:
            module_obj.phagcn_minlen = phagcn_minlen
            module_obj.hasOptions = True
        if phagcn_params:
            module_obj.phagcn_params = phagcn_params
            module_obj.hasOptions = True
        if diamond_params:  
            module_obj.diamond_params = diamond_params
            module_obj.hasOptions = True
        if genomad_db:
            module_obj.genomad_db = genomad_db
            module_obj.hasOptions = True
        if genomad_params:
            module_obj.genomad_params = genomad_params
            module_obj.hasOptions = True

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("viral-taxonomy", module_obj)
        logging.info(f"End module run")

@cli.command(
    'viral-host',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Viral Host module'
)
@click.option('--fasta', required=True, default=None)
@click.option('--outdir', required=True, default=None)
@click.option('--CHERRY-params', required=False, default=None)
@click.option('--PhaTYP-params', required=False, default=None)
@click.option('--iphop-cutoff', required=False, default=None)
@click.option('--iphop-params', required=False, default=None)
def run_viral_host(fasta, outdir, CHERRY_params, PhaTYP_params, iphop_cutoff, iphop_params):
        logging.info(f"Running module: viral-host")
        logging.info(f"fasta: {fasta}, outdir: {outdir}")
        
        module_obj = ViralHostModule()
        module_obj.name = "viral-host"
        # Set the attributes of the module object
        module_obj.fasta = fasta 
        module_obj.outdir = outdir 

        if CHERRY_params:
            module_obj.CERRY_params = CHERRY_params
            module_obj.hasOptions = True
        if PhaTYP_params:
            module_obj.PhaTYP_params = PhaTYP_params
            module_obj.hasOptions = True
        if iphop_cutoff:
            module_obj.iphop_cutoff = iphop_cutoff
            module_obj.hasOptions = True
        if iphop_params:
            module_obj.iphop_params = iphop_params
            module_obj.hasOptions = True

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("viral-host", module_obj)
        logging.info(f"End module run")

@cli.command(
    'viral-community',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Viral Community module'
)
@click.option('--outdir', required=True, default=None)
@click.option('--datadir', required=False, default=None)
@click.option('--samplelist', required=False, default=None)
@click.option('--mpa-indexv', required=False, default=None)
@click.option('--mpa-params', required=False, default=None)
def run_viral_community(outdir, datadir, samplelist, mpa_indexv, mpa_params):
        logging.info(f"Running module: viral-community")
        logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = ViralCommunityModule()
        module_obj.name = "viral-community"
        # Set the attributes of the module object
        module_obj.outdir = outdir 
        module_obj.datadir = datadir
        module_obj.samplelist = samplelist

        if mpa_indexv:
            module_obj.mpa_indexv = mpa_indexv
            module_obj.hasOptions = True
        if mpa_params:
            module_obj.mpa_params = mpa_params
            module_obj.hasOptions = True

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("viral-community", module_obj)
        logging.info(f"End module run")

@cli.command(
    'viral-annotate',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Viral Annotate module'
)
@click.option('--outdir', required=True, default=None)
@click.option('--datadir', required=False, default=None)
@click.option('--samplelist', required=False, default=None)
@click.option('--eggNOG-params', required=False, default=None)
@click.option('--PhaVIP-params', required=False, default=None)
def run_viral_annotate(outdir, datadir, samplelist, eggNOG_params, PhaVIP_params):
        logging.info(f"Running module: viral-annotate")
        logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = ViralAnnotateModule()
        module_obj.name = "viral-annotate"
        # Set the attributes of the module object
        module_obj.outdir = outdir 
        module_obj.datadir = datadir
        module_obj.samplelist = samplelist

        if eggNOG_params:
            module_obj.eggNOG_params = eggNOG_params
            module_obj.hasOptions = True
        if PhaVIP_params:
            module_obj.PhaVIP_params = PhaVIP_params
            module_obj.hasOptions = True

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("viral-annotate", module_obj)
        logging.info(f"End module run")

@cli.command(
    'prok-community',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Prokaryotic Community module'
)
@click.option('--outdir', required=True, default=None)
@click.option('--datadir', required=False, default=None)
@click.option('--samplelist', required=False, default=None)
@click.option('--mpa-params', required=False, default=None)
@click.option('--mpa-indexv', required=False, default=None)
def run_prok_community(outdir, datadir, samplelist, mpa_params, mpa_indexv):
        logging.info(f"Running module: prok-community")
        logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = ProkaryoticCommunityModule()
        module_obj.name = "prok-community"
        # Set the attributes of the module object
        module_obj.outdir = outdir 
        module_obj.datadir = datadir
        module_obj.samplelist = samplelist

        if mpa_params:
            module_obj.mpa_params = mpa_params
            module_obj.hasOptions = True    
        if mpa_indexv:
            module_obj.mpa_indexv = mpa_indexv
            module_obj.hasOptions = True

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("prok-community", module_obj)
        logging.info(f"End module run")

# TBD ProkaryoticBinningModule


@cli.command(
    'prok-annotate',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Prokaryotic Annotate module'
)
@click.option('--outdir', required=True, default=None)
@click.option('--datadir', required=False, default=None)
@click.option('--samplelist', required=False, default=None)
def run_prok_annotate(outdir, datadir, samplelist):
        logging.info(f"Running module: prok-annotate")
        logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = ProkaryoticAnnotateModule()
        module_obj.name = "prok-annotate"
        # Set the attributes of the module object
        module_obj.outdir = outdir 
        module_obj.datadir = datadir
        module_obj.samplelist = samplelist

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("prok-annotate", module_obj)
        logging.info(f"End module run")

@cli.command(
    'end-to-end',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the End-To-End module'
)
@click.option('--outdir', required=True, default=None)
@click.option('--datadir', required=False, default=None)
@click.option('--samplelist', required=False, default=None)
def run_end_to_end(outdir, datadir, samplelist):
        logging.info(f"Running module: end-to-end")
        logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = EndToEndModule()
        module_obj.name = "end-to-end"
        # Set the attributes of the module object
        module_obj.outdir = outdir 
        module_obj.datadir = datadir
        module_obj.samplelist = samplelist

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("end-to-end", module_obj)
        logging.info(f"End module run")

@cli.command(
    'cluster-fast',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Cluster Fast module'
)
@click.option('--fasta', required=True, default=None)
@click.option('--outdir', required=True, default=None)
@click.option('--clustering-fast', required=False, default=None)
@click.option('--cdhit-params', required=False, default=None)
@click.option('--vOTU-ani', required=False, default=None)
@click.option('--vOTU-targetcov', required=False, default=None)
@click.option('--vOTU-querycov', required=False, default=None)
def run_cluster_fast(fasta, outdir, clustering_fast, cdhit_params, vOTU_ani, vOTU_targetcov, vOTU_querycov):
        logging.info(f"Running module: cluster-fast")
        logging.info(f"fasta: {fasta}, outdir: {outdir}")
        
        module_obj = ClusterFastModule()
        module_obj.name = "cluster-fast"
        # Set the attributes of the module object
        module_obj.fasta = fasta 
        module_obj.outdir = outdir 

        if clustering_fast:
            module_obj.clustering_fast = clustering_fast
            module_obj.hasOptions = True
        if cdhit_params:
            module_obj.cdhit_params = cdhit_params
            module_obj.hasOptions = True
        if vOTU_ani:
            module_obj.vOTU_ani = vOTU_ani
            module_obj.hasOptions = True
        if vOTU_targetcov:
            module_obj.vOTU_targetcov = vOTU_targetcov
            module_obj.hasOptions = True
        if vOTU_querycov:
            module_obj.vOTU_querycov = vOTU_querycov
            module_obj.hasOptions = True

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("cluster-fast", module_obj)
        logging.info(f"End module run")

@cli.command(
    'checkv-pyhmmer',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the CheckV PyHMMER module'
)
@click.option('--fasta', required=True, default=None)
@click.option('--outdir', required=True, default=None)
@click.option('--checkv-original', required=False, default=None)
@click.option('--checkv-params', required=False, default=None)
@click.option('--checkv-database', required=False, default=None)
def run_checkv_pyhmmer(fasta, outdir, checkv_original, checkv_params, checkv_database):
        logging.info(f"Running module: checkv-pyhmmer")
        logging.info(f"fasta: {fasta}, outdir: {outdir}")
        
        module_obj = CheckVPyHMMERModule()
        module_obj.name = "checkv-pyhmmer"
        # Set the attributes of the module object
        module_obj.fasta = fasta 
        module_obj.outdir = outdir 

        if checkv_original:
            module_obj.checkv_original = checkv_original
            module_obj.hasOptions = True
        if checkv_params:
            module_obj.checkv_params = checkv_params
            module_obj.hasOptions = True
        if checkv_database:
            module_obj.checkv_database = checkv_database
            module_obj.hasOptions = True

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("checkv-pyhmmer", module_obj)
        logging.info(f"End module run")

@cli.command(
    'setup-database',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Setup Database module'
)
@click.option('--fasta', required=True, default=None)
@click.option('--outdir', required=True, default=None)
@click.option('--PhaBox2-db', required=False, default=None)
@click.option('--genomad-db', required=False, default=None)
@click.option('--checkv-db', required=False, default=None)
@click.option('--eggNOG-db', required=False, default=None)
@click.option('--eggNOG-db-params', required=False, default=None)
@click.option('--virsorter2-db', required=False, default=None)
@click.option('--iphop-db', required=False, default=None)
@click.option('--humann-db', required=False, default=None)
def run_setup_database(fasta, outdir, PhaBox2_db, genomad_db, checkv_db, eggNOG_db, eggNOG_db_params, virsorter2_db, iphop_db, humann_db):
        logging.info(f"Running module: setup-database")
        logging.info(f"fasta: {fasta}, outdir: {outdir}")
        
        module_obj = SetupDatabaseModule()
        module_obj.name = "setup-pyhdatabasemmer"
        # Set the attributes of the module object
        module_obj.fasta = fasta 
        module_obj.outdir = outdir 

        if PhaBox2_db:  
            module_obj.PhaBox2_db = PhaBox2_db
            module_obj.hasOptions = True
        if genomad_db:
            module_obj.genomad_db = genomad_db
            module_obj.hasOptions = True
        if checkv_db:
            module_obj.checkv_db = checkv_db
            module_obj.hasOptions = True
        if eggNOG_db:
            module_obj.eggNOG_db = eggNOG_db
            module_obj.hasOptions = True
        if eggNOG_db_params:
            module_obj.eggNOG_db_params = eggNOG_db_params
            module_obj.hasOptions = True
        if virsorter2_db:   
            module_obj.virsorter2_db = virsorter2_db
            module_obj.hasOptions = True
        if iphop_db:
            module_obj.iphop_db = iphop_db
            module_obj.hasOptions = True
        if humann_db:
            module_obj.humann_db = humann_db
            module_obj.hasOptions = True

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("setup-database", module_obj)
        logging.info(f"End module run")

