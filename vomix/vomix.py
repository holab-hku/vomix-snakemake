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
@click.option('--decontam-host', prompt='Decontamination host', default=True, 
required=True)
@click.option('--outdir', prompt='Output directory', required=True, default=None)
@click.option('--datadir', prompt='Data directory', required=True, default=None)
@click.option('--samplelist', prompt='Sample list file', required=True, default=None)
def run_preprocess(decontam_host, outdir, datadir, samplelist):
        logging.info(f"Running module: preprocess")
        logging.info(f"decontamHost: {decontam_host}, outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = PreProcessingModule()
        module_obj.name = "preprocess"
        # Set the attributes of the module object
        module_obj.decontamHost = decontam_host
        module_obj.outdir = outdir 
        module_obj.datadir = datadir
        module_obj.samplelist = samplelist

        vomix_actions_instance = vomix_actions()
        out = vomix_actions_instance.run_module("preprocess", module_obj)
        logging.info(f"End module run: {out}")

@cli.command(
    'assembly',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Assembly & Co-assembly module'
)
@click.option('--assembler', prompt='Assembler', default="megahit", 
required=True)
@click.option('--outdir', prompt='Output directory', required=True, default=None)
@click.option('--datadir', prompt='Data directory', required=True, default=None)
@click.option('--samplelist', prompt='Sample list file', required=True, default=None)
def run_assembly(assembler, outdir, datadir, samplelist):
        logging.info(f"Running module: assembly")
        logging.info(f"assembler: {assembler}, outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = AssemblyCoAssemblyModule()
        module_obj.name = "assembly"
        # Set the attributes of the module object
        module_obj.assembler = assembler
        module_obj.outdir = outdir 
        module_obj.datadir = datadir
        module_obj.samplelist = samplelist

        vomix_actions_instance = vomix_actions()
        out = vomix_actions_instance.run_module("assembler", module_obj)
        logging.info(f"End module run: {out}")

@cli.command(
    'viral-identify',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Viral Identify module'
)
@click.option('--outdir', prompt='Output directory', required=True, default=None)
@click.option('--datadir', prompt='Data directory', required=True, default=None)
@click.option('--samplelist', prompt='Sample list file', required=True, default=None)
def run_viral_identify(outdir, datadir, samplelist):
        logging.info(f"Running module: viral-identify")
        logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = ViralIdentifyModule()
        module_obj.name = "viral-identify"
        # Set the attributes of the module object
        module_obj.outdir = outdir 
        module_obj.datadir = datadir
        module_obj.samplelist = samplelist

        vomix_actions_instance = vomix_actions()
        out = vomix_actions_instance.run_module("viral-identify", module_obj)
        logging.info(f"End module run: {out}")


@cli.command(
    'viral-taxonomy',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Viral Taxonomy module'
)
@click.option('--fasta', prompt='Fasta', required=True, default=None)
@click.option('--outdir', prompt='Output directory', required=True, default=None)
def run_viral_taxonomy(fasta, outdir):
        logging.info(f"Running module: viral-taxonomy")
        logging.info(f"fasta: {fasta}, outdir: {outdir}")
        
        module_obj = ViralTaxonomyModule()
        module_obj.name = "viral-taxonomy"
        # Set the attributes of the module object
        module_obj.fasta = fasta 
        module_obj.outdir = outdir 

        vomix_actions_instance = vomix_actions()
        out = vomix_actions_instance.run_module("viral-taxonomy", module_obj)
        logging.info(f"End module run: {out}")

@cli.command(
    'viral-host',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Viral Host module'
)
@click.option('--fasta', prompt='Fasta', required=True, default=None)
@click.option('--outdir', prompt='Output directory', required=True, default=None)
def run_viral_host(fasta, outdir):
        logging.info(f"Running module: viral-host")
        logging.info(f"fasta: {fasta}, outdir: {outdir}")
        
        module_obj = ViralHostModule()
        module_obj.name = "viral-host"
        # Set the attributes of the module object
        module_obj.fasta = fasta 
        module_obj.outdir = outdir 

        vomix_actions_instance = vomix_actions()
        out = vomix_actions_instance.run_module("viral-host", module_obj)
        logging.info(f"End module run: {out}")

@cli.command(
    'viral-community',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Viral Community module'
)
@click.option('--outdir', prompt='Output directory', required=True, default=None)
@click.option('--datadir', prompt='Data directory', required=True, default=None)
@click.option('--samplelist', prompt='Sample list file', required=True, default=None)
def run_viral_community(outdir, datadir, samplelist):
        logging.info(f"Running module: viral-community")
        logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = ViralCommunityModule()
        module_obj.name = "viral-community"
        # Set the attributes of the module object
        module_obj.outdir = outdir 
        module_obj.datadir = datadir
        module_obj.samplelist = samplelist

        vomix_actions_instance = vomix_actions()
        out = vomix_actions_instance.run_module("viral-community", module_obj)
        logging.info(f"End module run: {out}")

@cli.command(
    'viral-annotate',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Viral Annotate module'
)
@click.option('--outdir', prompt='Output directory', required=True, default=None)
@click.option('--datadir', prompt='Data directory', required=True, default=None)
@click.option('--samplelist', prompt='Sample list file', required=True, default=None)
def run_viral_annotate(outdir, datadir, samplelist):
        logging.info(f"Running module: viral-annotate")
        logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = ViralAnnotateModule()
        module_obj.name = "viral-annotate"
        # Set the attributes of the module object
        module_obj.outdir = outdir 
        module_obj.datadir = datadir
        module_obj.samplelist = samplelist

        vomix_actions_instance = vomix_actions()
        out = vomix_actions_instance.run_module("viral-annotate", module_obj)
        logging.info(f"End module run: {out}")

@cli.command(
    'prok-community',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Prokaryotic Community module'
)
@click.option('--outdir', prompt='Output directory', required=True, default=None)
@click.option('--datadir', prompt='Data directory', required=True, default=None)
@click.option('--samplelist', prompt='Sample list file', required=True, default=None)
def run_prok_community(outdir, datadir, samplelist):
        logging.info(f"Running module: prok-community")
        logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = ProkaryoticCommunityModule()
        module_obj.name = "prok-community"
        # Set the attributes of the module object
        module_obj.outdir = outdir 
        module_obj.datadir = datadir
        module_obj.samplelist = samplelist

        vomix_actions_instance = vomix_actions()
        out = vomix_actions_instance.run_module("prok-community", module_obj)
        logging.info(f"End module run: {out}")

# TBD ProkaryoticBinningModule


@cli.command(
    'prok-annotate',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Prokaryotic Annotate module'
)
@click.option('--outdir', prompt='Output directory', required=True, default=None)
@click.option('--datadir', prompt='Data directory', required=True, default=None)
@click.option('--samplelist', prompt='Sample list file', required=True, default=None)
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
        out = vomix_actions_instance.run_module("prok-annotate", module_obj)
        logging.info(f"End module run: {out}")

@cli.command(
    'end-to-end',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the End-To-End module'
)
@click.option('--outdir', prompt='Output directory', required=True, default=None)
@click.option('--datadir', prompt='Data directory', required=True, default=None)
@click.option('--samplelist', prompt='Sample list file', required=True, default=None)
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
        out = vomix_actions_instance.run_module("end-to-end", module_obj)
        logging.info(f"End module run: {out}")

@cli.command(
    'cluster-fast',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Cluster Fast module'
)
@click.option('--fasta', prompt='Fasta', required=True, default=None)
@click.option('--outdir', prompt='Output directory', required=True, default=None)
def run_cluster_fast(fasta, outdir):
        logging.info(f"Running module: cluster-fast")
        logging.info(f"fasta: {fasta}, outdir: {outdir}")
        
        module_obj = ClusterFastModule()
        module_obj.name = "cluster-fast"
        # Set the attributes of the module object
        module_obj.fasta = fasta 
        module_obj.outdir = outdir 

        vomix_actions_instance = vomix_actions()
        out = vomix_actions_instance.run_module("cluster-fast", module_obj)
        logging.info(f"End module run: {out}")

@cli.command(
    'checkv-pyhmmer',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the CheckV PyHMMER module'
)
@click.option('--fasta', prompt='Fasta', required=True, default=None)
@click.option('--outdir', prompt='Output directory', required=True, default=None)
def run_checkv_pyhmmer(fasta, outdir):
        logging.info(f"Running module: checkv-pyhmmer")
        logging.info(f"fasta: {fasta}, outdir: {outdir}")
        
        module_obj = CheckVPyHMMERModule()
        module_obj.name = "checkv-pyhmmer"
        # Set the attributes of the module object
        module_obj.fasta = fasta 
        module_obj.outdir = outdir 

        vomix_actions_instance = vomix_actions()
        out = vomix_actions_instance.run_module("checkv-pyhmmer", module_obj)
        logging.info(f"End module run: {out}")

@cli.command(
    'setup-database',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Setup Database module'
)
@click.option('--fasta', prompt='Fasta', required=True, default=None)
@click.option('--outdir', prompt='Output directory', required=True, default=None)
def run_setup_database(fasta, outdir):
        logging.info(f"Running module: setup-database")
        logging.info(f"fasta: {fasta}, outdir: {outdir}")
        
        module_obj = SetupDatabaseModule()
        module_obj.name = "setup-pyhdatabasemmer"
        # Set the attributes of the module object
        module_obj.fasta = fasta 
        module_obj.outdir = outdir 

        vomix_actions_instance = vomix_actions()
        out = vomix_actions_instance.run_module("setup-database", module_obj)
        logging.info(f"End module run: {out}")


# @cli.command(
#     'run',
#     context_settings=dict(ignore_unknown_options=True),
#     short_help='Run a vomix module'
# )

# @click.argument('module')
# @click.option('--useLastOptions', prompt='Re-use last options from the previous run?', is_flag=True, default=False, callback=useLastOptionsCheck)
# @click.option('--fasta', prompt='Fasta', required=False, default=None)
# @click.option('--outdir', prompt='Output directory', required=False, default=None)
# @click.option('--datadir', prompt='Data directory', required=False, default=None)
# @click.option('--assembler', prompt='Assembler', default='megahit', required=False)
# @click.option('--samplelist', prompt='Sample list file', required=False, default=None)
# @click.option('--decontam-host', prompt='Decontamination host', default=True, required=False)
# @click.option('--binning-consensus', prompt='Binning consensus', default=False , required=False)
# @click.option('--email', prompt='Email')
# def run_module(module, uselastoptions, splits, fasta, outdir, datadir, assembler, samplelist, decontam_host, binning_consensus, email):

#     # Get the module object from the modules dictionary
#     module_obj = modules.get(module)

#     # Check if the module exists and its name matches the provided module parameter
#     if not module_obj or getattr(module_obj, "name", None) != module:
#         logging.error(f"Unknown or mismatched module: {module}")
#         sys.exit(1)

#     module_options = {
#         "name": module,
#         "splits": int(splits),
#         "fasta": fasta,
#         "outdir": outdir,
#         "datadir": datadir,
#         "assembler": assembler,
#         "samplelist": samplelist,
#         "decontam_host": decontam_host,
#         "binning_consensus": binning_consensus,
#         "email": email,
#     }

#     for key, value in module_options.items():
#         if hasattr(module_obj, key):
#             setattr(module_obj, key, value)

#     if module in modules_list:
#         if outdir == None:
#             logging.info("Re-using last options from the previous run")
#             out = vomix_actions.run_last_module(module)
#             logging.info(f"End module run: {out}")
#         else:
#             logging.info(f"Running module: {module}")
#             logging.info(f"outdir: {outdir}, datadir: {datadir}, assembler: {assembler}, samplelist: {samplelist}, decontam_host: {decontam_host}, binning_consensus: {binning_consensus}")
            
#             out = vomix_actions.run_module(module, module_obj)
#         logging.info(f"End module run: {out}")

#     # test case
#     elif module == "test":
#         if outdir == None:
#             logging.info("Re-using last options from the previous run")
#             out = vomix_actions.run_last_module(module)
#         else:
#             logging.info(f"[TEST] Running module: {module}")
#             logging.info(f"outdir: {outdir}, datadir: {datadir}, assembler: {assembler}, samplelist: {samplelist}, decontam_host: {decontam_host}, binning_consensus: {binning_consensus}")

#             out = vomix_actions.run_module(module, module_obj)
#         logging.info(f"End module run. Output: {out}")

#     else:
#         logging.error(f"Unknown module: {module}")
#         sys.exit(1)

