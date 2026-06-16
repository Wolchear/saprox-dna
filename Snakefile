import pandas as pd
from snakemake.utils import min_version
min_version("6.0")

from workflow.lib.utils import get_path

configfile: "config/config.yaml"

samples = pd.read_csv("config/samples.tsv", sep="\t")
BARCODE_IDS = samples["barcode"].tolist()


STAT_TABLES = expand(
    "{tables_dir}/{table}.tsv",
    tables_dir = get_path(config['qc'], 'seq_stats'),
    table = [
        'raw_stats',
        'trimmed_adapters_stats',
        'trimmed_barcodes_stats',
        'trimmed_primers_stats',
        "filtered_barcodes_stats"
    ]
)

JOINS_PLOTS = expand(
    "{tables_dir}/{sub_dir}/{status}/{table}.png",
    tables_dir = get_path(config['qc'], 'contamination_plots'),
    status = ['before', 'after'],
    sub_dir = [ 'trimmed_adapters', 'trimmed_barcodes', 'trimmed_primers'],
    table = BARCODE_IDS
)

NANOPLOTS = expand(
    "{nanoplots_dir}/{barcode}/{barcode}_NanoPlot-report.html",
    nanoplots_dir = get_path(config['qc'], 'nanoplots'),
    barcode = BARCODE_IDS
)

FILTERED_BARCODES = expand(
    "{filtered_dir}/{barcode}.fastq.gz",
    filtered_dir = get_path(config['data'], 'filtered_barcodes'),
    barcode = BARCODE_IDS
)

EXCLUDE_BARCODES = ["barcode07", "barcode18"]

ANALYSIS_BARCODES = [
    b for b in BARCODE_IDS
    if b not in EXCLUDE_BARCODES
]

CONSESUS_FILES = expand(
    "{consensus_dir}/{barcode}/consensusfile.fasta",
    consensus_dir = get_path(config['output'], 'amplicon_sort'),
    barcode = ANALYSIS_BARCODES
)

KRAKEN_REPORTS = expand(
    "{kraken_dir}/{barcode}/report.tsv",
    kraken_dir = get_path(config['output'], 'amplicons_classification'),
    barcode = ANALYSIS_BARCODES
)

rule all:
    input:
        STAT_TABLES,
        JOINS_PLOTS,
        NANOPLOTS,
        FILTERED_BARCODES,
        KRAKEN_REPORTS

RULES_DIR = get_path(config['workflow'], "rules")

module process_raw_data:
    snakefile: f"{RULES_DIR}/process_raw_data.smk"
    config: config
use rule * from process_raw_data


module amplicon_classification:
    snakefile: f"{RULES_DIR}/amplicon_classification.smk"
    config: config
use rule * from amplicon_classification