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
    table = ['raw_stats', 'trimmed_adapters_stats', 'trimmed_barcodes_stats']
)

JOINS_PLOTS = expand(
    "{tables_dir}/{sub_dir}/{status}/{table}.png",
    tables_dir = get_path(config['qc'], 'contamination_plots'),
    status = ['before', 'after'],
    sub_dir = [ 'trimmed_adapters', 'trimmed_barcodes'],
    table = BARCODE_IDS
)

rule all:
    input:
        STAT_TABLES,
        JOINS_PLOTS,


RULES_DIR = get_path(config['workflow'], "rules")

module process_raw_data:
    snakefile: f"{RULES_DIR}/process_raw_data.smk"
    config: config
use rule * from process_raw_data