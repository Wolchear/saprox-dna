import pandas as pd
from snakemake.utils import min_version
min_version("6.0")

from workflow.lib.utils import get_path

configfile: "config/config.yaml"

BARCODE_IDS = (
    [f"barcode{i:02d}" for i in range(1, 15)] +
    [f"barcode{i:02d}" for i in range(17, 25)] +
    ["unclassified"]
)

BARCODES = expand(
    "{barcode_dir}/{barcode}.fastq.gz",
    barcode_dir = get_path(config['data'], 'barcodes'),
    barcode=BARCODE_IDS
)

rule all:
    input:
        BARCODES

RULES_DIR = get_path(config['workflow'], "rules")

module process_raw_data:
    snakefile: f"{RULES_DIR}/process_raw_data.smk"
    config: config
use rule * from process_raw_data