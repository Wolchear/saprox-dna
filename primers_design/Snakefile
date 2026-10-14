import pandas as pd
from snakemake.utils import min_version
min_version("6.0")

from workflow.lib.utils import get_path

configfile: "config/config.yaml"

samples = pd.read_csv("config/samples.csv")
TREES = set(samples["tree"].tolist())


SUBUNITS_RAW = expand(
    "{raw_fata}/{subunit}.fasta.gz",
    raw_fata = get_path(config['data'], 'raw'),
    subunit = ["SSU", "LSU"]
)

ORGANISMS_RAW = expand(
    "{raw_fata}/{subunit}_{organism}.fasta.gz",
    raw_fata = get_path(config['data'], 'raw'),
    subunit = ["SSU", "LSU"],
    organism = ['fungi', 'trees']
)

rule all:
    input:
        SUBUNITS_RAW,
        ORGANISMS_RAW
        
RULES_DIR = get_path(config['workflow'], "rules")


module raw_data:
    snakefile: f"{RULES_DIR}/raw_data.smk"
    config: config
use rule * from raw_data