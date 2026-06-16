import pandas as pd
from workflow.lib.utils import get_path

samples = pd.read_csv("config/samples.tsv", sep="\t")
BARCODE_IDS = samples["barcode"].tolist()

FILTERED_BARCODES_DIR= get_path(config['data'], 'filtered_barcodes')
CONSENSUS_DIR = get_path(config['output'], 'amplicon_sort')

EXTERNAL = get_path(config['workflow'], 'external')
SCRIPTS = get_path(config['workflow'], 'scripts')


rule sort_amplicons:
    input:
        f"{FILTERED_BARCODES_DIR}/{{barcode}}.fastq.gz"
    output:
        f"{CONSENSUS_DIR}/{{barcode}}/consensusfile.fasta",
    threads: max(1, config['max_threads'])
    conda:
        '../envs/amplicon_sorter.yml'
    log:
        f"logs/amplicon_classification/sort_amplicons/{{barcode}}.log"
    params:
        script = f"{EXTERNAL}/amplicon_sorter/amplicon_sorter.py",
        outdir = lambda wc: f"{CONSENSUS_DIR}/{wc.barcode}"
    shell:
        r"""
        python3 {params.script} \
            -i {input} \
            -o {params.outdir} \
            -np {threads} \
            --compressed \
            > {log} 2>&1
        """

