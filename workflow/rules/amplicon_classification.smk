import pandas as pd
from workflow.lib.utils import get_path

samples = pd.read_csv("config/samples.tsv", sep="\t")
BARCODE_IDS = samples["barcode"].tolist()

FILTERED_BARCODES_DIR= get_path(config['data'], 'filtered_barcodes')
CONSENSUS_DIR = get_path(config['output'], 'amplicon_sort')
AMPLICONS_CLASSIFICATION_DIR = get_path(config['output'], 'amplicons_classification')
ADJUSTED_REPORTS_DIR = get_path(config['output'], 'adjusted_kraken_reports')

EXTERNAL = get_path(config['workflow'], 'external')
SCRIPTS = get_path(config['workflow'], 'scripts')

KRAKEN_LINK = config['kraken_db']
KRAKEN_DB_DIR = get_path(config['data'], 'kraken_db')


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
            -maxr 1000000000000 \
            --allreads \
            -sg 99.0 \
            -ss 99.0 \
            -sc 99.0 \
            > {log} 2>&1
        """

rule merge_unique_consensus:
    input:
        f"{CONSENSUS_DIR}/{{barcode}}/consensusfile.fasta"
    output:
        f"{CONSENSUS_DIR}/{{barcode}}/merged.fasta",
    threads: 1
    conda:
        '../envs/amplicon_sorter.yml'
    log:
        f"logs/amplicon_classification/merge_uniqu_consensus/{{barcode}}.log"
    params:
        unique = lambda wc: f"{CONSENSUS_DIR}/{wc.barcode}/{wc.barcode}_nogroup_unique.fasta.gz"
    shell:
        r"""
        set -euo pipefail

        : > {log}

        cat {input} > {output}

        if [ -s "{params.unique}" ]; then
            echo "Adding unique file: {params.unique}" >> {log}
            zcat "{params.unique}" >> {output}
        else
            echo "No unique file for {wildcards.barcode}: {params.unique}" >> {log}
        fi
        """

rule get_db:
    output:
        temp(f"{KRAKEN_DB_DIR}/k2_pluspf_16_GB_20260226.tar.gz")
    threads: 1
    conda:
        '../envs/kraken.yml'
    log:
        f"logs/amplicon_classification/get_db.log"
    params:
        link = KRAKEN_LINK
    shell:
        """
        wget -O {output} {params.link} > {log} 2>&1
        """

rule unpack_db:
    input:
        f"{KRAKEN_DB_DIR}/k2_pluspf_16_GB_20260226.tar.gz"
    output:
        f"{KRAKEN_DB_DIR}/hash.k2d"
    threads: 1
    conda:
        '../envs/kraken.yml'
    log:
        f"logs/amplicon_classification/unpack_db.log"
    params:
        db = KRAKEN_DB_DIR
    shell:
        """
        tar -xzf {input} -C {params.db} > {log} 2>&1
        """

rule classify_amplicons:
    input:
        db = f"{KRAKEN_DB_DIR}/hash.k2d",
        amplicon_fasta = f"{CONSENSUS_DIR}/{{barcode}}/merged.fasta"
    output:
        classified = f"{AMPLICONS_CLASSIFICATION_DIR}/{{barcode}}/classified.fasta",
        unclassified = f"{AMPLICONS_CLASSIFICATION_DIR}/{{barcode}}/unclassified.fasta",
        report = f"{AMPLICONS_CLASSIFICATION_DIR}/{{barcode}}/report.tsv",
        kraken_output = f"{AMPLICONS_CLASSIFICATION_DIR}/{{barcode}}/kraken_output.tsv"
    threads: max(1, config['max_threads'])
    conda:
        '../envs/kraken.yml'
    log:
        f"logs/amplicon_classification/classify_amplicons/{{barcode}}.log"
    params:
        db = KRAKEN_DB_DIR
    shell:
        r"""
        kraken2 \
            --db {params.db} \
            --threads {threads} \
            --classified-out {output.classified} \
            --unclassified-out {output.unclassified} \
            --output {output.kraken_output} \
            --confidence 0 \
            --memory-mapping \
            --report {output.report} \
            {input.amplicon_fasta} > {log} 2>&1
        """

rule adjust_report:
    input:
        report = f"{AMPLICONS_CLASSIFICATION_DIR}/{{barcode}}/report.tsv",
        kraken_output = f"{AMPLICONS_CLASSIFICATION_DIR}/{{barcode}}/kraken_output.tsv"
    output:
        f"{ADJUSTED_REPORTS_DIR}/{{barcode}}_report.tsv",
    threads: 1
    conda:
        '../envs/kraken.yml'
    log:
        f"logs/amplicon_classification/adjust_report/{{barcode}}.log"
    params:
        script = f"{SCRIPTS}/adjust_weights_to_kraken_report.py"
    shell:
        r"""
        python3 {params.script} \
            --kraken-output {input.kraken_output} \
            --report {input.report} \
            --out {output} > {log} 2>&1
        """