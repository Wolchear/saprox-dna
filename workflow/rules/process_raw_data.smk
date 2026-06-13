import pandas as pd
from workflow.lib.utils import get_path

NANOPORE_DATA_DIR = config['raw_nanopore_data']
RAW_BARCODES_DIR = get_path(config['data'], 'barcodes')
TRIMMED_BARCODES_DIR = get_path(config['data'], 'trimmed')

SEQ_STATS = get_path(config['qc'], 'seq_stats')
REPORTS_DIR = get_path(config['qc'], 'trimming_reports')

samples = pd.read_csv("config/samples.tsv", sep="\t")
BARCODE_IDS = samples["barcode"].tolist()

rule merge_barcodes:
    output:
        f"{RAW_BARCODES_DIR}/{{barcode}}.fastq.gz"
    threads: 1
    log:
        f"logs/process_raw_data/merge_barcodes/{{barcode}}.log"
    params:
        failed_dir = lambda wc: f"{NANOPORE_DATA_DIR}/fastq_fail/{wc.barcode}",
        passed_dir = lambda wc: f"{NANOPORE_DATA_DIR}/fastq_pass/{wc.barcode}"
    shell:
        r"""
            : > {output}
            : > {log}

            if [ -d "{params.passed_dir}" ]; then
                find "{params.passed_dir}" -name "*.fastq.gz" -type f -exec cat {{}} + >> {output} 2>> {log}
            fi

            if [ -d "{params.failed_dir}" ]; then
                find "{params.failed_dir}" -name "*.fastq.gz" -type f -exec cat {{}} + >> {output} 2>> {log}
            fi
        """

rule trim_primers:
    input:
        f"{RAW_BARCODES_DIR}/{{barcode}}.fastq.gz"
    output:
        fasta = f"{TRIMMED_BARCODES_DIR}/{{barcode}}.fastq.gz",
        json = f"{REPORTS_DIR}/{{barcode}}.json",
        html = f"{REPORTS_DIR}/{{barcode}}.html",
    threads: max(1, config['max_threads'])
    conda:
        '../envs/process_raw_data.yml'
    log:
        f"logs/process_raw_data/trim_primers/{{barcode}}.log"
    shell:
        r"""
        fastplong -i {input} -o {output.fasta} \
            --thread {threads} \
            -d 0.20 \
            --trimming_extension 0 \
            --length_required 500 \
            --length_limit 2500 \
            --verbose \
            --json {output.json} \
            --html {output.html} \
            > {log} 2>&1
        """


TABLE_INPUT = {
    'raw_stats': RAW_BARCODES_DIR,
    'trimmed_stats': TRIMMED_BARCODES_DIR
}

rule get_seqkit_stats:
    input:
        expand(
            "{barcode_dir}/{barcode}.fastq.gz",
            barcode_dir = lambda wc: TABLE_INPUT[wc.table],
            barcode = BARCODE_IDS
        )
    output:
        f"{SEQ_STATS}/{{table}}.tsv"
    threads: max(1, config['max_threads'])
    conda:
        '../envs/process_raw_data.yml'
    log:
        f"logs/process_raw_data/get_seqkit_stats/{{table}}.log"
    shell:
        r"""
        seqkit stats \
            -j {threads} \
            --all \
            --tabular \
            --basename \
            -o {output} \
            {input} \
            > {log} 2>&1
        """