import pandas as pd
from workflow.lib.utils import get_path

NANOPORE_DATA_DIR = config['raw_nanopore_data']
RAW_BARCODES_DIR = get_path(config['data'], 'barcodes')
TRIMMED_BARCODES_DIR = get_path(config['data'], 'trimmed')

FA2TSV_DIR = get_path(config['qc'], 'fa2tsv')
LOCATE_DIR = get_path(config['qc'], 'locate')
JOINS_DIR = get_path(config['qc'], 'joins')

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

rule get_seqkit_fx2tab:
    input:
        f"data/{{fasta_dir}}/{{barcode}}.fastq.gz"
    output:
        f"{FA2TSV_DIR}/{{fasta_dir}}/{{barcode}}.tsv"
    threads: max(1, config['max_threads'] // 2)
    wildcard_constraints:
        fasta_dir="merged_barcodes|trimmed_barcodes"
    log:
        f"logs/process_raw_data/get_seqkit_fx2tab/{{fasta_dir}}/{{barcode}}.log"
    conda:
        '../envs/process_raw_data.yml'
    shell:
        r"""
        {{
            seqkit fx2tab -n -l -i {input} \
            | sort -k1,1 
        }} > {output} 2> {log}
        """

ADAPTERS = {
    "merged_barcodes": "adapters.fasta",
}

rule get_seqkit_locate:
    input:
        f"data/{{fasta_dir}}/{{barcode}}.fastq.gz"
    output:
        f"{LOCATE_DIR}/{{fasta_dir}}/{{barcode}}.tsv"
    threads: max(1, config['max_threads'] // 2)
    wildcard_constraints:
        fasta_dir="merged_barcodes"
    log:
        f"logs/process_raw_data/get_seqkit_locate/{{fasta_dir}}/{{barcode}}.log"
    conda:
        '../envs/process_raw_data.yml'
    params:
        adapters_file = lambda wc: f"config/{ADAPTERS[wc.fasta_dir]}"
    shell:
        r"""
        {{
            seqkit locate -m 3 \
                -f {params.adapters_file} \
                {input} \
            | awk 'NR>1 {{print $1,$2,$3,$4,$5,$6,$7}}' OFS='\t' \
            | sort -k1,1 
        }} > {output} 2> {log}
        """

rule merge_raw_hits:
    input:
        hits = f"{LOCATE_DIR}/merged_barcodes/{{barcode}}.tsv",
        length = f"{FA2TSV_DIR}/merged_barcodes/{{barcode}}.tsv"
    output:
        f"{JOINS_DIR}/merged_barcodes/{{barcode}}.tsv"
    threads: 1
    log:
        f"logs/process_raw_data/merge_raw_hits/merged_barcodes/{{barcode}}.log"
    conda:
        '../envs/process_raw_data.yml'
    shell:
        r"""
        {{
            join -t $'\t' -1 1 -2 1 \
                {input.length} \
                {input.hits} \
            | awk 'BEGIN{{
                    OFS="\t";
                    print "seqID","read_len","patternName","pattern","strand","start","end","matched","rel_start","rel_end"
                }}
                {{
                    read_len=$2;
                    start=$6;
                    end=$7;
                    rel_start=start/read_len*100;
                    rel_end=end/read_len*100;
                    print $1,read_len,$3,$4,$5,start,end,$8,rel_start,rel_end
                }}'
        }} > {output} 2> {log}
        """

rule trim_adapters:
    input:
        f"{RAW_BARCODES_DIR}/{{barcode}}.fastq.gz"
    output:
        fasta = f"{TRIMMED_BARCODES_DIR}/{{barcode}}.fastq.gz",
        json = f"{REPORTS_DIR}/{{barcode}}.json",
        html = f"{REPORTS_DIR}/{{barcode}}.html",
    threads: max(1, config['max_threads'] // 2)
    conda:
        '../envs/process_raw_data.yml'
    log:
        f"logs/process_raw_data/trim_adapters/{{barcode}}.log"
    params:
        fwd = config['adapters']['5'],
        rev = config['adapters']['3']
    shell:
        r"""
        fastplong -i {input} -o {output.fasta} \
            --thread {threads} \
            -d 0.20 \
            -s {params.fwd} \
            -e {params.rev} \
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