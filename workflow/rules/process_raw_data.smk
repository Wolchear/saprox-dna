import pandas as pd
from workflow.lib.utils import get_path

NANOPORE_DATA_DIR = config['raw_nanopore_data']
RAW_BARCODES_DIR = get_path(config['data'], 'barcodes')
TRIMMED_ADAPTERS_DIR = get_path(config['data'], 'trimmed_adapters')
TRIMMED_BARCODES_DIR = get_path(config['data'], 'trimmed_barcodes')
TRIMMED_PRIMERS_DIR = get_path(config['data'], 'trimmed_primers')
FILTERED_BARCODES_DIR= get_path(config['data'], 'filtered_barcodes')

FA2TSV_DIR = get_path(config['qc'], 'fa2tsv')
LOCATE_DIR = get_path(config['qc'], 'locate')
JOINS_DIR = get_path(config['qc'], 'joins')
JOINS_PLOTS_DIR = get_path(config['qc'], 'contamination_plots')
STATUS = "before|after"
FASTA_DIRS = "trimmed_adapters|trimmed_barcodes|trimmed_primers"
SEQ_STATS = get_path(config['qc'], 'seq_stats')
REPORTS_DIR = get_path(config['qc'], 'trimming_reports')

NANOPLOTS_DIR = get_path(config['qc'], 'nanoplots')

samples = pd.read_csv("config/samples.tsv", sep="\t")
BARCODE_IDS = samples["barcode"].tolist()

BARCODE_FW = dict(zip(samples["barcode"], samples["Forward"]))
BARCODE_REV = dict(zip(samples["barcode"], samples["Reverse"]))

SCRIPTS = get_path(config['workflow'], 'scripts')

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

INPUT_FASTA = {
    "before":{
        "trimmed_adapters": "merged_barcodes",
        "trimmed_barcodes": "trimmed_adapters",
        "trimmed_primers": "trimmed_barcodes"
    },
    "after": {
        "trimmed_adapters": "trimmed_adapters",
        "trimmed_barcodes": "trimmed_barcodes",
        "trimmed_primers": "trimmed_primers"
    }
}

rule get_seqkit_fx2tab:
    input:
        lambda wc: f"data/{INPUT_FASTA[wc.status][wc.fasta_dir]}/{{barcode}}.fastq.gz"
    output:
        f"{FA2TSV_DIR}/{{fasta_dir}}/{{status}}/{{barcode}}.tsv"
    threads: max(1, config['max_threads'] // 2)
    wildcard_constraints:
        fasta_dir=FASTA_DIRS,
        status = STATUS
    log:
        f"logs/process_raw_data/get_seqkit_fx2tab/{{fasta_dir}}/{{status}}/{{barcode}}.log"
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
        "trimmed_adapters": "adapters.fasta",
        "trimmed_barcodes": "barcodes.fasta",
        "trimmed_primers": "primers.fasta"
}

MISMATCH = {
        "trimmed_adapters": "3",
        "trimmed_barcodes": "1",
        "trimmed_primers": "1"
}


rule get_seqkit_locate:
    input:
        lambda wc: f"data/{INPUT_FASTA[wc.status][wc.fasta_dir]}/{{barcode}}.fastq.gz"
    output:
        f"{LOCATE_DIR}/{{fasta_dir}}/{{status}}/{{barcode}}.tsv"
    threads: max(1, config['max_threads'] // 2)
    wildcard_constraints:
        fasta_dir=FASTA_DIRS,
        status = STATUS
    log:
        f"logs/process_raw_data/get_seqkit_locate/{{fasta_dir}}/{{status}}/{{barcode}}.log"
    conda:
        '../envs/process_raw_data.yml'
    params:
        adapters_file = lambda wc: f"config/{ADAPTERS[wc.fasta_dir]}",
        m = lambda wc: MISMATCH[wc.fasta_dir]
    shell:
        r"""
        {{
            seqkit locate -m {params.m} \
                -f {params.adapters_file} \
                {input} \
            | awk 'NR>1 {{print $1,$2,$3,$4,$5,$6,$7}}' OFS='\t' \
            | sort -k1,1 
        }} > {output} 2> {log}
        """

rule merge_hits:
    input:
        hits = f"{LOCATE_DIR}/{{fasta_dir}}/{{status}}/{{barcode}}.tsv",
        length = f"{FA2TSV_DIR}/{{fasta_dir}}/{{status}}/{{barcode}}.tsv"
    output:
        f"{JOINS_DIR}/{{fasta_dir}}/{{status}}/{{barcode}}.tsv"
    threads: 1
    wildcard_constraints:
        fasta_dir=FASTA_DIRS,
        status = STATUS
    log:
        f"logs/process_raw_data/merge_raw_hits/{{fasta_dir}}/{{status}}/{{barcode}}.log"
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

rule plot_joins:
    input:
        f"{JOINS_DIR}/{{fasta_dir}}/{{status}}/{{barcode}}.tsv"
    output:
        f"{JOINS_PLOTS_DIR}/{{fasta_dir}}/{{status}}/{{barcode}}.png"
    threads: 1
    wildcard_constraints:
        fasta_dir=FASTA_DIRS,
        status = STATUS
    log:
        f"logs/process_raw_data/plot_joins/{{fasta_dir}}/{{status}}/{{barcode}}.log"
    conda:
        '../envs/process_raw_data.yml'
    params:
        script = f"{SCRIPTS}/plot_body_coverage.py"
    shell:
        """
        python3 {params.script} --input {input} --output {output} > {log} 2>&1
        """
            
            
rule trim_ends:
    input:
        f"{RAW_BARCODES_DIR}/{{barcode}}.fastq.gz"
    output:
        fasta = f"{TRIMMED_ADAPTERS_DIR}/{{barcode}}.fastq.gz",
        json = f"{REPORTS_DIR}/trim_ends/{{barcode}}.json",
        html = f"{REPORTS_DIR}/trim_ends/{{barcode}}.html",
    threads: max(1, config['max_threads'] // 2)
    conda:
        '../envs/process_raw_data.yml'
    log:
        f"logs/process_raw_data/trim_ends/{{barcode}}.log"
    params:
        adapters_file =  f"config/adapters.fasta"
    shell:
        r"""
        fastplong -i {input} -o {output.fasta} \
            --thread {threads} \
            -d 0.20 \
            -a {params.adapters_file} \
            --length_required 500 \
            --length_limit 2500 \
            --verbose \
            --json {output.json} \
            --html {output.html} \
            > {log} 2>&1
        """

rule trim_barcodes:
    input:
        f"{TRIMMED_ADAPTERS_DIR}/{{barcode}}.fastq.gz",
    output:
        fasta = f"{TRIMMED_BARCODES_DIR}/{{barcode}}.fastq.gz",
        json = f"{REPORTS_DIR}/barcodes/{{barcode}}.json",
        html = f"{REPORTS_DIR}/barcodes/{{barcode}}.html",
    threads: max(1, config['max_threads'] // 2)
    conda:
        '../envs/process_raw_data.yml'
    log:
        f"logs/process_raw_data/trim_barcodes/{{barcode}}.log"
    params:
        fw = lambda wc: BARCODE_FW[wc.barcode],
        rev = lambda wc: BARCODE_REV[wc.barcode]
    shell:
        r"""
        fastplong -i {input} -o {output.fasta} \
            --thread {threads} \
            -d 0.20 \
            -s {params.fw} \
            -e {params.rev} \
            --length_required 500 \
            --length_limit 2500 \
            --verbose \
            --json {output.json} \
            --html {output.html} \
            > {log} 2>&1
        """

rule trim_primers:
    input:
        f"{TRIMMED_BARCODES_DIR}/{{barcode}}.fastq.gz"
    output:
        fasta = f"{TRIMMED_PRIMERS_DIR}/{{barcode}}.fastq.gz",
        json = f"{REPORTS_DIR}/primers/{{barcode}}.json",
        html = f"{REPORTS_DIR}/primers/{{barcode}}.html",
    threads: max(1, config['max_threads'] // 2)
    conda:
        '../envs/process_raw_data.yml'
    log:
        f"logs/process_raw_data/trim_primers/{{barcode}}.log"
    params:
        adapters_file =  f"config/primers_to_trim.fasta"
    shell:
        r"""
        fastplong -i {input} -o {output.fasta} \
            --thread {threads} \
            -d 0.20 \
            -a {params.adapters_file} \
            --length_required 500 \
            --trimming_extension 0 \
            --length_limit 2500 \
            --verbose \
            --json {output.json} \
            --html {output.html} \
            > {log} 2>&1
        """


TABLE_INPUT = {
    'raw_stats': RAW_BARCODES_DIR,
    'trimmed_adapters_stats': TRIMMED_ADAPTERS_DIR,
    'trimmed_barcodes_stats': TRIMMED_BARCODES_DIR,
    'trimmed_primers_stats': TRIMMED_PRIMERS_DIR,
    'filtered_barcodes_stats': FILTERED_BARCODES_DIR
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

rule plot_nanoplots:
    input:
        fasta = f"{TRIMMED_PRIMERS_DIR}/{{barcode}}.fastq.gz"
    output:
        f"{NANOPLOTS_DIR}/{{barcode}}/{{barcode}}_NanoPlot-report.html"
    threads: max(1, config['max_threads'])
    conda:
        '../envs/nanoplot.yml'
    log:
        f"logs/process_raw_data/plot_nanoplots/{{barcode}}.log"
    params:
        prefix = lambda wc: f"{wc.barcode}_",
        title = lambda wc: f"{wc.barcode}_trimmed_primers",
        outdir = lambda wc: f"{NANOPLOTS_DIR}/{wc.barcode}"
    shell:
        r"""
        NanoPlot \
            --fastq {input.fasta} \
            --outdir {params.outdir} \
            --threads {threads} \
            --plots dot \
            --N50 \
            --title {params.title} \
            --prefix {params.prefix} \
            > {log} 2>&1
        """

rule final_filtering:
    input:
        f"{TRIMMED_PRIMERS_DIR}/{{barcode}}.fastq.gz"
    output:
        fasta = f"{FILTERED_BARCODES_DIR}/{{barcode}}.fastq.gz",
        json = f"{REPORTS_DIR}/final_filtering/{{barcode}}.json",
        html = f"{REPORTS_DIR}/final_filtering/{{barcode}}.html",
    threads: max(1, config['max_threads'] // 2)
    conda:
        '../envs/process_raw_data.yml'
    log:
        f"logs/process_raw_data/final_filtering/{{barcode}}.log"
    shell:
        r"""
        fastplong -i {input} -o {output.fasta} \
            --disable_adapter_trimming \
            --thread {threads} \
            --mean_qual 11 \
            --length_required 500 \
            --length_limit 1500 \
            --verbose \
            --json {output.json} \
            --html {output.html} \
            > {log} 2>&1
        """