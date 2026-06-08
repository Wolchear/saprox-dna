from workflow.lib.utils import get_path

NANOPORE_DATA_DIR = config['raw_nanopore_data']
BARCODES_DIR = get_path(config['data'], 'barcodes')

rule merge_barcodes:
    output:
        f"{BARCODES_DIR}/{{barcode}}.fastq.gz"
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