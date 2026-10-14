import pandas as pd
from workflow.lib.utils import get_path

RAW_DATA = get_path(config['data'], 'raw')

SUBUNIT_TO_LINK = {
    "SSU": config['data']['ssu_raw'],
    "LSU": config['data']['lsu_raw']
}

ORGANISM_TO_QUERY = {
    "fungi": "Fungi",
    "trees": '|'.join(config['data']['trees'])
}

rule fetch_raw:
    output:
        f"{RAW_DATA}/{{subunit}}.fasta.gz"
    wildcard_constraints:
        subunit="SSU|LSU"
    threads: 1
    conda: '../envs/raw_data.yml'
    log:
        f"logs/raw_data/fetch_raw/{{subunit}}.log"
    params:
        link=lambda wildcards: SUBUNIT_TO_LINK[wildcards.subunit]
    shell:
        r"""
        wget -O {output} \
            -o {log:q} \
            {params.link}
        """

rule grep_target_organisms:
    input:
        f"{RAW_DATA}/{{subunit}}.fasta.gz"
    output:
        f"{RAW_DATA}/{{subunit}}_{{organism}}.fasta.gz"
    threads: 1
    log:
        f"logs/raw_data/grep_target_organisms/{{subunit}}_{{organism}}.log"
    params:
        query=lambda wc: (
            ";(" + ORGANISM_TO_QUERY[wc.organism] + ");"
        )
    shell:
        r"""
        zcat {input:q} |
            awk -v pattern={params.query:q} '
                /^>/ {{
                    keep = ($0 ~ pattern)
                }}
                keep
            ' |
            gzip -c > {output:q} 2> {log:q}
        """
