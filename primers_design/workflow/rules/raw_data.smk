import pandas as pd
from workflow.lib.utils import get_path

RAW_DATA = get_path(config['data'], 'raw')

SUBUNIT_TO_LINK = {
    "SSU": config['data']['ssu_raw'],
    "LSU": config['data']['lsu_raw']
}

rule fetch_raw:
    output:
        f"{RAW_DATA}/{{subunit}}.fasta.gz"
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