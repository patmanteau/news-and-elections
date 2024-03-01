# -*- coding: utf-8 -*-
import logging
from pathlib import Path
import pandas as pd
import json

def list_entities(input_files, output_file):
    """ Collects Wikidata entity IDs for all entities in the given input files and
        writes them to the output file.
    """
    logger = logging.getLogger(__name__)
    
    dataframes = []
    for input_file in input_files:
        logger.info(f'reading entities from {input_file}')
        dataframes.append(pd.read_json(open(input_file, "r", encoding="utf8"), lines=True))

    df = pd.concat(dataframes)
    logger.info(f'sorting, deduplicating and writing entities to {output_file}')
    with open(output_file, "w") as f_out:
            json.dump(df["wd_id"].sort_values().drop_duplicates().to_list(), f_out)


log_fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
logging.basicConfig(level=logging.INFO, format=log_fmt)

assert len(snakemake.input) >= 1
assert len(snakemake.output) == 1

list_entities(snakemake.input, snakemake.output[0])
