# -*- coding: utf-8 -*-
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv
import pandas as pd
import jsonlines

def merge_wikidata(entities_filepath, wikidata_filepath, output_filepath):
    """ Merges Wikidata information into an entity extract file.
    """
    logger = logging.getLogger(__name__)
    # logger.info('making final data set from raw data')
    df_entities = pd.read_json(open(entities_filepath, "r", encoding="utf8"), lines=True)
    df_wd_catalog = pd.read_json(open(wikidata_filepath, "r", encoding="utf8"), lines=True)
    
    df = df_entities.join(df_wd_catalog.loc[:, ['wd_id', 'wd_label', 'wd_description']].set_index('wd_id'), on='wd_id')

    with jsonlines.open(output_filepath, "w") as f:
        f.write_all(df.to_dict(orient="records"))

log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=log_fmt)

assert len(snakemake.input) == 2
assert len(snakemake.output) == 1

merge_wikidata(snakemake.input[0], snakemake.input[1], snakemake.output[0])