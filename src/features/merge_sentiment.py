# -*- coding: utf-8 -*-
import logging
import pandas as pd
import jsonlines

def merge_sentiment(sentiment_filepath, entities_filepath, output_filepath):
    """ Merges Wikidata information into an entity extract file.
    """
    logger = logging.getLogger(__name__)
    # logger.info('making final data set from raw data')
    df_sentiments = pd.read_json(open(sentiment_filepath, "r", encoding="utf8"), lines=True)
    df_entities = pd.read_json(open(entities_filepath, "r", encoding="utf8"), lines=True)
    
    df = df_entities.join(df_sentiments.loc[:, ["url", "shorttext_sentiment", "fulltext_sentiment"]].set_index("url"), on="url")

    with jsonlines.open(output_filepath, "w") as f:
        f.write_all(df.to_dict(orient="records"))

log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=log_fmt)

assert len(snakemake.input) == 2
assert len(snakemake.output) == 1

merge_sentiment(snakemake.input[0], snakemake.input[1], snakemake.output[0])