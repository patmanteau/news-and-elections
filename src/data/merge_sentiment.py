# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv
import pandas as pd
import jsonlines

@click.command()
@click.argument('sentiment_filepath', type=click.Path(exists=True))
@click.argument('entities_filepath', type=click.Path(exists=True))
@click.argument('output_filepath', type=click.Path())
def main(sentiment_filepath, entities_filepath, output_filepath):
    """ Merges Wikidata information into an entity extract file.
    """
    logger = logging.getLogger(__name__)
    # logger.info('making final data set from raw data')
    df_sentiments = pd.read_json(open(sentiment_filepath, "r", encoding="utf8"), lines=True)
    df_entities = pd.read_json(open(entities_filepath, "r", encoding="utf8"), lines=True)
    
    df = df_entities.join(df_sentiments.loc[:, ["url", "shorttext_sentiment", "fulltext_sentiment"]].set_index("url"), on="url")

    with jsonlines.open(output_filepath, "w") as f:
        f.write_all(df.to_dict(orient="records"))

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main()
