# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv
from tqdm import tqdm
import pandas as pd
import jsonlines
import sys

# sys.path.append('../models/GerVADER/')
from vaderSentimentGER import SentimentIntensityAnalyzer

@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.argument('output_filepath', type=click.Path())
def main(input_filepath, output_filepath):
    """ Runs data processing scripts to turn raw data from (../raw) into
        cleaned data ready to be analyzed (saved in ../processed).
    """
    logger = logging.getLogger(__name__)

    df_articles = pd.read_json(open(input_filepath, "r", encoding='utf8'), lines=True)
    df_articles['fulltext_sentiment'] = 'empty'
    df_articles['shorttext_sentiment'] = 'empty'

    analyzer = SentimentIntensityAnalyzer()

    def classification(vs: dict) -> str:
        if vs["compound"] >= 0.05:
            return 'positive'
        elif vs["compound"] <= -0.05:
            return 'negative'
        else:
            return 'neutral'

    # Create new `pandas` methods which use `tqdm` progress
    # (can use tqdm_gui, optional kwargs, etc.)
    tqdm.pandas()

    #df_articles['fulltext_sentiment'] =  df_articles.progress_apply(lambda row: classification(analyzer.polarity_scores(row['fulltext'])), axis=1)
    df_articles['shorttext_sentiment'] = df_articles.progress_apply(lambda row: classification(analyzer.polarity_scores(row['shorttext'])), axis=1)

    df_articles = df_articles.drop(columns=['fulltext', 'shorttext', 'date', 'title'])

    with jsonlines.open(output_filepath, 'w') as f:
        f.write_all(df_articles.to_dict('records'))


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main()
