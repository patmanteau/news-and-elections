# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv
from germansentiment import SentimentModel
from tqdm import tqdm
import pandas as pd
import jsonlines

@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.argument('output_filepath', type=click.Path())
def main(input_filepath, output_filepath):
    """ Runs data processing scripts to turn raw data from (../raw) into
        cleaned data ready to be analyzed (saved in ../processed).
    """
    logger = logging.getLogger(__name__)
    model = SentimentModel('mdraw/german-news-sentiment-bert')

    df_articles = pd.read_json(open(input_filepath, "r", encoding='utf8'), lines=True)

    n = 50

    df_articles['shorttext_sentiment'] = 'empty'
    df_articles['fulltext_sentiment'] = 'empty'
    df_articles_partitions = [df_articles[i:i+n] for i in range(0, df_articles.shape[0], n)]

    for partition in tqdm(df_articles_partitions):
        partition.loc[:, 'shorttext_sentiment'] = model.predict_sentiment(partition['shorttext'].to_list())
        partition.loc[:, 'fulltext_sentiment'] = model.predict_sentiment(partition['fulltext'].to_list())

    df_articles = pd.concat(df_articles_partitions)
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
