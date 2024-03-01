# -*- coding: utf-8 -*-
import logging
import pandas as pd


def mentions_per_section_per_sentiment_per_week(input_filepath, output_filepath):
    """ Calculates the number of mentions per week for each entity.
    """
    logger = logging.getLogger(__name__)

    df = pd.read_json(open(input_filepath, "r", encoding="utf8"), lines=True)
    df['tstamp'] = df['tstamp'].apply(lambda x: pd.Timestamp(x, unit='s', tz='Europe/Berlin'))
    df['section'] = df['sections'].apply(lambda x: x[0] if x else None)

    # Group by day and entity, then count occurrences
    mentions = df.groupby([pd.Grouper(key='tstamp', freq='W'), 'section', 'shorttext_sentiment', 'wd_id']).size().reset_index(name='count')
    mentions.sort_values(by='tstamp', inplace=True)
    mentions.rename(columns={'shorttext_sentiment': 'sentiment'}, inplace=True)
    mentions.to_json(output_filepath, orient='records', date_format='iso', lines=True)
    # with jsonlines.open(output_filepath, "w") as f:
    #     f.write_all(mentions.to_dict(orient="records"))


log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=log_fmt)

assert len(snakemake.input) == 1
assert len(snakemake.output) == 1

mentions_per_section_per_sentiment_per_week(snakemake.input[0], snakemake.output[0])