# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv
import pandas as pd
import jsonlines

@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.argument('output_filepath', type=click.Path())
def main(input_filepath, output_filepath):
    """ Calculates the number of mentions per week for each entity.
    """
    logger = logging.getLogger(__name__)

    df = pd.read_json(open(input_filepath, "r", encoding="utf8"), lines=True)
    df['tstamp'] = df['tstamp'].apply(lambda x: pd.Timestamp(x, unit='s', tz='Europe/Berlin'))

    # Group by day and entity, then count occurrences
    mentions = df.groupby([pd.Grouper(key='tstamp', freq='W'), 'wd_id']).size().reset_index(name='count')
    mentions.sort_values(by='tstamp', inplace=True)
    mentions.to_json(output_filepath, orient="records", date_format='iso', lines=True)
    # with jsonlines.open(output_filepath, "w") as f:
    #     f.write_all(mentions.to_dict(orient="records"))


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main()
