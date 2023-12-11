# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv
import pandas as pd
import json

@click.command()
@click.argument('input_files', nargs=-1, type=click.Path(exists=True))
@click.argument('output_filepath', nargs=1, type=click.Path())
def main(input_files, output_filepath):
    """ Collects Wikidata entity IDs for all entities in the given input files and
        writes them to the output file.
    """
    logger = logging.getLogger(__name__)
    
    dataframes = []
    for input_file in input_files:
        logger.info(f'collecting entities from {input_file}')
        dataframes.append(pd.read_json(open(input_file, "r", encoding="utf8"), lines=True))

    df = pd.concat(dataframes)
    with open(output_filepath, "w") as f_out:
            json.dump(df["wd_id"].sort_values().drop_duplicates().to_list(), f_out)


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main()
