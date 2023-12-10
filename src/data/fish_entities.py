# -*- coding: utf-8 -*-
import click
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
import logging
import jsonlines as jl
from pathlib import Path
from dotenv import find_dotenv, load_dotenv
import spacy
import pandas as pd

@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.argument('output_filepath', type=click.Path())
def main(input_filepath, output_filepath):
    """ Runs data processing scripts to turn raw data from (../raw) into
        cleaned data ready to be analyzed (saved in ../processed).
    """
    logger = logging.getLogger(__name__)
    logger.info(f'extracting entities from {input_filepath}')

    articles = pd.read_json(open(input_filepath, "r", encoding="utf8"), lines=True)
    # print(articles.sample(5))

    nlp = spacy.load("de_core_news_lg")
    nlp.add_pipe("entityfishing", config={"language": "de", "api_ef_base": "http://127.0.0.1:8090/service"})

    def entity_list(doc) -> list[dict]:
        return [
            {
                "ner_text": ent.text,
                "label": ent.label_,
                "wd_id": ent._.kb_qid,
                "wd_url": ent._.url_wikidata,
                "score": ent._.nerd_score,
            }
            for ent in doc.ents
        ]

    def do_extract_entities(df, partition_info: dict):
        file_name = f"{output_filepath}/2_5_processed-{partition_info['number']}.jsonl"

        with jl.open(file_name, "w") as f_out:
            for text in df.loc[:, ["url", "fulltext"]].itertuples():
                doc = nlp(text.fulltext)
                entities = entity_list(doc)
                for entity in entities:
                    f_out.write(dict(url=text.url, **entity))

    ddf_text = dd.from_pandas(articles, chunksize=200)

    with ProgressBar():
        ddf_text.map_partitions(do_extract_entities, meta=("entities", "object")).compute()

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main()