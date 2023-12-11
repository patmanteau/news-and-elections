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
import tempfile
import shutil

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

    # write partitioned jsonl files into a temp dir
    with tempfile.TemporaryDirectory() as tmpdirname:
        logger.info(f'working directory is {tmpdirname}')
        
        def do_extract_entities(df, partition_info: dict):
            file_name = f"{tmpdirname}/entities-{partition_info['number']}.jsonl"

            with jl.open(file_name, "w") as f_out:
                for text in df.loc[:, ["tstamp", "url", "fulltext", "sections"]].itertuples():
                    doc = nlp(text.fulltext)
                    entities = entity_list(doc)
                    for entity in entities:
                        f_out.write(dict(url=text.url, tstamp=text.tstamp, sections=text.sections, **entity))

        ddf_text = dd.from_pandas(articles, chunksize=200)

        with ProgressBar():
            ddf_text.map_partitions(do_extract_entities, meta=("entities", "object")).compute()

        # read the partitioned jsonl files and write them into a single jsonl file
        with open(output_filepath, "wb") as f_out:
            for chunk in Path(tmpdirname).glob("*.jsonl"):
                with open(chunk,"rb") as f_in:
                    shutil.copyfileobj(f_in, f_out)

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main()
