# -*- coding: utf-8 -*-
import aiohttp
import asyncio
import click
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
import logging
import jsonlines as jl
from pathlib import Path
from dotenv import find_dotenv, load_dotenv
import spacy
import pandas as pd
import numpy as np
import tempfile
import shutil
import requests

def chunked_iterable(iterable, chunk_size):
    """Yield successive chunks from iterable."""
    for i in range(0, len(iterable), chunk_size):
        yield iterable[i:i + chunk_size]

async def call_disambiguate(session, text, url, tstamp, sections):
    req = {
        "text": text,
        "shortText": "",
        "termVector": [],
        "language": {
            "lang": "de"
        },
        "entities": [],
        "mentions": [
            "ner",
            "wikipedia"
        ],
        "nbest": False,
        "sentence": False,
        "minSelectorScore": 0.3,
    }

    async with session.post("http://127.0.0.1:8090/service/disambiguate", json=req) as response:
        doc = await response.json()
        entities = []
        for entity in doc.get("entities", []):
            entities.append(dict(
                url=url,
                tstamp=tstamp,
                sections=sections,
                ner_text=entity.get("text"),
                label="",
                wd_id=entity.get("wikidataId"),
                wd_url=f"https://www.wikidata.org/wiki/{entity.get('wikidataId')}",
                score=entity.get("confidence_score"),
            ))
        return entities
    
async def do_async_restfisher(input_filepath, output_filepath):
    logger = logging.getLogger(__name__)
    logger.info(f'extracting entities from {input_filepath}')

    articles = pd.read_json(open(input_filepath, "r", encoding="utf8"), lines=True)

    # write partitioned jsonl files into a temp dir
    with tempfile.TemporaryDirectory() as tmpdirname:
        chunks = chunked_iterable(articles, 200)
        for chunk_number, chunk in enumerate(chunks):
            async with aiohttp.ClientSession() as session:
                file_name = f"{tmpdirname}/entities-{chunk_number}.jsonl"
                with jl.open(file_name, "w") as f_out:
                    logger.info(f'working on chunk {chunk_number} and writing to {file_name}')
                    tasks = []
                    for text in chunk.loc[:, ["tstamp", "url", "fulltext", "sections"]].itertuples():
                        tasks.append(asyncio.ensure_future(call_disambiguate(session, text.fulltext, text.url, text.tstamp, text.sections)))
                    logger.info(f'waiting for {len(tasks)} tasks to complete')
                    entities = await asyncio.gather(*tasks)
                    f_out.write_all(entities)
                        
        # read the partitioned jsonl files and write them into a single jsonl file
        with open(output_filepath, "wb") as f_out:
            for chunk in Path(tmpdirname).glob("*.jsonl"):
                with open(chunk,"rb") as f_in:
                    shutil.copyfileobj(f_in, f_out)

@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.argument('output_filepath', type=click.Path())
def async_restfisher(input_filepath, output_filepath):
    asyncio.run(do_async_restfisher(input_filepath, output_filepath))

@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.argument('output_filepath', type=click.Path())
def restfisher(input_filepath, output_filepath):
    logger = logging.getLogger(__name__)
    logger.info(f'extracting entities from {input_filepath}')

    articles = pd.read_json(open(input_filepath, "r", encoding="utf8"), lines=True)

    # write partitioned jsonl files into a temp dir
    with tempfile.TemporaryDirectory() as tmpdirname:
        logger.info(f'working directory is {tmpdirname}')
        
        def do_extract_entities(df, partition_info: dict):
            file_name = f"{tmpdirname}/entities-{partition_info['number']}.jsonl"

            with jl.open(file_name, "w") as f_out:
                for text in df.loc[:, ["tstamp", "url", "fulltext", "sections"]].itertuples():
                    req = {
                        "text": text.fulltext,
                        "shortText": "",
                        "termVector": [],
                        "language": {
                            "lang": "de"
                        },
                        "entities": [],
                        "mentions": [
                            "ner",
                            "wikipedia"
                        ],
                        "nbest": False,
                        "sentence": False,
                        "minSelectorScore": 0.3,
                    }
                    response = requests.post("http://127.0.0.1:8090/service/disambiguate", json=req)
                    if response.ok:
                        doc = response.json()
                        entities = doc.get("entities", [])
                        for entity in entities:
                            f_out.write(dict(
                                url=text.url,
                                tstamp=text.tstamp,
                                sections=text.sections,
                                ner_text=entity.get("text"),
                                label="",
                                wd_id=entity.get("wikidataId"),
                                wd_url=f"https://www.wikidata.org/wiki/{entity.get('wikidataId')}",
                                score=entity.get("confidence_score"),
                            ))

        ddf_text = dd.from_pandas(articles, chunksize=200)

        with ProgressBar():
            ddf_text.map_partitions(do_extract_entities, meta=("entities", "object")).compute()

        # read the partitioned jsonl files and write them into a single jsonl file
        with open(output_filepath, "wb") as f_out:
            for chunk in Path(tmpdirname).glob("*.jsonl"):
                with open(chunk,"rb") as f_in:
                    shutil.copyfileobj(f_in, f_out)


@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.argument('output_filepath', type=click.Path())
def spacyfisher(input_filepath, output_filepath):
    """ Runs data processing scripts to turn raw data from (../raw) into
        cleaned data ready to be analyzed (saved in ../processed).
    """
    logger = logging.getLogger(__name__)
    logger.info(f'extracting entities from {input_filepath}')

    articles = pd.read_json(open(input_filepath, "r", encoding="utf8"), lines=True)

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

@click.group()
def cli():
    pass


cli.add_command(spacyfisher)
# cli.add_command(restfisher)
# cli.add_command(async_restfisher)


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    cli()
