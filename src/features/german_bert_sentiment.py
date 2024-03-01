# -*- coding: utf-8 -*-
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv
from germansentiment import SentimentModel
from tqdm import tqdm
import pandas as pd
import jsonlines


def detect_germanbert_sentiment(input_filepath, output_filepath):
    logger = logging.getLogger(__name__)
    model = SentimentModel("mdraw/german-news-sentiment-bert")

    df_articles = pd.read_json(open(input_filepath, "r", encoding="utf8"), lines=True)

    n = 50
    logger.info(f"analyzing sentiment of {df_articles.shape[0]} articles, partition size: {n}")

    df_articles["shorttext_sentiment"] = "empty"
    df_articles["fulltext_sentiment"] = "empty"
    df_articles_partitions = [
        df_articles[i : i + n] for i in range(0, df_articles.shape[0], n)
    ]

    for partition in tqdm(df_articles_partitions, desc="predicting sentiment", unit="partitions"):
        partition.loc[:, "shorttext_sentiment"] = model.predict_sentiment(
            partition["shorttext"].to_list()
        )
        partition.loc[:, "fulltext_sentiment"] = model.predict_sentiment(
            partition["fulltext"].to_list()
        )

    df_articles = pd.concat(df_articles_partitions)
    df_articles = df_articles.drop(columns=["fulltext", "shorttext", "date", "title"])

    with jsonlines.open(output_filepath, "w") as f:
        f.write_all(df_articles.to_dict("records"))


log_fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
logging.basicConfig(level=logging.INFO, format=log_fmt)

assert len(snakemake.input) == 1
assert len(snakemake.output) == 1

detect_germanbert_sentiment(snakemake.input[0], snakemake.output[0])
