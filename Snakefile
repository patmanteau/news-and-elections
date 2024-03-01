import os
from dotenv import load_dotenv
load_dotenv()  # take environment variables from .env.

SOURCES = ["tagesschau_full", "zeit_full", "bild_full"]

rule all:
    input: 
        expand("data/processed/6_mentions_per_week_{src}.json.gz", src=SOURCES),
        expand("data/processed/6_mentions_per_section_per_week_{src}.json.gz", src=SOURCES),
        expand("data/processed/6_mentions_per_section_per_sentiment_per_week_{src}.json.gz", src=SOURCES),
        expand("data/processed/6_topics_{src}.json.gz", src=SOURCES)

rule check_download:
    input: expand("data/raw/1_scrape_{src}.jsonl.gz", src=SOURCES)

# Assume scraped files exist, don't run scrapers by default
# rule run_scrapers:
#     output: expand("data/raw/1_scrape_{src}.jsonl", src=SOURCES)
#     shell: "cd src/data/newsscrape && scrapy crawl tagesspider -O ../../../{output}"

# rule ensure_scrapes:
#     output: expand("data/raw/1_scrape_{src}.jsonl", src=SOURCES)
#     run: print(f"One or more of {output} could not be found, please check your Zenodo download.")

rule preprocess_tagesschau:
    input: "data/raw/1_scrape_tagesschau_{typ}.jsonl"
    output: "data/interim/1_5_preprocessed_tagesschau_{typ}.jsonl"
    script: "src/data/preprocessors/passthrough.py"

rule preprocess_bild:
    input: "data/raw/1_scrape_bild_{typ}.jsonl"
    output: "data/interim/1_5_preprocessed_bild_{typ}.jsonl"
    script: "src/data/preprocessors/passthrough.py"
    
rule preprocess_zeit:
    input: "data/raw/1_scrape_zeit_{typ}.jsonl"
    output: "data/interim/1_5_preprocessed_zeit_{typ}.jsonl"
    shell: "grep -v '\"sections\": \[\"Seite 1\"\]' {input} > {output}"
    
rule link_entities:
    input: "data/interim/1_5_preprocessed_{source}.jsonl"
    output: "data/interim/2_entities_{source}.jsonl"
    script: "src/features/link_entities.py"

rule detect_germanbert_sentiment:
    input: "data/interim/1_5_preprocessed_{source}.jsonl"
    output: "data/interim/2_germanbert_sentiment_{source}.jsonl"
    script: "src/features/german_bert_sentiment.py"

rule detect_gervader_sentiment:
    input: "data/interim/1_5_preprocessed_{source}.jsonl"
    output: "data/interim/2_vader_sentiment_{source}.jsonl"
    shell: 
        "python3 src/data/gervader_sentiment.py {input} {output}"

rule list_entities: 
    input: expand("data/interim/2_entities_{src}.jsonl", src=SOURCES)
    output: "data/interim/3_catalog.json"
    script: "src/data/list_entities.py"

rule scrape_wikidata:
    input: "data/interim/3_catalog.json"
    output: "data/interim/4_wikidata-catalog.jsonl"
    shell: "cd src/data/wdscrape && scrapy crawl wdspider -a input_file=../../../{input} -a output_file=../../../{output} -O ../../../{output}"

rule merge_wikidata:
    input: "data/interim/2_entities_{source}.jsonl", "data/interim/4_wikidata-catalog.jsonl"
    output: "data/interim/5_{source}.jsonl"
    script: "src/features/merge_wikidata.py"

rule merge_sentiment:
    input: "data/interim/2_germanbert_sentiment_{source}.jsonl", "data/interim/5_{source}.jsonl"
    output: "data/processed/5_5_{source}.jsonl"
    script: "src/features/merge_sentiment.py"

rule mentions_per_week:
    input: "data/interim/5_{source}.jsonl"
    output: "data/processed/6_mentions_per_week_{source}.jsonl"
    script: "src/features/mentions_per_week.py"

rule mentions_per_section_per_week:
    input: "data/interim/5_{source}.jsonl"
    output: "data/processed/6_mentions_per_section_per_week_{source}.jsonl"
    script: "src/features/mentions_per_section_per_week.py"
    
rule mentions_per_section_per_sentiment_per_week:
    input: "data/processed/5_5_{source}.jsonl"
    output: "data/processed/6_mentions_per_section_per_sentiment_per_week_{source}.jsonl"
    script: "src/features/mentions_per_section_per_sentiment_per_week.py"

rule merge_wikidata_2:
    input:
        "data/processed/6_mentions_per_week_{source}.jsonl", "data/interim/4_wikidata-catalog.jsonl"
    output:
        "data/processed/6_5_mentions_per_week_{source}.jsonl"
    script: "src/features/merge_wikidata.py"

rule merge_wikidata_3:
    input:
        "data/processed/6_mentions_per_section_per_week_{source}.jsonl", "data/interim/4_wikidata-catalog.jsonl"
    output:
        "data/processed/6_5_mentions_per_section_per_week_{source}.jsonl"
    script: "src/features/merge_wikidata.py"

rule merge_wikidata_4:
    input:
        "data/processed/6_mentions_per_section_per_sentiment_per_week_{source}.jsonl", "data/interim/4_wikidata-catalog.jsonl"
    output:
        "data/processed/6_5_mentions_per_section_per_sentiment_per_week_{source}.jsonl"
    script: "src/features/merge_wikidata.py"

rule extract_topics:
    input: "data/interim/1_5_preprocessed_{source}.jsonl"
    output: "data/interim/2_topics_{source}.jsonl"
    shell: "python3 src/data/extract_bertopic.py {input} {output}"

rule merge_topics:
    input: "data/interim/2_topics_{source}.jsonl"
    output: "data/processed/6_topics_{source}.jsonl"
    shell:  "python3 src/data/merge_bertopic.py {input} {output}"

rule jsonl_to_json:
    input: "{filename}.jsonl"
    output: "{filename}.json"
    shell: "jq -s '.' {input} > {output}"

rule gzip:
    input: "{filename}"
    output: "{filename}.gz"
    wildcard_constraints:
        filename=".*(?<!.gz)$" # don't gzip gzipped files to avoid infinite recursion
    shell: "gzip -c {input} > {output}"

# rule gunzip:
#     input: "{filename}.gz"
#     output: "{filename}"
#     wildcard_constraints:
#         filename=".*(?<!.gz)$" # don't gunzip gzipped files to avoid infinite recursion
#     shell: "gzip -dc {input} > {output}"
