#SOURCES = ["tagesschau"]

SOURCES = ["tagesschau_full"]
# SRCTYPE = ["small"] # ["full"]

rule all:
    input: 
        expand("data/processed/6_mentions_per_week_{src}.json.gz", src=SOURCES),
        expand("data/processed/6_mentions_per_section_per_week_{src}.json.gz", src=SOURCES),
        expand("data/processed/6_mentions_per_section_per_sentiment_per_week_{src}.json.gz", src=SOURCES),
        expand("data/processed/6_topics_{src}.json.gz", src=SOURCES)
        # expand("data/processed/2_sentiment_{src}.json", src=SOURCES),

# $(INTERIM_DATA_DIR)/2_entities_%.jsonl: src/data/fish_entities.py $(RAW_DATA_DIR)/1_scrape_%.jsonl
# 	$(PYTHON_INTERPRETER) $^ $@

# Assume scraped files exist, don't run scrapers by default
# rule run_scrapers:
#     output: expand("data/raw/1_scrape_{src}.jsonl", src=SOURCES)
#     shell: "cd src/data/newsscrape && scrapy crawl tagesspider -O ../../../{output}"

rule preprocess_tagesschau:
    input: "data/raw/1_scrape_tagesschau_{typ}.jsonl"
    output: "data/interim/1_5_preprocessed_tagesschau_{typ}.jsonl"
    script: "src/data/preprocessors/passthrough.py"
    # shell: "python3 src/data/preprocessors/tagesschau.py {input} {output}"

rule preprocess_bild:
    input: "data/raw/1_scrape_bild_{typ}.jsonl"
    output: "data/interim/1_5_preprocessed_bild_{typ}.jsonl"
    script: "src/data/preprocessors/passthrough.py"
    
rule preprocess_zeit:
    input: "data/raw/1_scrape_zeit_{typ}.jsonl"
    output: "data/interim/1_5_preprocessed_zeit_{typ}.jsonl"
    shell: "grep -v '\"sections\": \[\"Seite 1\"\]' {input} > {output}"
    # shell: "python3 src/data/preprocessors/zeit.py {input} {output}"

rule link_entities:
    input: "data/interim/1_5_preprocessed_{source}.jsonl"
    output: "data/interim/2_entities_{source}.jsonl"
    shell: 
        "python3 src/data/fish_entities.py spacyfisher {input} {output}"

rule link_entities_rest:
    input: "data/interim/1_5_preprocessed_{source}.jsonl"
    output: "data/interim/2_entities_rest_{source}.jsonl"
    shell: 
        "python3 src/data/fish_entities.py async-restfisher {input} {output}"

# rule link_entities:
#     input: expand("data/interim/2_entities_{source}.jsonl", source=SOURCES)

rule detect_germanbert_sentiment:
    input: "data/interim/1_5_preprocessed_{source}.jsonl"
    output: "data/interim/2_germanbert_sentiment_{source}.jsonl"
    shell: 
        "python3 src/data/german_bert_sentiment.py {input} {output}"

rule detect_gervader_sentiment:
    input: "data/interim/1_5_preprocessed_{source}.jsonl"
    output: "data/interim/2_vader_sentiment_{source}.jsonl"
    shell: 
        "python3 src/data/gervader_sentiment.py {input} {output}"

rule list_found_entities: 
    input: expand("data/interim/2_entities_{src}.jsonl", src=SOURCES)
    output: "data/interim/3_catalog.json"
    shell: "python3 src/data/list_entities.py {input} {output}"

rule scrape_wikidata:
    input: "data/interim/3_catalog.json"
    output: "data/interim/4_wikidata-catalog.jsonl"
    shell: "cd src/data/wdscrape && scrapy crawl wdspider -a input_file=../../../{input} -a output_file=../../../{output} -O ../../../{output}"

rule merge_wikidata:
    input: "data/interim/2_entities_{source}.jsonl", "data/interim/4_wikidata-catalog.jsonl"
    output: "data/interim/5_{source}.jsonl"
    shell: "python3 src/data/merge_wikidata.py {input} {output}"

rule merge_sentiment:
    input: "data/interim/2_germanbert_sentiment_{source}.jsonl", "data/interim/5_{source}.jsonl"
    output: "data/processed/5_5_{source}.jsonl"
    shell: "python3 src/data/merge_sentiment.py {input} {output}"

rule mentions_per_week:
    input: "data/interim/5_{source}.jsonl"
    output: "data/processed/6_mentions_per_week_{source}.jsonl"
    shell: "python3 src/data/mentions_per_week.py {input} {output}"

rule mentions_per_section_per_week:
    input: "data/interim/5_{source}.jsonl"
    output: "data/processed/6_mentions_per_section_per_week_{source}.jsonl"
    shell: "python3 src/data/mentions_per_section_per_week.py {input} {output}"

rule mentions_per_section_per_sentiment_per_week:
    input: "data/processed/5_5_{source}.jsonl"
    output: "data/processed/6_mentions_per_section_per_sentiment_per_week_{source}.jsonl"
    shell: "python3 src/data/mentions_per_section_per_sentiment_per_week.py {input} {output}"

rule merge_wikidata_2:
    input:
        expand("data/processed/6_mentions_per_week_{src}.jsonl", src=SOURCES), "data/interim/4_wikidata-catalog.jsonl"
    output:
        expand("data/processed/6_5_mentions_per_week_{src}.jsonl", src=SOURCES),
    shell: "python3 src/data/merge_wikidata.py {input} {output}"

rule merge_wikidata_3:
    input:
        expand("data/processed/6_mentions_per_section_per_week_{src}.jsonl", src=SOURCES), "data/interim/4_wikidata-catalog.jsonl"
    output:
        expand("data/processed/6_5_mentions_per_section_per_week_{src}.jsonl", src=SOURCES),
    shell: "python3 src/data/merge_wikidata.py {input} {output}"

rule merge_wikidata_4:
    input:
        expand("data/processed/6_mentions_per_section_per_sentiment_per_week_{src}.jsonl", src=SOURCES), "data/interim/4_wikidata-catalog.jsonl"
    output:
        expand("data/processed/6_5_mentions_per_section_per_sentiment_per_week_{src}.jsonl", src=SOURCES), 
    shell: "python3 src/data/merge_wikidata.py {input} {output}"

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
    shell: "gzip -c {input} > {output}"
