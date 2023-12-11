SOURCES = ["tagesschau", "zeit"]

rule all:
    input: 
        expand("data/processed/6_mentions_per_week_{src}.json", src=SOURCES),
        expand("data/processed/6_mentions_per_section_per_week_{src}.json", src=SOURCES),

# $(INTERIM_DATA_DIR)/2_entities_%.jsonl: src/data/fish_entities.py $(RAW_DATA_DIR)/1_scrape_%.jsonl
# 	$(PYTHON_INTERPRETER) $^ $@

# Assume scraped files exist, don't run scrapers by default
# rule run_scrapers:
#     output: expand("data/raw/1_scrape_{src}.jsonl", src=SOURCES)
#     shell: "cd src/data/newsscrape && scrapy crawl tagesspider -O ../../../{output}"

rule link_entities:
    input: "data/raw/1_scrape_{source}.jsonl"
    output: "data/interim/2_entities_{source}.jsonl"
    shell: 
        "python3 src/data/fish_entities.py {input} {output}"

# rule link_entities:
#     input: expand("data/interim/2_entities_{source}.jsonl", source=SOURCES)

rule catalog: 
    input: expand("data/interim/2_entities_{src}.jsonl", src=SOURCES)
    output: "data/interim/3_catalog.json"
    shell: "python3 src/data/list_entities.py {input} {output}"

rule wikidata_catalog:
    input: "data/interim/3_catalog.json"
    output: "data/interim/4_wikidata-catalog.jsonl"
    shell: "cd src/data/wdscrape && scrapy crawl wdspider -a input_file=../../../{input} -a output_file=../../../{output} -O ../../../{output}"

rule merge_wikidata:
    input: "data/interim/2_entities_{source}.jsonl", "data/interim/4_wikidata-catalog.jsonl"
    output: "data/processed/5_{source}.jsonl"
    shell: "python3 src/data/merge_wikidata.py {input} {output}"

rule mentions_per_week:
    input: "data/processed/5_{source}.jsonl"
    output: "data/processed/6_mentions_per_week_{source}.json"
    shell: "python3 src/data/mentions_per_week.py {input} {output}"

rule mentions_per_section_per_week:
    input: "data/processed/5_{source}.jsonl"
    output: "data/processed/6_mentions_per_section_per_week_{source}.json"
    shell: "python3 src/data/mentions_per_section_per_week.py {input} {output}"
