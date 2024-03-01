import scrapy
import json
import jsonlines as jl
from pathlib import Path

class WdspiderSpider(scrapy.Spider):
    name = "wdspider"
    allowed_domains = ["www.wikidata.org"]
    # start_urls = ["https://www.wikidata.org/wiki/Special:EntityData/Q42.json"]


    def start_requests(self):
        existing = []

        existing_entity_ids = set()
        if Path(self.output_file).exists():
            with jl.open(self.output_file, "r") as wd_items:
                existing_entity_ids = set(
                    [
                        item["wd_id"]
                        for item in wd_items.iter(type=dict, skip_invalid=True)
                    ]
                )
        self.logger.info(f"found {len(existing_entity_ids)} existing items")

        with open(self.input_file, "r") as f_in:
            req_entity_ids = json.load(f_in)

        for req_entity_id in req_entity_ids:
            if req_entity_id not in existing_entity_ids:
                yield scrapy.Request(
                    url=f"https://www.wikidata.org/wiki/Special:EntityData/{req_entity_id}.json",
                    callback=self.parse,
                    cb_kwargs=dict(entity_id=req_entity_id),
                )
            else:
                self.logger.info(f"I already know {req_entity_id}, skipping")

    def parse(self, response, entity_id):
        res = json.loads(response.body)
        item = {}
        item["wd_url"] = response.url
        item["wd_id"] = entity_id

        # get correct id if redirected
        id_ = list(res["entities"].keys())[0]
        if id_ != entity_id:
            item["wd_redirect"] = id_

        # try to get the german label first
        if res["entities"][id_]["labels"].get("de"):
            item["wd_label"] = res["entities"][id_]["labels"]["de"]["value"]
        elif res["entities"][id_]["labels"].get("en"):
            item["wd_label"] = res["entities"][id_]["labels"]["en"]["value"]
        else:
            item["wd_label"] = id_

        item["wd_description"] = (
            res["entities"][id_]["descriptions"]["en"]["value"]
            if res["entities"][id_]["descriptions"].get("en")
            else id_
        )

        yield item
