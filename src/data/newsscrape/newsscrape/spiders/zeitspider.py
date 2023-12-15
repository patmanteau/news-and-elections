import datetime
import scrapy
import os

from hashlib import blake2b
h = blake2b(digest_size=10)

from dotenv import load_dotenv

load_dotenv()  # take environment variables from .env.

from scrapy_playwright.page import PageMethod
from playwright.sync_api import expect 
from ..items import NewsscrapeItem

class ZeitspiderSpider(scrapy.Spider):
    name = "zeitspider"
    allowed_domains = ["www.zeit.de", "meine.zeit.de"]
    # start_urls = ["https://www.zeit.de"]
    start_year = None
    end_year = None
    start_issue = None
    end_issue = None

    @classmethod
    def update_settings(cls, settings):
        super().update_settings(settings)
        settings.set(
            "DOWNLOAD_HANDLERS",
            {
                "http": "scrapy_playwright.handler.ScrapyPlaywrightDownloadHandler",
                "https": "scrapy_playwright.handler.ScrapyPlaywrightDownloadHandler",
            }, 
            priority="spider"
        )

    def __init__(self, start_year: int = None, end_year: int = None, start_issue: int = None, end_issue: int = None, *args, **kwargs):
        super(ZeitspiderSpider, self).__init__(*args, **kwargs)

        if not all((start_year, end_year, start_issue, end_issue)):
            raise ValueError("Please provide start_year, end_year, start_issue and end_issue")

        self.start_year = int(start_year)
        self.end_year = int(end_year)
        self.start_issue = int(start_issue)
        self.end_issue = int(end_issue)


    def start_requests(self):
        year = self.start_year
        issue = self.start_issue
        urls = []

        while year <= self.end_year and issue <= self.end_issue:
            urls.append(f"https://www.zeit.de/{year}/{issue}/index")
            issue += 1
            if issue == 56:
                year += 1
                issue = 1

        for url in urls:
            yield scrapy.Request(
                url=url,
                callback=self.parse,
                meta=dict(
                    playwright=True,
                    # playwright_include_page=True,
                    playwright_context="auth",
                    playwright_context_kwargs={
                        "storage_state": "state.json"
                    },
                    # playwright_page_methods=[
                    #     PageMethod('wait_for_load_state', state = 'domcontentloaded'),
                    #     PageMethod("evaluate", "window.scrollBy(0, document.body.scrollHeight)"),
                    # ]                ),
                ),
            )

        
    
    def parse(self, response):
        # page = response.meta["playwright_page"]
        
        articles = response.css("article")

        for article in articles:
            item = NewsscrapeItem()
            
            large_title = article.css(".teaser-large__title::text").get()
            small_title = article.css(".teaser-small__title::text").get()

            if large_title:
                item["title"] = large_title.strip()
            elif small_title:
                item["title"] = small_title.strip()
            else:
                item["title"] = None

            item["shorttext"] = article.css("div > p::text").get()
            if item["shorttext"]:
                item["shorttext"] = item["shorttext"].strip()
            
            item["url"] = article.css("a::attr(href)").get()
            if item["url"]:
                yield scrapy.Request(
                    item["url"],
                    self.parse_full_item,
                    cb_kwargs={"item": item},
                    meta=dict(
                        playwright=True,
                        playwright_include_page=True,
                        # playwright_page=page,
                        playwright_context="auth",
                        playwright_context_kwargs={
                            "storage_state": "state.json"
                        },
                        errback=self.errback_close_page,
                        # playwright_page_methods=[
                        #     PageMethod('wait_for_load_state', state = 'domcontentloaded'),
                        #     PageMethod("evaluate", "window.scrollBy(0, document.body.scrollHeight)"),
                        # ],
                    ),
                )
            else:
                item["fulltext"] = ""
                item["date"] = ""
                
                yield item


    async def parse_full_item(self, response, item):
        page = response.meta["playwright_page"]
        # h.update(response.url.encode("utf-8"))
        # urlhash = h.hexdigest()
        # await page.screenshot(path=f"{urlhash}.png", full_page=True)

        komplettansicht_url = response.xpath("//a[@class='article-toc__fullview']/@href").get()
        self.logger.info(f"Komplettansicht URL: {komplettansicht_url}")
        
        
        if komplettansicht_url:
            yield scrapy.Request(
                komplettansicht_url,
                self.parse_full_item,
                cb_kwargs={"item": item},
                meta=dict(
                    playwright=True,
                    playwright_include_page=True,
                    # playwright_page=page,
                    playwright_context="auth",
                    playwright_context_kwargs={
                        "storage_state": "state.json"
                    },
                    errback=self.errback_close_page,
                    # playwright_page_methods=[
                    #     PageMethod('wait_for_load_state', state = 'domcontentloaded'),
                    #     PageMethod("evaluate", "window.scrollBy(0, document.body.scrollHeight)"),
                    # ],
                ),
            )
            await page.close()
        else:
            item["fulltext"] = " ".join( 
                [p.strip() for p in await page.locator(".paragraph").all_inner_texts()]
            )
            item["tstamp"] = await page.locator(".metadata__date > time").get_attribute("datetime")
                 
            yield item
            await page.close()

    async def errback_close_page(self, failure):
        page = failure.request.meta["playwright_page"]
        await page.close()