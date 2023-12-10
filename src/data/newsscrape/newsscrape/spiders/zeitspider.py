import datetime
import scrapy
import os

from dotenv import load_dotenv

load_dotenv()  # take environment variables from .env.

from scrapy_playwright.page import PageMethod
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
        yield scrapy.Request(
            "https://www.zeit.de/zustimmung",
            callback=self.accept_ads,
            meta=dict(
                playwright=True,
                playwright_include_page=True,
            ),
        )
        
        
    async def accept_ads(self, response):
        page = response.meta["playwright_page"]
        await page.frame_locator("iframe[title='SP Consent Message']").get_by_role("button").click()
        await page.screenshot(path="example.png", full_page=True)
        await page.close()

        yield scrapy.Request(
            "https://meine.zeit.de/anmelden",
            callback=self.login,
            meta=dict(
                playwright=True,
                playwright_include_page=True,
            ),
        )


    async def login(self, response):
        username = os.getenv("ZEIT_USERNAME")
        password = os.getenv("ZEIT_PASSWORD")
        if not username and password:
            raise ValueError("Please provide ZEIT_USERNAME and ZEIT_PASSWORD in .env")

        page = response.meta["playwright_page"]
        await page.fill("#login_email", username)
        await page.fill("#login_pass", password)
        await page.click("input[type='submit']")
        await page.close()

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
                    # playwright_page_methods=[
                    #     # PageMethod("click", selector="a[aria-controls^='truncated']")
                    #     PageMethod("click", selector="a[aria-expanded='false']")
                    # ]
                ),
            )


    def parse(self, response):
        # page = response.meta["playwright_page"]
        
        # truncators = page.get_by_text("Weitere Artikel anzeigen").all()
        # while truncators:
        #     for truncator in truncators:
        #         truncator.click()
        #     truncators = page.get_by_text("Weitere Artikel anzeigen").all()
        

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
                    ),
                )
            else:
                item["fulltext"] = None
                item["date"] = None
                yield item

            # item["tstamp"] = article.css("::attr(data-teaserdate)").get()
            # item["title"] = article.css(".teaser-right__headline::text").get().strip()
            # # naive datetime string
            # item["date"] = article.css(".teaser-right__date::text").get().strip()
            # item["shorttext"] = (
            #     article.css(".teaser-right__shorttext::text").get().strip()
            # )

            # # try to get the article's url and scrape its full text
            # item["url"] = article.css(".teaser-right__link::attr(href)").get()
            # if item["url"]:
            #     yield scrapy.Request(
            #         f'https://www.tagesschau.de{item["url"]}',
            #         self.parse_full_text,
            #         cb_kwargs={"item": item},
            #     )
            # else:
            #     item["fulltext"] = None
            #     yield item

            yield item

    def parse_full_item(self, response, item):
        item["fulltext"] = " ".join(
            [p.strip() for p in response.css("article > p").css("::text").getall()]
        )
        yield item
