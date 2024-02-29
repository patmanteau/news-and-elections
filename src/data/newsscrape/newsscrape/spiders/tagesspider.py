import scrapy

# from scrapy.linkextractors import LinkExtractor
from scrapy.spiders import CrawlSpider, Rule
from scrapy.loader import ItemLoader
import datetime

from ..items import NewsscrapeItem


class TagesspiderSpider(CrawlSpider):
    name = "tagesspider"
    allowed_domains = ["www.tagesschau.de"]
    # start_urls = ["https://www.tagesschau.de/archiv?datum=2023-06-22"]
    # start_urls = list(urls())

    def start_requests(self):
        # base = datetime.datetime.today()
        # start date on commandline is mandatory, ISO format
        # end date on commandline is optional, default is today
        if hasattr(self, "end_date"):
            base = datetime.datetime.strptime(self.end_date, "%Y-%m-%d")
        else:
            base = datetime.datetime.today()
        td = base - datetime.datetime.strptime(self.start_date, "%Y-%m-%d")
        date_list = reversed(
            [base - datetime.timedelta(days=x) for x in range(td.days)]
        )
        for date in date_list:
            yield scrapy.Request(
                url=f"https://www.tagesschau.de/archiv?datum={date.strftime('%Y-%m-%d')}",
                callback=self.parse,
            )

    def parse(self, response):
        next_ = response.css("li.next")
        if next_:
            # print(next_)
            next_page = next_.css(".paginierung__liste--link::attr(href)").get()
            if next_page:
                yield scrapy.Request(
                    url=f"https://www.tagesschau.de/archiv/{next_page}",
                    callback=self.parse
                )
        articles = response.css(".teaser-right")
        for article in articles:
            item = NewsscrapeItem()
            # item = {}
            # unix timestamp?
            item["tstamp"] = article.css("::attr(data-teaserdate)").get()
            item["title"] = article.css(".teaser-right__headline::text").get().strip()
            # naive datetime string
            item["date"] = article.css(".teaser-right__date::text").get().strip()
            item["shorttext"] = (
                article.css(".teaser-right__shorttext::text").get().strip()
            )

            # try to get the article's url and scrape its full text
            item["url"] = article.css(".teaser-right__link::attr(href)").get()
            if item["url"]:
                item["sections"] = item["url"][1:].split('/')[:-1]
                yield scrapy.Request(
                    f'https://www.tagesschau.de{item["url"]}',
                    self.parse_full_text,
                    cb_kwargs={"item": item},
                )
            else:
                item["fulltext"] = None
                item["sections"] = []
                yield item

    def parse_full_text(self, response, item):
        item["fulltext"] = " ".join(
            [p.strip() for p in response.css("article > p").css("::text").getall()]
        )
        yield item
