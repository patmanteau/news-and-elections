import click
import datetime
import jsonlines as jl
from bs4 import BeautifulSoup
import time
from playwright.sync_api import sync_playwright, BrowserContext
from tqdm import tqdm
# from dotenv import load_dotenv

# load_dotenv()  # take environment variables from .env.

# Constants for the website, login details, and throttle time
# BASE_URL = "https://www.zeit.de/issue/{issue}"
# LOGIN_URL = "https://meine.zeit.de/anmelden"
# SUCCESS_URL = "https://www.zeit.de/konto"
# USERNAME = os.environ["ZEIT_USERNAME"]
# PASSWORD = os.environ["ZEIT_PASSWORD"]
# THROTTLE_TIME = 2  # seconds between requests

# HEADERS = {
#     "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.2.1 Safari/605.1.15"
# }


# # Function to handle authentication
# def authenticate(session):
#     # Find csrf_token
#     response = session.get(LOGIN_URL)
#     login_page = BeautifulSoup(response.text)

#     csrf_token = None
#     for inp in login_page.find_all("input"):
#         if inp.has_attr("name") and inp["name"] == "csrf_token":
#             csrf_token = inp["value"]

#     # print(login_page)
#     print(f"{csrf_token=}")
#     print(response)

#     payload = {
#         "email": USERNAME,
#         "pass": PASSWORD,
#         "entry_service": "sonstige",
#         "product_id": "sonstige",
#         "return_url": None,
#         "permanent": True,
#         "csrf_token": csrf_token,
#     }
#     print(payload)
#     response = session.post(LOGIN_URL, data=payload, headers=HEADERS)

#     print(response.text)
#     if response.url == SUCCESS_URL:
#         print("Login successful!")
#     else:
#         print("Login failed.")


# # Function to throttle requests
# def throttle():
#     time.sleep(THROTTLE_TIME)


@click.command
@click.argument("start_year", type=int)
@click.argument("end_year", type=int)
@click.argument("start_issue", type=int)
@click.argument("end_issue", type=int)
def list_(
    start_year: int, end_year: int, start_issue: int, end_issue: int
) -> list[str]:
    issues = list_issues(start_year, end_year, start_issue, end_issue)
    print(f"{issues=}")

    articles = []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        context = browser.new_context()

        page = context.new_page()
        page.goto("https://www.zeit.de/index")
        page.frame_locator('iframe[title="SP Consent Message"]').get_by_role(
            "button", name="Zustimmen und weiter"
        ).click()

        print("Listing issues...")
        with jl.open("zeit_article_list.jsonl", "w") as f:
            for issue in tqdm(issues):
                listed_articles = list_articles(context, issue)
                articles.extend(listed_articles)
                f.write_all(listed_articles)

        context.close()
        browser.close()

    return articles


def list_issues(
    start_year: int, end_year: int, start_issue: int, end_issue: int
) -> list[str]:
    urls = []
    year = start_year
    issue = start_issue
    
    assert start_year <= end_year

    while year < end_year or (issue <= end_issue):
        urls.append(f"https://www.zeit.de/{year}/{issue:02}/index")
        issue += 1
        if issue == 56:
            year += 1
            issue = 1

    return urls


def scrape_article(context: BrowserContext, article: dict) -> dict:
    page = context.new_page()

    retries = 0
    while retries < 10:
        try:
            response = page.goto(article.get("url"))
            if not response.ok:
                raise ValueError("not response.ok")
            break
        except:
            retries += 1
            time.sleep(1)

    article_page = BeautifulSoup(page.content(), features="lxml")

    # Komplettansicht available?
    full_page = article_page.select_one("a[data-ct-label='all']")
    if full_page:
        new_url = full_page["href"]
        article["url"] = new_url
        return scrape_article(context, article)
    else:
        paragraphs = [
            " ".join(p.stripped_strings)
            for p in article_page.select("p.paragraph.article__item")
        ]

        date = article_page.select_one("time")["datetime"]
        tstamp = datetime.datetime.fromisoformat(date).timestamp()
        title = article_page.select_one("span.article-heading__title")
        scraped = {
            "tstamp": str(int(tstamp)),
            "title": title.string if title else "",
            "date": date,
            "shorttext": " ".join(paragraphs),
            "url": article.get("url"),
            "sections": article.get("sections"),
            "fulltext": " ".join(paragraphs),
        }
        return scraped


def list_articles(context: BrowserContext, issue_url: str) -> list[str]:
    page = context.new_page()
    response = page.goto(issue_url)

    if not response.ok:
        return []
    else:
        issue_page = BeautifulSoup(page.content(), features="lxml")
        articles = []

        # Find all articles and section headers, in order
        def is_article_or_section_header(element):
            if element.name not in ["article", "h2"]:
                return False
            elif element.name == "article":
                return True
            else:
                print(f"{element=}")
                return element.get("class") == "cp-area__headline"
            
        articles_and_headers = issue_page.find_all(["article", "h2"])
        current_section = ""
        for element in articles_and_headers:
            match element.name:
                case "h2":
                    current_section = element.string
                case "article":
                    if current_section != "Ressorts":
                        articles.append({
                            "url": element.a["href"],
                            "sections": [current_section] if current_section else "",
                        })
        return articles
        

@click.command
def scrape():
    articles = []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        context = browser.new_context()

        page = context.new_page()
        page.goto("https://www.zeit.de/index")
        page.frame_locator('iframe[title="SP Consent Message"]').get_by_role(
            "button", name="Zustimmen und weiter"
        ).click()

        context.close()
        context = browser.new_context(storage_state="newsscrape/state.json")

        # load article urls from file and drop sub ads
        with jl.open("zeit_article_list.jsonl") as f_in:
            articles = list([a for a in f_in if a.get("url").find("premium.zeit.de") == -1])

        count = 0
        with jl.open("zeit_articles.jsonl", "w", flush=True) as f_out:
            print("Scraping articles...")
            for article in tqdm(articles):
                try:
                    scraped = scrape_article(context, article)
                except:
                    print(f"Error scraping {article.get('url')}")
                    continue
                # print(scraped)
                f_out.write(scraped)
                # restart browser every 200 articles to avoid running out of memory
                count += 1
                if count > 200:
                    context.close()
                    browser.close()
                    browser = playwright.chromium.launch(headless=True)
                    context = browser.new_context(storage_state="newsscrape/state.json")
                    count = 0
                
        context.close()
        browser.close()
    

@click.group()
def cli():
    pass


cli.add_command(scrape)
cli.add_command(list_)

if __name__ == "__main__":
    cli()
