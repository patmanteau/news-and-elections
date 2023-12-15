import asyncio
import os
from playwright.async_api import Playwright, async_playwright, expect

from dotenv import load_dotenv

load_dotenv()  # take environment variables from .env.


async def run(playwright: Playwright) -> None:
    browser = await playwright.chromium.launch(headless=False)
    context = await browser.new_context()
    page = await context.new_page()
    await page.goto("https://www.zeit.de/index")
    await page.frame_locator("iframe[title=\"SP Consent Message\"]").get_by_role("button", name="Zustimmen und weiter").click()
    await page.get_by_role("button", name="Nutzermenü").click()
    await page.get_by_role("link", name="Anmelden").click()
    await page.get_by_placeholder("E-Mail-Adresse").click()
    await page.get_by_placeholder("E-Mail-Adresse").fill(os.getenv("ZEIT_USERNAME"))
    await page.get_by_placeholder("E-Mail-Adresse").click()
    await page.get_by_placeholder("Passwort").fill(os.getenv("ZEIT_PASSWORD"))
    await page.get_by_role("button", name="Anmelden").click()

    # ---------------------
    storage = await context.storage_state(path="state.json")

    await context.close()
    await browser.close()


async def main() -> None:
    async with async_playwright() as playwright:
        await run(playwright)


asyncio.run(main())