from selenium import webdriver
from selenium.webdriver.chrome.service import Service as ChromeService
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException, NoSuchElementException
from webdriver_manager.chrome import ChromeDriverManager

driverAddress = "/Users/amirho3in/Documents/Stockholm University/Thesis/chromedriver-mac-arm64/chromedriver"

import time
import math
import os
import shutil

import logging 
logging.basicConfig(filename="info.log", 
					format='%(asctime)s %(message)s', 
					filemode='w')

logger=logging.getLogger()
logger.setLevel(logging.ERROR)

download_folder = os.path.abspath("./download")

if not os.path.exists(download_folder):
    os.makedirs(download_folder)

chrome_options = Options()
chrome_options.add_experimental_option("prefs", {
    "download.default_directory": download_folder,
    "profile.content_settings.exceptions.automatic_downloads.*.setting": 1,
    "profile.default_content_settings.popups": 0
})
chrome_options.add_argument("--headless")
chrome_options.add_argument("enable-automation")
chrome_options.add_argument("--no-sandbox")
chrome_options.add_argument("--disable-gpu")

service = ChromeService(service=ChromeService(ChromeDriverManager().install()))
# service = ChromeService(executable_path=driverAddress)
driver = webdriver.Chrome(service=service, options=chrome_options)

try:
    driver.get('https://portal.gdc.cancer.gov/exploration?facetTab=cases&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.primary_site%22%2C%22value%22%3A%5B%22ovary%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.program.name%22%2C%22value%22%3A%5B%22TCGA%22%5D%7D%7D%2C%7B%22content%22%3A%7B%22field%22%3A%22genes.is_cancer_gene_census%22%2C%22value%22%3A%5B%22true%22%5D%7D%2C%22op%22%3A%22in%22%7D%5D%7D&genesTable_size=100&searchTableTab=genes')

    try:
        wait = WebDriverWait(driver, 10)
        accept_button = wait.until(EC.element_to_be_clickable((By.XPATH, "//button[@data-test='modal-cancel-button']")))
        accept_button.click()
        print("Acceptance modal has been acknowledged.")
        logger.info("Acceptance modal has been acknowledged.")
    except TimeoutException:
        print("No acceptance modal appeared.")
        logger.info("No acceptance modal appeared.")

    rows_xpath = "//table[@id='genes-table']/tbody[1]/tr"
    rows = WebDriverWait(driver, 15).until(EC.presence_of_all_elements_located((By.XPATH, rows_xpath)))

    oldElement = driver.find_element(By.XPATH, "//table[@id='genes-table']/tbody[1]/tr[1]/td[2]/a[1]")

    links = {}

    for row in rows:
        title = row.find_element(By.XPATH, './/td[2]/a[1]')
        link = row.find_element(By.XPATH, './/td[5]/span[1]/div[1]/span[1]/span[1]/a[1]')
        links[title.text] = link.get_attribute('href')[:42] + "cases_size=100&" + link.get_attribute('href')[42:]

    next_xpath = "//button[contains(text(),'›')]"
    next = WebDriverWait(driver, 15).until(EC.presence_of_element_located((By.XPATH, next_xpath)))
    driver.execute_script("arguments[0].click();", next)


    WebDriverWait(driver, 10).until(EC.staleness_of(oldElement))

    dropdown_xpath = "//table[@id='genes-table']/tbody[1]/tr"
    rows = WebDriverWait(driver, 10).until(EC.presence_of_all_elements_located((By.XPATH, "//table[@id='genes-table']/tbody[1]/tr")))
    
    for row in rows:
        title = row.find_element(By.XPATH, './/td[2]/a[1]')
        link = row.find_element(By.XPATH, './/td[5]/span[1]/div[1]/span[1]/span[1]/a[1]')
        links[title.text] = link.get_attribute('href')[:42] + "cases_size=100&" + link.get_attribute('href')[42:]

    for link in links:
        new_folder = os.path.abspath("./" + link)
        if not os.path.exists(new_folder):
            os.makedirs(new_folder)
        try:
            driver.get(links[link])
            amount_xpath = "(//html[1]/body[1]/div[2]/div[1]/div[3]/span[1]/div[1]/div[2]/div[1]/span[1]/div[2]/div[2]/div[1]/div[2]/div[1]/div[1]/div[1]/div[1]/strong[3])"
            amount = WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.XPATH, amount_xpath)))
            num = amount.text.replace(",", "")

            pages = math.ceil(int(num) / 100)

            for page in range(1, pages+1):
                try:
                    time.sleep(1)
                    if page == 1:
                        download_button = WebDriverWait(driver, 30).until(EC.presence_of_element_located((By.XPATH, "//button[@class=' button css-1a92y13']")))

                    download_button = driver.find_element(By.XPATH, "//button[@class=' button css-1a92y13']")
                    driver.execute_script("arguments[0].click();", download_button)

                    time.sleep(3)

                    files = os.listdir(download_folder)
                    old_file_name = files[0]
                    new_file_name = str(page) + ".tsv"
                    source_file_path = os.path.join(download_folder, old_file_name)
                    destination_file_path = os.path.join(new_folder, new_file_name)
                    shutil.move(source_file_path, destination_file_path)
                    
                    oldElement = driver.find_element(By.XPATH, "//html[1]/body[1]/div[2]/div[1]/div[3]/span[1]/div[1]/div[2]/div[1]/span[1]/div[2]/div[2]/div[1]/div[2]/div[1]/div[1]/div[2]/table[1]/tbody[1]/tr[1]/td[2]/a[1]")
                    time.sleep(1)

                    if page < pages:
                        next_page_button = driver.find_element(By.XPATH, "//button[text()='›']")
                        driver.execute_script("arguments[0].click();", next_page_button)
                        WebDriverWait(driver, 30).until(EC.staleness_of(oldElement))


                except TimeoutException as e:
                    print(f"{link} - A timeout occurred on page download: {e}")
                    logger.error(f"{link} - A timeout occurred on page download: {e}")
                    break 
        except Exception as e:
            print(f"{link} - An error occurred while processing): {e}")
            logger.error(f"{link} - An error occurred while processing: {e}")
            continue 


except TimeoutException as e:
    print(f"A timeout occurred: {e}")
    logger.error(f"A timeout occurred: {e}")
except NoSuchElementException as e:
    print(f"Element not found: {e}")
    logger.error(f"Element not found: {e}")
finally:
    driver.quit()