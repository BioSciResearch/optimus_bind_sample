'''
Example: https://github.com/benlindsay/nfl-dash/blob/master/src/data/make_raw_data.py
'''

# Scrape
import requests
import bs4
import lxml.etree as xml
import re

# Download
from urllib.request import urlretrieve
import sys

# import zipfile
import os

# move to top
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv


def skempi_scraper():
    ''' Scrapes SKEMPI download page for the purposes of:
          • Finding path of downloads
          • Checking for updates
          • Downloading old versions for reproducability
        probably overkill, but I'm putting it forth as an option.

        Output:
            current = largest key in versions dictionary
            versions = {ver:[csv link, pdb link], ...}
                        download links

        Notes:
          • url to SKEMPI page is hardcoded
          • v1.1 is html, not a downloadable csv
          • v1.0 & 1.1 have no tgz of pdb, all included in v2+
    '''
    # URLs
    home = 'https://life.bsc.es'
    URL = home+'/pid/skempi2/database/index'

    # Get HTML: header prevents getting blacklisted
    webPage = bs4.BeautifulSoup(requests.get(
        URL,
        stream=True,
        headers={"UserAgent": "Mozilla/5.0 (X11; Linux x86_64)" +
                 "AppleWebKit/537.36 (KHTML, like Gecko)" +
                 "Chrome/65.0.3325.183 Safari/537.36"}).text,
        "lxml")

    # Exctract potential download links:
    links = webPage.findAll(name='a', attrs={"href": re.compile(".csv|tgz")})

    # Find most recent verion:
    versions = dict()  # format = {ver:[csv link, pdb link], ...}
    for a in links:
        ver = eval(re.match('SKEMPI v[.]?(\d+\.\d+)', a.text)[1])
        href = a.get('href')
        if ver not in versions:
            versions[ver] = [None, None]
        versions[ver]['tgz' in href] = home+href if href[0] == '/' else href
    current = max(versions.keys())
    return current, versions


def progress_hook(downloaded, chunk_size, total_size):
    '''A hook to report int–percent progress of a downloads.
    Input:
      • downloaded: number of data chunks downloaded so far
      • chunk_size: bytes per chunk
      • total_size: total bytes of file

    Output:
      • Terminal response of download progress

    Code modified from:
      • https://stackoverflow.com/questions/2028517/python-urllib2-progress-hook
      • https://github.com/jfconavarrete/kaggle-facebook/blob/master/src/data/download_dataset.py
      • https://blog.shichao.io/2012/10/04/progress_speed_indicator_for_urlretrieve_in_python.html
    '''
    percent = downloaded * chunk_size * 100 // total_size  # %chunks downloaded
    total = round(total_size/chunk_size)  # total number of chunks
    response = f"Downloaded {downloaded} of {total} bytes {percent}%\r"
    sys.stdout.write(response)
    if downloaded >= total:
        sys.stdout.write('\n')
    sys.stdout.flush()


def skempi_download(v=None):
    '''MESSY, just wanted to get it working...

        Input:
          • v = choice of version ['1.0','2.0']
          * vkey = choice of version ['1.0','2.0']
          * versions or URLs
          * download path?

        Created:
          • raw/Skempi_{v}.csv
          • raw/Skempi_{v}.tgz
          • raw/PDBs/ .pdb and .mapping files

        Outputs:
          • progress reports
    '''
    # Init variables
    # ADD AS INPUTS WHEN PROPERLY FORMATTED
    vkey, versions = skempi_scraper()  # should be called out side of function
    # os.system('ls')
    raw_path = 'data/raw/'  # '../../data/raw/'    # make file path output at the bottom and make input to this function!
    print(versions)

    # pdb's must be downloaded first as versions only change csv.
    # download may not exist otherwise.
    if not os.path.isfile(f'{raw_path}skempi_{vkey}.tgz'):
        pdb_path, _ = urlretrieve(versions[vkey][1],
                                  f'{raw_path}skempi_{vkey}.tgz',
                                  reporthook=progress_hook)
        print(pdb_path)
    else:
        pdb_path = f'{raw_path}skempi_{vkey}.tgz'

    # extract PDBs
    if not os.path.isdir(f'{raw_path}PDBs'):
        os.system(f'tar -xvzf {pdb_path} --directory {raw_path}')
        os.remove(pdb_path)
        print('extracted')

    # Input version
    if v:
        vkey = v

    # Download requested version of csv.
    # Check if file exists
    csv_path, _ = urlretrieve(versions[vkey][0],
                              f'{raw_path}skempi_{vkey}.csv',
                              reporthook=progress_hook)
    print(csv_path)


# rename file to fetch_dataset?
# also test out Git-lfs
# cron or crontab – https://desktop.arcgis.com/en/arcmap/10.4/analyze/executing-tools/scheduling-a-python-script-to-run-at-prescribed-times.htm

logging.basicConfig(level=logging.INFO)
@click.command()
# @click.argument('input_filepath', type=click.Path(exists=True))
# @click.argument('output_filepath', type=click.Path())
def main():  # main(input_filepath, output_filepath):
    """ Runs data processing scripts to turn raw data from (../raw) into
        cleaned data ready to be analyzed (saved in ../processed).
    """
    skempi_download()

    logger = logging.getLogger(__name__)
    logger.info('making raw data set')


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main()
