# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv


@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.argument('output_filepath', type=click.Path())
def main(input_filepath, output_filepath):
    """ Runs data processing scripts to turn raw data from (../raw) into
        cleaned data ready to be analyzed (saved in ../processed).
    """
    logger = logging.getLogger(__name__)
    logger.info('making final data set from raw data')


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main()



#rough draft for now...
#(use inheritence for more tailored use later) or (set these as functions and dont us a class)?
#probabaly overkill for what is likely to be static data anyway, oh well

#run skempi.check_latest on regular basis
#cron or crontab – https://desktop.arcgis.com/en/arcmap/10.4/analyze/executing-tools/scheduling-a-python-script-to-run-at-prescribed-times.htm
class skempi: 
    def __init__(self):
        self.version = None
        self.data = None
        self.pdb = 'path(s) and status?' #store raw pdbs in raw

    def __repr__(self):
        if self.data== None:
            return 'No data initialized'
        if self.version== None:
            return 'No data, hence no version'
    

    def check_version(self):
        #try
            #check if file exists get version from name
        #except
            #file missing...
            #download new
        pass
        
    def check_latest(self):
        #if up to date
            #pass
        #else
            #set path
            #download new
        pass
    
    def download(self):
        #run wget from py?
        #save.replace file as "skempi v#.#" to raw
        #data=extract pandas dataframe
        #wget pdb and send to raws 
        pass


#make more universal class and use inheritence... 
class zemu: 
    def __init__(self):
        self.version = None
        self.data = None  # make empty dataframe?
        self.pdb = 'path(s) and status?' #

    def check_existence(self):
        #try
            #data empty?
            #do pdbs and folder exist?
        #else:
            #call download
        pass
    
    def download(self):
        #run wget from py?
        #save.replace file as "skempi v#.#" to raw
        #data=extract pandas dataframe
        #wget pdb and send to raws 
        pass