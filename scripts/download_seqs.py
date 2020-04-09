from urllib.request import urlretrieve
import pandas as pd
import os
import gzip
import shutil

def gunzip(infile):
    '''
    Unzip a file with .gz extension. Will remove extension in outfile.
    If the file does not have a .gz extension, do not unzip.
    '''
    if not infile.endswith('.gz'):
        return
    with gzip.open(infile, 'rb') as f_in:
        with open(infile.rstrip('.gz'), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(infile)

def retrieve_file(address, outfile):
    '''
    If url, download. If file, move it to the new location.
    Note the path of a file should be the absolute path because
    snakemake will likely be run in --directory mode.
    '''
    try:
        local_path, headers = urlretrieve(address, outfile)
    #not url
    except ValueError:
        shutil.copy(address, outfile)

ann_df = pd.read_csv(os.path.join(snakemake.config['parameter_dir'],
snakemake.config['seqs_and_annotations'])).set_index('organism', drop = False)

download_dict = snakemake.params['to_download']
outdir = snakemake.params['outdir']
for i in download_dict:
    address = ann_df.loc[snakemake.wildcards.org, i]
    if address.endswith('.gz'):
        gzipped_file = os.path.join(outdir, download_dict[i]) + '.gz'
        retrieve_file(address, gzipped_file)
        gunzip(gzipped_file)
    else:
        retrieve_file(address, os.path.join(outdir, download_dict[i]))
