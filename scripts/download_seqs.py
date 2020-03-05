from urllib.request import urlretrieve
import pandas as pd
import os
import gzip
import shutil

print('downloading seqs!!!')
#infiles = snakemake.params['infiles']
def gunzip(infile):
    '''Unzip a file with .gz extension. Will remove extension in outfile'''
    with gzip.open(infile, 'rb') as f_in:
        with open(infile.rstrip('.gz'), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(infile)

ann_df = pd.read_csv(os.path.join(snakemake.config['parameter_dir'],
snakemake.config['seqs_and_annotations'])).set_index('organism', drop = False)

#ann_df = pd.read_csv(snakemake.params['ann_file']).set_index('organism', drop = False)
#can you pass a dict from params?
download_dict = snakemake.params['to_download']
outdir = snakemake.params['outdir']
for i in download_dict:
    address = ann_df.loc[snakemake.wildcards.org, i]
    gzipped_file = os.path.join(outdir, download_dict[i]) + '.gz'
    local_path, headers = urlretrieve(address, gzipped_file)
    gunzip(gzipped_file)
