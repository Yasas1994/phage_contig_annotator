import glob
import os
import sys
import argparse
from utils import * #call_genes, search_hmms, search_hmms_hhsuite, parse_hmmsearch, generate_plots_and_gff, run_combine
from pathlib import Path
import subprocess
import configparser


def download_dbs(path, extract_path):
    import requests
    import zipfile
    url =  "https://nextcloud.uni-greifswald.de/index.php/s/w2pgjQXdifsCtGA/download/databases.zip"
    # Ensure the download directory exists
    os.makedirs(path, exist_ok=True)

    # Download the zip file
    response = requests.get(url)
    if response.status_code == 200:
        zip_file_path = os.path.join(path, 'database.zip')
        with open(zip_file_path, 'wb') as file:
            file.write(response.content)

        # Ensure the extraction directory exists
        os.makedirs(extract_path, exist_ok=True)

        # Extract the zip file
        with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
            zip_ref.extractall(path)

        # Clean up the downloaded zip file
        os.remove(zip_file_path)

        logger.info(f"File downloaded and extracted to {extract_path}")
        return True
    else:
        logger.info(f"Failed to download the file from {url}. Status code: {response.status_code}")
        return False


libpath=os.path.dirname(os.path.realpath(__file__))
def main():
    parser = argparse.ArgumentParser(description='A pipeline to annotate genes in phage contigs with phrogs and VOGs\n')

    parser.add_argument("--cpus",
                        type = int,
                        required=False,
                        default=8,
                        help='set the number of cpus to use')

    parser.add_argument("--type",
                        required=False,
                        default='contigs',
                        help='define the input type default:contigs, options:[contigs,proteins]')

    subparsers = parser.add_subparsers(help='mode', dest='command')
    parser1 = subparsers.add_parser('runall', help='run the entire pipeline')
    parser1.add_argument("-i","--input",
                        type=is_valid_file_path,
                        required=True,
                        help="path to input fasta file with putative ")
    parser1.add_argument("--contigs",
                required=False,
                action='store_true',
                default=False,
                help='predict proteins')
    parser1.add_argument("-o","--output", 
                        type=is_valid_dir,
                        required=True,
                        help='path to output dir')
    parser1.add_argument("-db",
                        required=False,
                        type=is_valid_dir,
                        help='use a custom hmm database')
    #run custom 
    parser2 = subparsers.add_parser('custom', help='run a custom pipeline')
    parser2.add_argument('--run_prodigal',
                        required=False,
                        action='store_true',
                        default=False,
                        help='runs prodigal')
    parser2.add_argument('--run_hmmer',
                        required=False,
                        action='store_true',
                        default=False,
                        help='runs hmmer, input should contain valid protein sequences')
    parser2.add_argument('--run_combine',
                        required=False,
                        action='store_true',
                        default=False,
                        help='runs combine.py, merges all results to a single file')
    parser2.add_argument('--run_blast',
                        required=False,
                        action='store_true',
                        default=False,
                        help='runs blastp, input should contain valid protein sequences')
    parser2.add_argument('--run_diamond',
                        required=False,
                        action='store_true',
                        default=False,
                        help='runs diamond blastp,input should contain valid protein sequences')
    parser2.add_argument('--run_trnascan',
                        required=False,
                        action='store_true',
                        default=False,
                        help='runs trnascan, input should contain valid protein sequences')
    parser2.add_argument("-i","--input",
                        type=is_valid_file_path,
                        required=True,
                        help="path to input fasta file with putative ")
    parser2.add_argument("-o","--output", 
                        type=is_valid_dir,
                        required=True,
                        help='path to output dir')
    parser2.add_argument("-db",
                        required=False,
                        type=is_valid_dir,
                        help='use a custom hmm database')
    parser3 = subparsers.add_parser('download_db', help='run a custom pipeline')
    parser3.add_argument("-p","--path",
                        type=is_valid_dir,
                        required=False,
                        help="path to store the hmm database")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        
        parser.print_help()
        sys.exit(1)

    logger = get_logger(quiet=False)
    config_path = os.path.join(os.path.dirname(__file__), 'config.ini')
    config = configparser.ConfigParser()
    config.read(config_path)

    logger.info('phage_contig_annotator is a simple pipeline to add PHROG annotations to phage contigs')

    if args.command == 'download_db':
        logger.info('downloading PHROG database')
        if args.path == '':
            args.path=f'{libpath}/databases'
        if download_dbs(path=args.path):
            config.set('Databases','dbroot',args.path)
        sys.exit()

    if args.db:
        db_dir = args.db
    else:
        db_dir = f'{libpath}/databases'

    #search all subdirectories and determine how many dbs are avaiable
    hmmerdb = dict()
    blastdb = dict()
    diamonddb = dict()
    meta_path = str
    for path in glob.glob(f"{db_dir}/*"):

        if 'hmmerdb' in path:
            logger.info('hmmerdb found')
            db_name = path.split('/')[-1].split('_')[0]
            hmmerdb[db_name] = f"{path}"

        if 'blastdb' in path:
            logger.info('blastdb found')
            db_name = dbname(path)
            blastdb = f"{path}/{db_name}"

        if 'diamond' in path:
            logger.info('diamond_db found')
            db_name = dbname(path)
            diamonddb = f"{path}/{db_name}"

        if 'meta' in path:
            logger.info('metadata file found')
            meta_path = glob.glob(f"{path}/phrog_annot_v4.tsv")[0]

    #tmp dirs and file names
    fname = args.input.split('/')[-1].rsplit('.',1)[0]
    tmp_dir = f'{os.path.abspath(args.output)}/{fname}'

    fna_dir = f'{tmp_dir}/fna'
    prot_dir = f"{tmp_dir}/proteins" #this can be done inside the call genes funciton

    hmmsearch_dirs={}
    for key in hmmerdb.keys():
        hmmsearch_dirs[key] = f"{tmp_dir}/hmmsearch_{key}"

    plots_dir= f"{tmp_dir}/plots"
    combine_dir = f"{tmp_dir}/combined"

    #creating dir
    for i in hmmsearch_dirs.values():
        Path(i).mkdir(parents=True, exist_ok=True)

    Path(plots_dir).mkdir(parents=True, exist_ok=True)
    Path(combine_dir).mkdir(parents=True, exist_ok=True)

    if args.type != 'contigs':
        #create a symlobic link to the input file in tmp dir
        try:
            os.symlink(os.path.abspath(args.input), f"{tmp_dir}/proteins.faa")
        except:
            logger.warning('symlink exists')
        

    gff_path = f'{tmp_dir}/proteins.gff' #path to save .gff 
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)

    if args.command == 'runall':
        if args.contigs:
            
            call_genes(f'{args.input}', threads=args.cpus, tmp_dir=tmp_dir, trna=True)
        for key in hmmerdb.keys(): #run hmmsearch for all available databases 
            search_hmms(hmmout_dir=hmmsearch_dirs[key],proteins_dir=prot_dir, threads=args.cpus, db_dir=hmmerdb[key], tmp_dir=tmp_dir)
        generate_plots_and_gff(tmp_dir=tmp_dir,trna_dir=f'{tmp_dir}/trna.tsv',  hmmsearch_dir=f'{tmp_dir}/hmmsearch.txt', meta_dir=meta_path, gff_dir=gff_path,) 

    elif args.command == 'custom':
        if args.run_prodigal:
            call_genes(f'{args.input}', threads=args.cpus,tmp_dir=tmp_dir,trna=args.run_trnascan)
        if args.run_hmmer:
            for key in hmmerdb.keys():
                search_hmms(hmmout_dir=hmmsearch_dirs[key],proteins_dir=prot_dir, threads=args.cpus, db_dir=hmmerdb[key], tmp_dir=tmp_dir)

        if args.run_combine:
            code = run_combine(protein_gff=f"{tmp_dir}/proteins.gff", trna_tsv=f"{tmp_dir}/trna.tsv", hmmsearch=f"{tmp_dir}/hmmsearch.csv", annotation=meta_path,output=f"{tmp_dir}/{fname}.csv",out=combine_dir)
            if not code:
                logger.error('summarizing failed')

if __name__ == "__main__":
    main()