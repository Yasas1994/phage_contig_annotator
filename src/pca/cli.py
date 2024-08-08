import glob
import os
import sys
import io
import argparse
from pathlib import Path
import configparser
import logging
import requests
import tarfile
import tqdm
from importlib.resources import files
from pca.utils import (get_logger,
                       is_valid_file_path,
                       is_valid_dir,
                       dbname,
                       call_genes,
                       search_hmms,
                       search_hmms_hhsuite,
                       parse_hmmsearch,
                       generate_plots_and_gff,
                       run_combine)

# To-do
# add support for custom databases [.hmms & metadata]


class TqdmToLogger(io.StringIO):
    """
        Output stream for TQDM which will output to logger module instead of
        the StdOut.
    """
    logger = None
    level = None
    buf = ''

    def __init__(self, logger, level=None):
        super(TqdmToLogger, self).__init__()
        self.logger = logger
        self.level = level or logging.INFO

    def write(self, buf):
        self.buf = buf.strip('\r\n\t ')

    def flush(self):
        self.logger.log(self.level, self.buf)


def download_dbs(path: str, logger: str, force: bool = False) -> bool:
    if (not Path(path).joinpath('db_chkpt').exists()) or force:
        url = "https://nextcloud.uni-greifswald.de/index.php/s/ft8FAoQXscoj9eo/download/database.tar.gz"
        logger.info(url)
        # Ensure the download directory exists
        os.makedirs(path, exist_ok=True)
        # 1KB = 1024 bytes
        tar_file_path = os.path.join(path, 'database.tar.gz')
        with requests.get(url, stream=True) as response:
            response.raise_for_status()

            pbar = tqdm.tqdm(unit='B', unit_scale=True, unit_divisor=1000)

            with open(tar_file_path, "wb") as file:
                for chunk in response.iter_content(chunk_size=1024):
                    if chunk:
                        pbar.update(1024)
                        file.write(chunk)

            # Ensure the extraction directory exists
            os.makedirs(path, exist_ok=True)

            # Extract the tar.gz file
            with tarfile.open(tar_file_path, "r:gz") as tar:
                tar.extractall(path=path)
                # Clean up the downloaded zip file
            os.remove(tar_file_path)

            logger.info(f"database downloaded and extracted to {path}")
            Path(path).joinpath('db_chkpt').touch()
            return True

    else:
        logger.info(f"Skipping database download as checkpoint found at {path}")
        return True


def main():

    libpath = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser(description='A pipeline to annotate genes in phage contigs with phrogs and VOGs\n')

    subparsers = parser.add_subparsers(help='mode', dest='command')
    parser1 = subparsers.add_parser('runall', help='run the entire pipeline')

    parser1.add_argument("-i",
                         "--input",
                         type=is_valid_file_path,
                         required=True,
                         help="path to input fasta file with putative ")
    parser1.add_argument("--contigs",
                         required=False,
                         action='store_true',
                         default=False,
                         help='predict proteins')
    parser1.add_argument("-o",
                         "--output",
                         type=is_valid_dir,
                         required=True,
                         help='path to output dir')
    parser1.add_argument("-db",
                         required=False,
                         type=is_valid_dir,
                         help='use a custom hmm database')
    parser1.add_argument("--cpus",
                         type=int,
                         required=False,
                         default=8,
                         help='set the number of cpus to use')

    parser1.add_argument("--type",
                         required=False,
                         default='contigs',
                         help='define the input type default:contigs, options:[contigs, proteins]')
    # run custom
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
    parser2.add_argument("-i",
                         "--input",
                         type=is_valid_file_path,
                         required=True,
                         help="path to input fasta file with putative ")
    parser2.add_argument("-o",
                         "--output",
                         type=is_valid_dir,
                         required=True,
                         help='path to output dir')
    parser2.add_argument("-db",
                         required=False,
                         type=is_valid_dir,
                         help='use a custom hmm database')
    parser2.add_argument("--cpus",
                         type=int,
                         required=False,
                         default=8,
                         help='set the number of cpus to use')

    parser2.add_argument("--type",
                         required=False,
                         default='contigs',
                         help='define the input type default:contigs, options:[contigs, proteins]')

    parser3 = subparsers.add_parser('download_db', help='downloads the database')
    parser3.add_argument("-p",
                         "--path",
                         type=is_valid_dir,
                         required=False,
                         help="path to store the hmm database")
    parser3.add_argument("-f",
                         required=False,
                         help="override the checkpoint and re-download the database")
    args = parser.parse_args()

    if len(sys.argv) == 1:

        parser.print_help()
        sys.exit(1)

    logger = get_logger(quiet=False)

    config_path = files('pca.data').joinpath('config.ini')
    config = configparser.ConfigParser()
    config.read(config_path)

    logger.info('phage_contig_annotator is a simple pipeline to add PHROG annotations to phage contigs')

    if args.command == 'download_db':
        logger.info('downloading PHROG database')
        if args.path is None:
            args.path = os.path.join(libpath, 'databases')
            os.makedirs(args.path, exist_ok=True)

        if download_dbs(path=args.path, logger=logger):
            config.set('databases',
                       'dbroot',
                       os.path.realpath(args.path))
            config.set('databases',
                       'hmmdb',
                       os.path.join(os.path.realpath(args.path), 'hmmdb'))
            config.set('databases',
                       'meta',
                       os.path.join(os.path.realpath(args.path), 'meta'))
            with open(config_path, 'w') as configfile:
                config.write(configfile)
        sys.exit(0)

    if args.db:
        db_dir = args.db
    else:
        db_dir = config['databases']['dbroot']
        logger.info(db_dir)

    # search all subdirectories and determine how many dbs are avaiable
    hmmerdb = dict()
    blastdb = dict()
    diamonddb = dict()
    hmmerdb_meta_path = str
    blastdb_meta_path = str
    diamonddb_meta_path = str

    for path in glob.glob(os.path.join(db_dir, "*")):

        if 'hmmerdb' in path:
            logger.info('hmmerdb found')
            db_name = os.path.basename(path).split('_')[0]
            hmmerdb[db_name] = path

        if 'blastdb' in path:
            logger.info('blastdb found')
            db_name = os.path.basename(path).split('_')[0]
            blastdb[db_name] = path

        if 'diamond' in path:
            logger.info('diamond_db found')
            db_name = os.path.basename(path).split('_')[0]
            diamonddb[db_name] = path

        if 'meta' in path:
            logger.info('metadata file found')
            hmmerdb_meta_path = glob.glob(os.path.join(path, "*"))[0]

    # tmp dirs and file names
    fname = os.path.basename(args.input).rsplit('.', 1)[0]
    tmp_dir = os.path.join(os.path.abspath(args.output), fname)

    fna_dir = os.path.join(tmp_dir, 'fna')
    # this can be done inside the call genes funciton
    prot_dir = os.path.join(tmp_dir, 'proteins')
    plots_dir = os.path.join(tmp_dir, 'plots')
    combine_dir = os.path.join(tmp_dir, 'combined')

    hmmsearch_dirs = {}
    for key in hmmerdb.keys():
        hmmsearch_dirs[key] = os.path.join(tmp_dir, f"hmmsearch_{key}")

    # creating dirs
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)
    if args.type != 'contigs':
        # create a symlobic link to the input file in tmp dir
        try:
            os.symlink(os.path.abspath(args.input),
                       os.path.join(tmp_dir, 'proteins.faa'))
        except Exception as e:
            logger.warning(e)
    else:
        Path(fna_dir).mkdir(parents=True, exist_ok=True)
        Path(prot_dir).mkdir(parents=True, exist_ok=True)

    for i in hmmsearch_dirs.values():
        Path(i).mkdir(parents=True, exist_ok=True)

    Path(plots_dir).mkdir(parents=True, exist_ok=True)
    Path(combine_dir).mkdir(parents=True, exist_ok=True)

    if args.command == 'runall':
        if args.type == 'contigs':

            call_genes(os.path.realpath(args.input),
                       threads=args.cpus,
                       tmp_dir=tmp_dir,
                       trna=True)

        for key in hmmerdb.keys():  # run hmmsearch for all available databases

            search_hmms(hmmout_dir=hmmsearch_dirs[key],
                        proteins_dir=prot_dir,
                        threads=args.cpus,
                        db_dir=hmmerdb[key],
                        tmp_dir=tmp_dir)

        generate_plots_and_gff(tmp_dir=tmp_dir,
                               trna_dir=os.path.join(tmp_dir, 'trna.gff'),
                               hmmsearch_dir=os.path.join(tmp_dir,
                                                          'hmmsearch.txt'),
                               meta_dir=hmmerdb_meta_path,
                               gff_dir=os.path.join(tmp_dir,
                                                    'proteins.gff'))

    elif args.command == 'custom':
        if args.run_prodigal:
            call_genes(args.input,
                       threads=args.cpus,
                       tmp_dir=tmp_dir,
                       trna=args.run_trnascan)
        if args.run_hmmer:
            for key in hmmerdb.keys():
                search_hmms(hmmout_dir=hmmsearch_dirs[key],
                            proteins_dir=prot_dir,
                            threads=args.cpus,
                            db_dir=hmmerdb[key],
                            tmp_dir=tmp_dir)

        if args.run_combine:
            code = run_combine(protein_gff=os.path.join(tmp_dir,
                                                        'proteins.gff'),
                               trna_tsv=os.path.join(tmp_dir,
                                                     'trna.gff'),
                               hmmsearch=os.path.join(tmp_dir,
                                                      'hmmsearch.csv'),
                               annotation=hmmerdb_meta_path,
                               output=os.path.join(tmp_dir,
                                                   f"{fname}.csv"),
                               out=combine_dir)
            if not code:
                logger.error('summarizing failed')
                sys.exit(1)


if __name__ == "__main__":
    main()
