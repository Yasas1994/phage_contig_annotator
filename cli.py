import glob
import os
import sys
import argparse
from utils import call_genes, search_hmms, search_hmms_hhsuite, parse_hmmsearch, generate_plots

libpath=os.path.dirname(os.path.realpath(__file__))
def dir_path(path):
    '''checks path and creates if absent'''
    if os.path.isdir(path):
        return path
    else:
        os.mkdir(path)
        return path


def file_path(path):
    '''checks if file is present'''
    if os.path.isfile(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"ERROR:{path} is not a valid file")

parser = argparse.ArgumentParser(description='Annotate and visualize phage contigs\n')

parser.add_argument("-i","--input",
                    type=file_path,
                    required=True,
                    help="path to input fasta file with putative ")
parser.add_argument("-o","--output", 
                    type=dir_path,
                    required=True,
                    help='path to output dir')

args = parser.parse_args()


tmp_dir = f'{args.output}/tmp'
hmmsearch_dir = f'{args.output}/tmp/hmmsearch.txt'
meta_dir = f'{libpath}/databases/meta/phrog_annot_v3.tsv'
gff_dir = f'{args.output}/tmp/proteins.gff'
trna_dir = f'{args.output}/tmp/proteins_trnascan.tsv'
call_genes(f'{args.input}', f'{args.output}', 12)
search_hmms(tmp_dir, threads=12, db_dir=f'{libpath}/databases/hmmerdb/hmm_db_with_annot')
#search_hmms_hhsuite(tmp_dir, threads=2, db_dir='../databases/hhdb/phrogs_with_annot')
generate_plots(tmp_dir, hmmsearch_dir,trna_dir, meta_dir, gff_dir)