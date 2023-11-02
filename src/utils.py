import bz2
import gzip
import logging
import math
import multiprocessing as mp
import os
import platform
import resource
import signal
import subprocess as sp
import sys
import time
import glob
from collections import Counter
from enum import Enum, auto
from pathlib import Path
import argparse
import psutil
import pandas as pd
import numpy as np



from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from BCBio import GFF
from dna_features_viewer import BiopythonTranslator

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')

this = os.path.abspath(__file__)

#Most things are shamelessly copied from checkv
logger = logging.getLogger(__name__)

class Compression(Enum):
    gzip = auto()
    bzip2 = auto()
    xz = auto()
    noncompressed = auto()


def max_mem_usage():
    """Return max mem usage (GB) of self and child processes"""
    max_mem_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    max_mem_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if platform.system() == "Linux":
        return (max_mem_self + max_mem_child) / float(1e6)
    else:
        return (max_mem_self + max_mem_child) / float(1e9)


def is_compressed(filepath):
    """Checks if a file is compressed (gzip, bzip2 or xz)"""
    with open(filepath, "rb") as fin:
        signature = fin.peek(8)[:8]
        if tuple(signature[:2]) == (0x1F, 0x8B):
            return Compression.gzip
        elif tuple(signature[:3]) == (0x42, 0x5A, 0x68):
            return Compression.bzip2
        elif tuple(signature[:7]) == (0xFD, 0x37, 0x7A, 0x58, 0x5A, 0x00, 0x00):
            return Compression.xz
        else:
            return Compression.noncompressed


def get_compressed_file_handle(path):
    filepath_compression = is_compressed(path)
    if filepath_compression == Compression.gzip:
        f = gzip.open(path, "rt")
    elif filepath_compression == Compression.bzip2:
        f = bz2.open(path, "rt")
    elif filepath_compression == Compression.xz:
        f = lzma.open(path, "rt")
    else:
        f = open(path, "r")
    return f


def get_logger(quiet):
    logger = logging.getLogger(__name__)
    if not quiet:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
    formatter = logging.Formatter(fmt="%(asctime)s : %(levelname)s : %(message)s",datefmt="%d-%b-%y %H:%M:%S")
    stream_handler = logging.StreamHandler(sys.stderr)
    stream_handler.setFormatter(formatter)
    logger.handlers.clear()
    logger.addHandler(stream_handler)
    return logger

def dbname(path,**kwargs):
    '''gets the database name by parsing the database filename'''
    for file in glob(path+'/*'):
        file_name = os.path.basename(file)
        db_name = file_name.split('.')[0]
    return db_name

def is_valid_dir(path,**kwargs):
    '''checks path and creates if absent'''
    if os.path.isdir(path):
        
        return path
    elif os.path.isfile(path):
  
        sys.exit('a file with outdir name exists',1)
    else:
        os.mkdir(path)
        return path


def is_valid_file_path(path,**kwargs):
    '''checks if file is present'''
    if os.path.isfile(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"ERROR:{path} is not a valid file")
    
def check_fasta(path, tmp_dir,**kwargs):
    checkpoint_file = os.path.join(tmp_dir, "input_validation_checkpoint")
    if not os.path.isfile(checkpoint_file):
        f = get_compressed_file_handle(path)
        fasta_parser = SeqIO.parse(f, "fasta")
        if not any(fasta_parser):
            f.close()
            sys.stderr.write("You input FASTA file is empty or not properly formatted.")
            sys.exit()
        else:
            f = get_compressed_file_handle(path)
            fasta_parser = SeqIO.parse(f, "fasta")
            seq_id_counter = Counter([record.id for record in fasta_parser])
            f.close()
            repeated_seq_ids = [i for i, j in seq_id_counter.items() if j > 1]
            if repeated_seq_ids:
                sys.stderr.write(
                    f"Please remove duplicated sequence IDs from the input FASTA file: {', '.join(repeated_seq_ids)}"
                )
                sys.exit()
            else:
                with open(checkpoint_file, "w") as fout:
                    pass


def check_executables(requirements,**kwargs):
    fails = 0
    for program in requirements:
        found = False
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path.strip('"'), program)
            if os.path.isfile(exe_file) and os.access(exe_file, os.X_OK):
                found = True
                break
        if not found:
            msg = f"Error: required program '{program}' not executable or not found on $PATH\n"
            sys.stderr.write(msg)
            fails += 1
    if fails > 0:
        sys.exit()


def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def terminate_tree(pid, including_parent=True):
    parent = psutil.Process(pid)
    for child in parent.children(recursive=True):
        child.terminate()
    if including_parent:
        parent.terminate()


def async_parallel(function, argument_list, threads,**kwargs):
    """Based on: https://gist.github.com/admackin/003dd646e5fadee8b8d6"""
    # threads = len(argument_list) ## why is this being defined again here?
    pool = mp.Pool(threads, init_worker)
    try:
        results = []
        for arguments in argument_list:
            p = pool.apply_async(function, args=arguments)
            results.append(p)
        pool.close()
        while True:
            if all(r.ready() for r in results):
                return [r.get() for r in results]
            time.sleep(1)
    except KeyboardInterrupt:
        # when you want to kill everything, including this program
        # https://www.reddit.com/r/learnpython/comments/7vwyez/how_to_kill_child_processes_when_using/dtw3oh4/
        pid = os.getpid()
        terminate_tree(pid)


def check_database(dbdir,**kwargs):
    """check existence of database blastp, diamond and hmm files"""
    if dbdir is None:
        if "COILDB" not in os.environ:
            msg = "Error: database dir not specified\nUse -d or set CHECKVDB environmental variable"
            sys.exit(msg)
        else:
            dbdir = os.environ["PCADB"]
    dbdir = os.path.abspath(dbdir)
    if not os.path.exists(dbdir):
        msg = f"Error: database dir not found '{dbdir}'"
        sys.exit(msg)
    files = [
        "blastdb/phrogs.*",
        "diamonddb/phrogs.*",
        "hmmdb/phrogs_with_annot*"
    ]
    for f in files:
        path = os.path.join(dbdir, f)
        if not glob.glob(path):
            msg = f"Error: database file not found '{path}'"
            sys.exit(msg)
    return dbdir


def read_fasta(path,**kwargs):
    """Read fasta file and yield (header, sequence)"""
    filepath_compression = is_compressed(path)
    if filepath_compression == Compression.gzip:
        f = gzip.open(path, "rt")
    elif filepath_compression == Compression.bzip2:
        f = bz2.open(path, "rt")
    elif filepath_compression == Compression.xz:
        f = lzma.open(path, "rt")
    else:
        f = open(path, "r")
    for record in SeqIO.parse(f, "fasta"):
        name = record.description
        seq = str(record.seq).upper()
        if name != "" and seq != "":
            yield name, seq
    f.close()



def run_prodigal(out,in_,**kwargs):
    cmd = "prodigal-gv "
    cmd += " -m "
    cmd += "-p meta "
    cmd += f"-i {in_}.fna "
    cmd += f"-a {out}.faa "
    cmd += f"-f gff "
    cmd += f"-o {out}.gff "
    cmd += "1> /dev/null "
    cmd += f"2> {out}.log"
    with open(f"{out}.cmd", "w") as file:
        file.write(cmd + "\n")
    p = sp.Popen(cmd, shell=True)
    return_code = p.wait()
    return return_code == 0

def run_trnascan(out,in_,**kwargs):
    cmd = "tRNAscan-SE "
    cmd += " -G "
    cmd += "-p meta "
    cmd += f"-o {out}.tsv "
    cmd += f"-j {out}.gff "
    cmd += f"-i {in_}.fna "
    cmd += "1> /dev/null "
    cmd += f"2> {out}.log"
    with open(f"{out}.cmd", "w") as file:
        file.write(cmd + "\n")
    p = sp.Popen(cmd, shell=True)
    return_code = p.wait()
    return return_code == 0

def run_diamond(out, db, faa, tmp, threads,**kwargs):
    cmd = "diamond blastp "
    cmd += "--outfmt 6 "
    cmd += "--evalue 1e-5 "
    cmd += "--query-cover 50 "
    cmd += "--subject-cover 50 "
    cmd += "-k 10000 "
    cmd += f"--query {faa} "
    cmd += f"--db {db} "
    cmd += f"--threads {threads} "
    cmd += f"> {out} "
    cmd += f"2> {tmp}.log"
    with open(f"{tmp}.cmd", "w") as file:
        file.write(cmd + "\n")
    p = sp.Popen(cmd, shell=True)
    return_code = p.wait()
    return return_code == 0

def run_ffindex_build(out,faa_dir,**kwargs):
    '''faa_dir must contain individual .faa files'''
    cmd = "ffindex_build"
    cmd += "-s"
    cmd += f"{out}.ffdata "
    cmd += f"{out}.ffindex "
    cmd += f"{faa_dir} "
    cmd += f"2> {out}.log "
    with open(f"{out}.cmd", "w") as file:
        file.write(cmd + "\n")
    p = sp.Popen(cmd, shell=True)
    return_code = p.wait()
    return return_code == 0

def run_ffindex_from_fasta(out,faa, **kwargs):
    '''expect input to be a single multi faa file'''
    cmd = "ffindex_from_fasta "
    cmd += "-s "
    cmd += f"{out}.ffdata "
    cmd += f"{out}.ffindex "
    cmd += f"{faa} "
    cmd += f"2> {out}.ffindex.log "
    with open(f"{out}.ffindex.cmd", "w") as file:
        file.write(cmd + "\n")
    p = sp.Popen(cmd, shell=True)
    return_code = p.wait()
    return return_code == 0  

def run_hhblits_omp(out, db, faa,threads=0, evalue=0.001,**kwargs):
    #-i ./ffindex/phage_msa  -d ./tmp/phage -e 0.001 -cpu 12 -z 3 -Z 3 -b 0 -B 0 -v 1 -M 50 -o outfile
    cmd = "hhblits_omp "
    cmd += "-i "
    cmd += f"{faa} "
    cmd += "-d "
    cmd += f"{db} "
    cmd += f"-e {evalue} "
    cmd += f"-blasttab {out}.tbl "
    cmd += f"-cpu {threads} "
    cmd += f"-z 3 "
    cmd += f"-Z 3 "
    cmd += f"-b 0 "
    cmd += f"-B 0 "
    cmd += "-M 50 "
    cmd += f"2> {out}.log "
    with open(f"{out}.cmd", "w") as file:
        file.write(cmd + "\n")
    p = sp.Popen(cmd, shell=True)
    return_code = p.wait()
    return return_code == 0


def run_hmmsearch(out, db, faa, threads=2, evalue=10, **kwargs):
    cmd = "hmmsearch "
    cmd += "--noali "
    cmd += "-o /dev/null "
    cmd += f"-E {evalue} "
    cmd += f"--tblout {out} "
    cmd += f"--cpu {threads} "
    cmd += f"{db} "
    cmd += f"{faa} "
    cmd += f"2> {out}.log "
    with open(f"{out}.cmd", "w") as file:
        file.write(cmd + "\n")
    p = sp.Popen(cmd, shell=True)
    return_code = p.wait()
    return return_code == 0

def run_combine(protein_gff, trna_tsv, hmmsearch, annotation,output,out, **kwargs):
    '''expect input to be a single multi faa file'''
    cmd = f"python {os.path.dirname(os.path.abspath(__file__))}/combine.py "
    cmd += f"{protein_gff} "
    cmd += f"{trna_tsv} "
    cmd += f"{hmmsearch} "
    cmd += f"{annotation} "
    cmd += f"{output} "
    cmd += f"2> {out}/combine.log "
    with open(f"{out}/combine.cmd", "w") as file:
        file.write(cmd + "\n")
    p = sp.Popen(cmd, shell=True)
    return_code = p.wait()
    return return_code == 0  

def search_hmms(hmmout_dir,proteins_dir, threads, db_dir, tmp_dir, **kwargs):
    tmp = f"{tmp_dir}/hmmsearch.txt"
    # list splits to process
    checkpoint_file = os.path.join(tmp_dir, "hmmsearch_chkpt")
    if not os.path.isfile(checkpoint_file):
        splits = []
        for file in os.listdir(db_dir):
            split = file.split(".")[0]
            out = os.path.join(hmmout_dir, f"{split}.hmmout")
            # file doesn't exist; add to list for processing
            if not os.path.exists(out):
                splits.append(split)
            # check if file is complete
            else:
                x = False
                with open(out) as subf:
                    for line in subf:
                        if line == "# [ok]\n":
                            x = True
                if not x:
                    splits.append(split)
        # run hmmer
        logger.info('running hhmsearch')
        args_list = []
        for split in splits:
            out = os.path.join(hmmout_dir, f"{split}.hmmout")
            hmmdb = os.path.join(db_dir, f"{split}.hmm")
            faa = os.path.join(proteins_dir+".faa")
            args_list.append([out, hmmdb, faa])
        results = async_parallel(run_hmmsearch, args_list, threads)
        if not all(results):
            num_fails = len(results) - sum(results)
            sys.exit(
                f"\nError: {num_fails} hmmsearch tasks failed. Program should be rerun."
            )
        # check outputs are complete
        logger.info('checking the outputs are complete')
        complete = []
        for file in os.listdir(hmmout_dir):
            if file.split(".")[-1] == "hmmout":
                x = False
                with open(os.path.join(hmmout_dir, file)) as subf:
                    for line in subf:
                        if line == "# [ok]\n":
                            x = True
                complete.append(x)
        num_fails = complete.count(False)
        if num_fails > 0:
            sys.exit(
                f"\nError: {num_fails}/80 hmmsearch tasks failed. Program should be rerun."
            )
        # cat output
        logger.info('gathering search results')
        with open(tmp, "a") as f:
            for file in os.listdir(hmmout_dir):
                if file.split(".")[-1] == "hmmout":
                    with open(os.path.join(hmmout_dir, file)) as subf:
                        for line in subf:
                            f.write(line)
        search_results=pd.DataFrame(parse_hmmsearch(tmp))
        search_results.to_csv(f"{tmp_dir}/hmmsearch.csv",index=False)
        Path(checkpoint_file).touch()
    else:
        logger.info('hhmsearch checkpoint found')


def search_hmms_hhsuite(tmp_dir, threads, db_dir,**kwargs):
    # make tmp
    hmm_dir = os.path.join(tmp_dir, "hhsuite")
    index_dir = os.path.join(hmm_dir,"index")
    hhblits_dir = os.path.join(hmm_dir,"hhblits")
    for dir in [hmm_dir,index_dir,hhblits_dir]:
        if not os.path.exists(dir):
            os.makedirs(dir)
    # list faa files
    all_proteins = os.path.join(tmp_dir, "proteins.faa")
    # build index
    logger.info('builing index for parallel execution of hhblits')
    index_file = os.path.join(index_dir,"index")
    result = run_ffindex_from_fasta(index_file,all_proteins)
    if not result:
        sys.exit(
            logger.error(f"\nError: building index for hhblits. Program should be rerun.")
        )
    # run hhblits
    logger.info('hmmdb search started')
    hhout = os.path.join(hhblits_dir, "hhblits")
    result = run_hhblits_omp(out=hhout, db=db_dir, faa=index_file,threads=threads, evalue=0.001)
    if not result:
        sys.exit(
            logger.error(f"\nError: hhbits failed to run. Program should be rerun.")
        )
    logger.info('hmmdb search ended')
    out = os.path.join(tmp_dir,'hhblits.tsv')
    #extracting search results
    logger.info('unpacking search results')
    with open(hhout+'.tbl.ffdata','r') as fh, open(out,'w') as wh:
        wh.write('query\ttarget\t#match/tLen\talnLen\t#mismatch\t#gapOpen\tqstart\tqend\ttstart\ttend\teval\tscore\n')
        for line in fh:
            wh.write(line.strip('\x00'))
    

def call_genes(in_fna, threads, tmp_dir, trna=True, **kwargs):
    #check is checkpoint is passed

    checkpoint_file = os.path.join(tmp_dir, "gene_calling_chkpt")
    if not os.path.isfile(checkpoint_file):
        logger.info('gene calling started')
        fna_dir = os.path.join(tmp_dir,"fna")
        proteins_dir = os.path.join(tmp_dir,"proteins")
        trna_dir = os.path.join(tmp_dir,"trna")
        Path(fna_dir).mkdir(parents=True, exist_ok=True)
        Path(proteins_dir).mkdir(parents=True, exist_ok=True)
        Path(trna_dir).mkdir(parents=True, exist_ok=True)
        # count seqs in fasta
        num_seqs = sum(1 for _ in read_fasta(in_fna))
        # split fna into equal sized chunks
        split_size = int(math.ceil(1.0 * num_seqs / threads))
        iteration = 1
        count = 0
        out = open(os.path.join(fna_dir, f"{iteration}.fna"), "w")
        for id, seq in read_fasta(in_fna):
            # check if new file should be opened
            if count == split_size:
                count = 0
                iteration += 1
                out = open(os.path.join(fna_dir, f"{iteration}.fna"), "w")
            # write seq to file
            out.write(">" + id + "\n" + seq + "\n")
            count += 1
        out.close()
        # call genes
        args_list = []
        for i in range(1, iteration + 1):
            in_ = os.path.join(fna_dir, str(i))
            out = os.path.join(proteins_dir, str(i))
            args_list.append([out, in_])
        results = async_parallel(run_prodigal, args_list, threads)
        if not all(results):
            num_fails = len(results) - sum(results)
            sys.exit(
                logger.error(f"\nError: {num_fails} prodigal tasks failed. Program should be rerun.")
            )
        # cat output faa
        # mapping = dict()
        with open(f"{proteins_dir}.faa", "w") as f:
            for i in range(1, iteration + 1):
                # avoid trying to read empty fasta file
                if i <= threads:
                    with open(os.path.join(proteins_dir, f"{i}.faa")) as subf:
                        j = 0
                        for line in subf:
                            #if line[0] == '>':
                                #j += 1
                                #linex = line.split('cov')[0] + f'{j}\n'
                                #mapping[line] = linex
                            f.write(line)
        # with open(f'{tmp}.pkl', 'wb') as f:
        #     pickle.dump(mapping, f)


        if trna:
            logger.info('calling trna genes')
            args_list = []
            for i in range(1, iteration + 1):
                out = os.path.join(trna_dir, str(i))
                in_ = os.path.join(fna_dir, str(i))
        
                args_list.append([out,in_])
            results = async_parallel(run_trnascan, args_list, threads)
            if not all(results):
                num_fails = len(results) - sum(results)
                sys.exit(
                    logger.error(f"\nError: {num_fails} tRNAscan-SE tasks failed. Program should be rerun.")
                )
    #cat output trnascan
            with open(f"{trna_dir}.tsv", "w") as f:
                for i in range(1, iteration + 1):
                    # avoid trying to read empty fasta file
                    if i <= threads:
                        with open(os.path.join(trna_dir, f"{i}.tsv")) as subf:
                            j = 0
                            for line in subf:
                                #if line[0] == '>':
                                    #j += 1
                                    #linex = line.split('cov')[0] + f'{j}\n'
                                    #mapping[line] = linex
                                f.write(line)
            with open(f"{trna_dir}.gff", "w") as f:
                    f.write('##gff-version  3\n')
                    for i in range(1, iteration + 1):
                        # avoid trying to read empty fasta file
                        if i <= threads:
                            with open(os.path.join(trna_dir, f"{i}.gff")) as subf:
                                j = 0
                                for index,line in enumerate(subf):
                                    if index > 0:
                                        f.write(line)
        #cat output gff
        with open(f"{proteins_dir}.gff", "w") as f:
            f.write('##gff-version  3\n')
            for i in range(1, iteration + 1):
                # avoid trying to read empty fasta file
                if i <= threads:
                    with open(os.path.join(proteins_dir, f"{i}.gff")) as subf:
                        j = 0
                        for index,line in enumerate(subf):
                            if index > 0:
                                f.write(line)
        Path(checkpoint_file).touch() #create the chkpt file
    else:
        logger.info('gene calling checkpoint found')

def parse_blastp(path, **kwargs):
    with open(path) as f:
        names = [
            "qname",
            "tname",
            "pid",
            "aln",
            "mis",
            "gap",
            "qstart",
            "qstop",
            "tstart",
            "tstop",
            "eval",
            "score",
        ]
        formats = [str, str, float, int, int, int, int, int, int, int, float, float]
        for line in f:
            values = line.split()
            yield dict([(names[i], formats[i](values[i])) for i in range(12)])


def parse_hmmsearch(path,**kwargs):
    with open(path) as f:
        names = [
            "qname",
            "qacc",
            "tname",
            "tacc",
            "eval",
            "score",
            "bias",
            "beval",
            "bscore",
            "bbias",
        ]
        formats = [str, str, str, str, float, float, float, float, float, float]
        for line in f:
            if not line.startswith("#"):
                values = line.split()
                try:
                    yield dict([(names[i], formats[i](values[i])) for i in range(10)])
                except:
                    print("skipping erroneous line")

#parse tRNAscan-SE output file
def parse_trna(path, **kwargs):
    with open(path) as f:
        names = [
            "qname",
            "trna_no",
            "begin",
            "end",
            "trna_type",
            "anticodon",
            "intron_begin",
            "intron_end",
            "score",
        ]
        formats = [str, int, int, int, str, str, int, int, float]
        for line in f:
            if not (line.startswith("Sequence") or line.startswith("Name") or line.startswith("-----")):
                values = line.split()
                try:
                    yield dict([(names[i], formats[i](values[i])) for i in range(9)])
                except Exception as e:
                    print(f"{e}/;skipping erroneous line")

def get_cordinates(x, **kwargs):
    if x['begin'] > x['end']:
        return pd.Series([x["qname"],x["trna_no"],x['end'],x['begin'],-1,x["trna_type"],x["score"]], index=["contig","trna_no","begin","end","strand","trna_type","score"])
    else:
        return pd.Series([x["qname"],x["trna_no"],x['begin'],x['end'],1,x["trna_type"],x["score"]], index=["contig","trna_no","begin","end","strand","trna_type","score"])

def create_feature(x, **kwargs):
    qualifiers = {
        "source": "tRNAscan-SE",
        "score": x["score"],
        "trna_type": x["trna_type"],
        "score" : x["score"],
        "label" : x["trna_type"]+"_tRNA",
        "ID": "trna_"+str(x["trna_no"]),
    }

    return (SeqFeature(FeatureLocation(x["begin"], x["end"]), type="tRNA",id=str(x["trna_no"]), strand=x["strand"], qualifiers=qualifiers))


def generate_plots_and_gff(tmp_dir, hmmsearch_dir, trna_dir ,meta_dir, gff_dir, **kwargs):    
    #check if tmp/plots exists, eles create the dir

    checkpoint_file = os.path.join(tmp_dir, "plotting_chkpt")
    if not os.path.isfile(checkpoint_file):

        plots_dir = os.path.join(tmp_dir, "plots")
        Path(plots_dir).mkdir(parents=True, exist_ok=True)

        #process hmmsearch results
        logger.info('processing hmm results')
        search_results=pd.DataFrame(parse_hmmsearch(hmmsearch_dir))
        trna=pd.DataFrame(parse_trna(trna_dir))
        if search_results.empty:
            sys.stderr.write('Exiting because hmmsearch returned zero matches!')
            sys.exit( )
        
        if not trna.empty:
            trna=trna.apply(lambda x : get_cordinates(x) , axis=1).sort_values(by=["contig","begin"])
        else:
            trna = None

        phrogs_anno=pd.read_table(meta_dir)
        phrogs_anno=phrogs_anno.fillna('unknown function')
        search_results=pd.DataFrame(parse_hmmsearch(hmmsearch_dir))
        phrogs_anno['phrog']=phrogs_anno['phrog'].apply(lambda x : f'phrog_{x}')
        results_with_annotate = search_results.merge(phrogs_anno, how='inner', left_on='tname', right_on='phrog')
        results_with_annotate['position'] = results_with_annotate['qname'].apply(lambda x: int(x.split('_')[-1]))
        results_with_annotate['contig'] = results_with_annotate['qname'].apply(lambda x: x.rsplit('_',1)[0])
        
        results_filtered = results_with_annotate.iloc[results_with_annotate.groupby('qname')['score'].idxmax()].query('score > 50')

        logger.info('generating annotation plots')
        gff_out_dir = os.path.join(tmp_dir, "annotations.gff")
        with open(gff_out_dir, "w") as out_handle:
        
            for i in GFF.parse(gff_dir):
                
                
                tmp = results_filtered.query(f"contig == '{i.id}'")

                for pos,feature in enumerate(i.features, start=1):
                    tmp_feature = tmp.query(f"position == {pos}")[["category","color","annot","phrog","score","eval"]]
                    if not tmp_feature.empty:
                        #print(tmp_feature["annot"])
                        feature.qualifiers.update({"label":tmp_feature["annot"].values[0]})
                        feature.qualifiers.update({"category":tmp_feature["category"].values[0]})
                        feature.qualifiers.update({"eval":tmp_feature["eval"].values[0]})
                        feature.qualifiers.update({"score":tmp_feature["score"].values[0]})
                        feature.qualifiers.update({"color":tmp_feature["color"].values[0]})
                        feature.qualifiers.update({"phrog":tmp_feature["phrog"].values[0]})
                    else:
                        feature.qualifiers.update({"label":"unknown function"})
                        feature.qualifiers.update({"color": "#c9c9c9"})
                    
                #create trna features
                if trna is not None:
                    tmp_trna = trna.query(f"contig == '{i.id}'")
                    tmp_trna=tmp_trna.reset_index(drop=True)
                    if not tmp_trna.empty:
                        tmp_trna["feature"]=tmp_trna.apply(lambda x : create_feature(x), axis=1)
                        i.features.extend(tmp_trna["feature"].to_list())
                #write updated to gff record to a file
                graphic_record = BiopythonTranslator().translate_record(i) 
                GFF.write([i], out_handle) #write gff records to the out file
                for feat in graphic_record.features: #turns off the labels of cds' without annotations
                    if feat.label == "unknown function":
                        feat.label =None
                fig, ax1 = plt.subplots(1, 1, figsize=(15, 4))
                fig.tight_layout(pad=2.5) 
                ax, _ = graphic_record.plot(ax=ax1, strand_in_label_threshold=7,annotate_inline=False,figure_height=3 )
                ax.set_title(i.id)
                out_name = os.path.join(plots_dir, f"{i.id}.png")
                ax.figure.savefig(out_name, bbox_inches='tight')
                plt.clf()
                plt.close("all")

        Path(checkpoint_file).touch()
    else:
        logger.info('plotting checkpoint found')

