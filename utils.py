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
import pickle
import numpy as np
import psutil
from Bio import SeqIO
import pandas as pd
import numpy as np
from BCBio import GFF
import pickle
from dna_features_viewer import BiopythonTranslator
import matplotlib.pyplot as plt
import os


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO,datefmt="%d-%b-%y %H:%M:%S")

#Most things are shamelessly copied from checkv


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


def get_logger_old(quiet):
    if not quiet:
        logging.basicConfig(level=logging.INFO, format="%(message)s")
    else:
        logging.basicConfig(level=logging.WARNING, format="%(message)s")
    return logging.getLogger()


def get_logger(quiet):
    logger = logging.getLogger(__name__)
    if not quiet:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
    formatter = logging.Formatter(fmt="%(message)s")
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.handlers.clear()
    logger.addHandler(stream_handler)
    return logger


def check_fasta(path, tmp_dir):
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


def check_executables(requirements):
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


def async_parallel(function, argument_list, threads):
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


def check_database(dbdir):
    """check existence of database blastp, diamond and hmm files"""
    if dbdir is None:
        if "COILDB" not in os.environ:
            msg = "Error: database dir not specified\nUse -d or set CHECKVDB environmental variable"
            sys.exit(msg)
        else:
            dbdir = os.environ["CHECKVDB"]
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


def read_fasta(path):
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



def run_prodigal(out):
    cmd = "prodigal "
    cmd += " -m "
    cmd += "-p meta "
    cmd += f"-i {out}.fna "
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


def run_diamond(out, db, faa, tmp, threads):
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

def run_ffindex_build(out,faa_dir):
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

def run_ffindex_from_fasta(out,faa):
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

def run_hhblits_omp(out, db, faa,threads=0, evalue=0.001):
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


def run_hmmsearch(out, db, faa, threads=0, evalue=10):
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

def search_hmms(tmp_dir, threads, db_dir):
    # make tmp
    hmm_dir = os.path.join(tmp_dir, "hmmsearch")
    if not os.path.exists(hmm_dir):
        os.makedirs(hmm_dir)
    # list faa files
    faa = [
        file
        for file in os.listdir(os.path.join(tmp_dir, "proteins"))
        if file.split(".")[-1] == "faa"
    ]
    # list splits to process
    splits = []
    for file in os.listdir(db_dir):
        split = file.split(".")[0]
        out = os.path.join(hmm_dir, f"{split}.hmmout")
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
    logging.info('running hhmsearch')
    args_list = []
    for split in splits:
        out = os.path.join(hmm_dir, f"{split}.hmmout")
        hmmdb = os.path.join(db_dir, f"{split}.hmm")
        faa = os.path.join(tmp_dir, "proteins.faa")
        args_list.append([out, hmmdb, faa])
    results = async_parallel(run_hmmsearch, args_list, threads)
    if not all(results):
        num_fails = len(results) - sum(results)
        sys.exit(
            f"\nError: {num_fails} hmmsearch tasks failed. Program should be rerun."
        )
    # check outputs are complete
    logging.info('checking the outputs are complete')
    complete = []
    for file in os.listdir(hmm_dir):
        if file.split(".")[-1] == "hmmout":
            x = False
            with open(os.path.join(hmm_dir, file)) as subf:
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
    logging.info('gathering search results')
    with open(os.path.join(tmp_dir, "hmmsearch.txt"), "w") as f:
        for file in os.listdir(hmm_dir):
            if file.split(".")[-1] == "hmmout":
                with open(os.path.join(hmm_dir, file)) as subf:
                    for line in subf:
                        f.write(line)


def search_hmms_hhsuite(tmp_dir, threads, db_dir):
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
    logging.info('builing index for parallel execution of hhblits')
    index_file = os.path.join(index_dir,"index")
    result = run_ffindex_from_fasta(index_file,all_proteins)
    if not result:
        sys.exit(
            logging.error(f"\nError: building index for hhblits. Program should be rerun.")
        )
    # run hhblits
    logging.info('hmmdb search started')
    hhout = os.path.join(hhblits_dir, "hhblits")
    result = run_hhblits_omp(out=hhout, db=db_dir, faa=index_file,threads=12, evalue=0.001)
    if not result:
        sys.exit(
            logging.error(f"\nError: hhbits failed to run. Program should be rerun.")
        )
    logging.info('hmmdb search ended')
    out = os.path.join(tmp_dir,'hhblits.tsv')
    #extracting search results
    logging.info('unpacking search results')
    with open(hhout+'.tbl.ffdata','r') as fh, open(out,'w') as wh:
        wh.write('query\ttarget\t#match/tLen\talnLen\t#mismatch\t#gapOpen\tqstart\tqend\ttstart\ttend\teval\tscore\n')
        for line in fh:
            wh.write(line.strip('\x00'))
    

def call_genes(in_fna, out_dir, threads):
    # make tmp dir
    logging.info('gene calling started')
    tmp = f"{out_dir}/tmp/proteins"
    if not os.path.exists(tmp):
        os.makedirs(tmp)
    # count seqs in fasta
    num_seqs = sum(1 for _ in read_fasta(in_fna))
    # split fna into equal sized chunks
    split_size = int(math.ceil(1.0 * num_seqs / threads))
    iteration = 1
    count = 0
    out = open(os.path.join(tmp, f"{iteration}.fna"), "w")
    for id, seq in read_fasta(in_fna):
        # check if new file should be opened
        if count == split_size:
            count = 0
            iteration += 1
            out = open(os.path.join(tmp, f"{iteration}.fna"), "w")
        # write seq to file
        out.write(">" + id + "\n" + seq + "\n")
        count += 1
    out.close()
    # call genes
    args_list = []
    for i in range(1, iteration + 1):
        out = os.path.join(tmp, str(i))
        args_list.append([out])
    results = async_parallel(run_prodigal, args_list, threads)
    if not all(results):
        num_fails = len(results) - sum(results)
        sys.exit(
            logging.error(f"\nError: {num_fails} prodigal tasks failed. Program should be rerun.")
        )
    # cat output faa
    # mapping = dict()
    with open(f"{tmp}.faa", "w") as f:
        for i in range(1, iteration + 1):
            # avoid trying to read empty fasta file
            if i <= threads:
                with open(os.path.join(tmp, f"{i}.faa")) as subf:
                    j = 0
                    for line in subf:
                        #if line[0] == '>':
                            #j += 1
                            #linex = line.split('cov')[0] + f'{j}\n'
                            #mapping[line] = linex
                        f.write(line)
    # with open(f'{tmp}.pkl', 'wb') as f:
    #     pickle.dump(mapping, f)

    #cat output gff
    with open(f"{tmp}.gff", "w") as f:
        f.write('##gff-version  3\n')
        for i in range(1, iteration + 1):
            # avoid trying to read empty fasta file
            if i <= threads:
                with open(os.path.join(tmp, f"{i}.gff")) as subf:
                    j = 0
                    for index,line in enumerate(subf):
                        if index > 0:
                            f.write(line)


def parse_blastp(path):
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


def parse_hmmsearch(path):
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
                yield dict([(names[i], formats[i](values[i])) for i in range(10)])

def generate_plots(tmp_dir, hmmsearch_dir, meta_dir,gff_dir):    
    #check if tmp/plots exists, eles create the dir
    
    plots_dir = os.path.join(tmp_dir, "plots")
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)
    #process hmmsearch results
    logging.info('processing hmm results')
    phrogs_anno=pd.read_table(meta_dir)
    phrogs_anno=phrogs_anno.fillna('unknown function')
    search_results=pd.DataFrame(parse_hmmsearch(hmmsearch_dir))
    phrogs_anno['phrog']=phrogs_anno['phrog'].apply(lambda x : f'phrog_{x}')
    results_with_annotate = search_results.merge(phrogs_anno, how='inner', left_on='tname', right_on='phrog')
    results_with_annotate['position'] = results_with_annotate['qname'].apply(lambda x: int(x.split('_')[-1]))
    results_with_annotate['contig'] = results_with_annotate['qname'].apply(lambda x: x.rsplit('_',1)[0])

    results_filtered = results_with_annotate.iloc[results_with_annotate.groupby('qname')['score'].idxmax()]


    logging.info('generating annotation plots')
    gff_out_dir = os.path.join(tmp_dir, "proteins_annot.gff")
    with open(gff_out_dir, "w") as out_handle:
    
        for i in GFF.parse(gff_dir):
            
            
            tmp = results_filtered.query(f"contig == '{i.id}'")
            for pos,feature in enumerate(i.features):
                
                tmp_feature = tmp.query(f"position == {pos}")[["category","color","annot"]]
                if not tmp_feature.empty:
                    #print(tmp_feature["annot"])
                    feature.qualifiers.update({"label":tmp_feature["annot"].values[0]})
                    feature.qualifiers.update({"color":tmp_feature["color"].values[0]})
                else:
                    feature.qualifiers.update({"label":"unknown function"})
                    feature.qualifiers.update({"color": '#c9c9c9'})
            graphic_record = BiopythonTranslator().translate_record(i)
            GFF.write([i], out_handle)
            for feat in graphic_record.features: #turns off the labels of cds' without annotations
                if feat.label == "unknown function":
                    feat.label =None
            fig, ax1 = plt.subplots(1, 1, figsize=(15, 5))
            fig.tight_layout(pad=2.5) 
            ax, _ = graphic_record.plot(ax=ax1, strand_in_label_threshold=7,annotate_inline=False,figure_height=3 )
            ax.set_title(i.id)
            out_name = os.path.join(plots_dir, f"{i.id}.png")
            ax.figure.savefig(out_name, bbox_inches='tight')
            plt.close()

