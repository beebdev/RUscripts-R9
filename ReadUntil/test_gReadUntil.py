import time
import errno
from socket import error as socket_error
import threading
import sys, os, re
from Bio import SeqIO
from io import StringIO
import string
import mlpy
import sklearn.preprocessing
import random
import math
import csv
import numpy as np
import array as ar
import configargparse
import shutil
import pickle
import multiprocessing
import subprocess
import re
import glob
import h5py
import platform

from ruutils import process_model_file, check_files, check_fasta, process_ref_fasta, squiggle_search2, go_or_no

import coloredlogs, logging
logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)

def get_seq_len(ref_fasta):
    seqlens = dict()
    for record in SeqIO.parse(ref_fasta, 'fasta'):
        seq = record.seq
        seqlens[record.id] = len(seq)
    return seqlens

def process_hdf5(arg):
    # Extract info from arguments
    filename, seqids, threedarray, procampres, seqlen, args = arg

    # Open hdf5 file
    hdf = h5py.File(filename, 'r')

    # Process each read in the hdf file
    for read in hdf['Analyses']['EventDetection_000']['Reads']:
        events = hdf['Analyses']['EventDetection_000']['Reads'][read]['Events'][()]
        event_collection = list()
        for event in events:
            event_collection.append(float(event[0]))
        
        # We ignore the first 50 events (Protein driver) and process the following 250 events
        squiggle = event_collection[50:300]

        # Search squiggle in reference squiggle
        squiggleres = squiggle_search2(squiggle, 0, 0, args, seqids, threedarray, seqlen)

        if True:
            # logger.info("[%s, %s, %s]", squiggleres[0], squiggleres[2], squiggleres[3])
            try:
                result = go_or_no(squiggleres[0], squiggleres[2], squiggleres[3], seqlen, args)
            except Exception as err:
                logger.error("%s", err)
    hdf.close()
    return (result, filename, squiggleres)

def mycallback(arg):
    (result, filename, squiggleres) = arg
    filetocheck = os.path.split(filename)
    sourcefile = filename

    if result == "Sequence":
        path1 = os.path.join(args.output_folder,'sequence')
        path2 = os.path.join(path1,'downloads')
        path3 = os.path.join(path2,'pass')
        path4 = os.path.join(path2,'fail')

        if not os.path.exists(path1):
            os.makedirs(path1)
        if not os.path.exists(path2):
            os.makedirs(path2)
        if not os.path.exists(path3):
            os.makedirs(path3)
        if not os.path.exists(path4):
            os.makedirs(path4)

        logger.info("[%s-%s @%s] \033[42;1mSequence found\033[0m\n[%s]", squiggleres[0], squiggleres[2], squiggleres[3], filename)
        if "pass" in filename:
            destfile = os.path.join(path3,filetocheck[1])
        else:
            destfile = os.path.join(path4,filetocheck[1])
        try:
            shutil.copy(sourcefile,destfile)
        except Exception as err:
            logger.error("File Copy Failed: %s", err)
    else:
        path1 = os.path.join(args.output_folder,'reject')
        path2 = os.path.join(path1,'downloads')
        path3 = os.path.join(path2,'pass')
        path4 = os.path.join(path2,'fail')

        if not os.path.exists(path1):
            os.makedirs(path1)
        if not os.path.exists(path2):
            os.makedirs(path2)
        if not os.path.exists(path3):
            os.makedirs(path3)
        if not os.path.exists(path4):
            os.makedirs(path4)

        logger.info("[%s-%s @%s] \033[41;1mNo match\033[0m\n[%s]", squiggleres[0], squiggleres[2], squiggleres[3], filename)
        if "pass" in filename:
            destfile = os.path.join(path3,filetocheck[1])
        else:
            destfile = os.path.join(path4,filetocheck[1])
        try:
            #os.symlink(sourcefile, destfile)
            shutil.copy(sourcefile,destfile)
        except Exception as err:
            logger.error("File copy failed: %s", err)

if __name__ == "__main__":

    print ("***********************************************************************************************")
    print ("**** This code will open a collection of reads and simulate read until on them. It will    ****")
    print ("**** copy reads into a secondary folder for subsequent processing by another analysis      ****")
    print ("**** package.                                                                              ****")
    print ("***********************************************************************************************")

    global oper

    oper = platform.system()
    if oper == 'Windows':  # MS
        oper = 'windows'
    else:
        oper = 'linux'  # MS


    ## Linux version
    if (oper == "linux"):
            config_file = os.path.join(os.path.sep, os.path.dirname(os.path.realpath('__file__')), 'amp.config')

    ## Windows version
    if (oper == "windows"):
            config_file = os.path.join(os.path.sep, os.path.dirname(os.path.realpath('__file__')), 'ampW.config')

    __version__ = "1.1"
    __date__ = "1st May 2016"

    # Argument parsing
    parser = configargparse.ArgParser(description='real_read_until: A program providing read until with the Oxford Nanopore minION device. This program will ultimately be driven by minoTour to enable selective remote sequencing. This program is heavily based on original code generously provided by Oxford Nanopore Technologies.')
    parser.add('-fasta', '--reference_fasta_file', type=str, dest='fasta', required=True, default=None, help="The fasta format file describing the reference sequence for your organism.")
    parser.add('-targets', nargs = '*', dest='targets', required=True, help = 'Positional IDs to enrich for in the form seqid:start-stop . Can be space seperated eg: J02459:10000-15000 J02459:35000-40000')
    parser.add('-procs', '--proc_num', type=int, dest='procs',required=True, help = 'The number of processors to run this on.')
    parser.add('-m', '--model', type=str, required=True, help = 'The appropriate template model file to use', dest='temp_model')
    parser.add('-log', '--log-file', type=str, dest='logfile', default='readuntil.log', help="The name of the log file that data will be written to regarding the decision made by this program to process read until.")
    parser.add('-w', '--watch-dir', type=str, required=True, default=None, help="The path to the folder containing the downloads directory with fast5 reads to analyse - e.g. C:\data\minion\downloads (for windows).", dest='watchdir')
    parser.add('-length', '--library-length', type=int, dest='length', required=False, default=0, help="Provide the average expected length of your library. This offset will be applied to reads that are likely to extend into your region of interest on either strand.")
    parser.add('-o', '--output', type=str, required=False, default='test_read_until_out', help="Path to a folder to symbolically place reads representing match and not match.", dest='output_folder')
    parser.add('-v', '--verbose-true', action='store_true', help="Print detailed messages while processing files.", default=False, dest='verbose')
    parser.add_argument('-ver', '--version', action='version', version=('%(prog)s version={version} date={date}').format(version=__version__,date=__date__))
    args = parser.parse_args()

    check_files((args.fasta, args.temp_model))
    check_fasta(args.fasta)

    if not os.path.isdir(args.watchdir):
        logger.error("Folder" + args.watchdir + " cannot be found. Exiting.")
        sys.exit()

    # Multiprocess setup
    p = multiprocessing.Pool(args.procs)
    manager = multiprocessing.Manager()
    procampres = manager.dict()
    fasta_file = args.fasta
    seqlen = get_seq_len(fasta_file)

    # Process model and reference fasta
    model_file = args.temp_model
    global model_kmer_means
    global kmer_len
    model_kmer_means, kmer_len = process_model_file(model_file)
    seqids, threedarray = process_ref_fasta(fasta_file, model_kmer_means, kmer_len)

    # Scrap filenames
    d = list()
    filenamecounter = 0
    for filename in glob.glob(os.path.join(args.watchdir, '*.fast5')):
        filenamecounter += 1
        d.append([filename, seqids, threedarray, procampres, seqlen, args])
    for filename in glob.glob(os.path.join(args.watchdir, "pass",'*.fast5')):
        filenamecounter += 1
        d.append([filename, seqids, threedarray, procampres, seqlen, args])
    for filename in glob.glob(os.path.join(args.watchdir, "fail",'*.fast5')):
        filenamecounter += 1
        d.append([filename, seqids, threedarray, procampres, seqlen, args])
    procdata = tuple(d)

    # Assign process hdf5 to processes
    logger.info("Start spawning hdf5 processes")
    results = []
    for x in (procdata):
        r = p.apply_async(process_hdf5, args=(x,), callback=mycallback)
        # print(r.get())
        results.append(r)
    for r in results:
        r.wait()
    
    logger.info("Read Until completed. Exiting..")