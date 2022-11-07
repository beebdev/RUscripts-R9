import os
from re import T
import sys
import csv
# import ctypes
import time
import numpy as np
# import multiprocessing as mp
import sklearn.preprocessing as skprep
from Bio import SeqIO
import dtw


def file_existance_check(files):
    print("Checking files...", file=sys.stderr)
    for file in files:
        if not os.path.isfile(file):
            print("One or more files cannot be found. Exiting..", file=sys.stderr)
            print(file, file=sys.stderr)
            sys.exit()
    print("File check complete.", file=sys.stderr)


def validate_fasta_file(fasta):
    """
    Args:
        - fasta: path to the fasta file
    Returns:
        - None
    Notes:
        * fasta is the reference sequence in base space
    """
    try:
        for record in SeqIO.parse(fasta, "fasta"):
            print("Fasta record: ", record.id, file=sys.stderr)
    except:
        print("Invalid Fasta. Exiting...", file=sys.stderr)
        sys.exit()


def get_seq_len(ref_fasta):
    """
    Args:
        - None
    Returns:
        - seqLens: Length of each seq in fasta
    """
    seqLens = dict()
    for record in SeqIO.parse(ref_fasta, "fasta"):
        seq = record.seq
        seqLens[record.id] = len(seq)
    return seqLens


def process_model_file(model_file):
    """
    Args:
        - model: filepath to the model file
    Returns:
        - model_kmers: dict with kmer string as key and current value as value
        - kmer_len: length of the kmer string (Example run includes the b'' characters)
    Notes:
        * A model file is the mapping from base to squiggle space
    """
    print("Processing model file: ", model_file, file=sys.stderr)
    model_kmers = dict()
    with open(model_file) as file:
        csv_reader = csv.reader(file, delimiter="\t")
        raw_data = list(csv_reader)
        header = raw_data[0]
        kmer_columnID = header.index("kmer")
        mean_columnID = header.index("level_mean")

        data = np.array(raw_data[1:])
        model_kmers = dict(data[:, [kmer_columnID, mean_columnID]])
        kmer_len = len(data[0][kmer_columnID])
        # print(kmer_len)
    return model_kmers, kmer_len


def process_ref_fasta(ref_fasta, model_kmer_means, kmer_len):
    '''
    Description: Generates the reference sequence in squiggle space.
    Args:
        - ref_fasta: path to the fasta file
        - model_kmer_means: kmer dict that has mapping from base to squiggle space
        - kmer_len: length of each kmer string
    Returns:
        - seqIDs: A list of seqIDs
        - threedarrayshared_array: A 3D array of [seqID][F or R][ref seq]
    '''
    kmer_means = {}
    kmer_len -= 3

    for record in SeqIO.parse(ref_fasta, "fasta"):
        kmer_means[record.id] = {}
        forward = []  # Forward reference
        reverse = []  # Reverse reference
        print("Record ID: ", record.id, file=sys.stderr)
        print("Record length: ", len(record.seq), "bases", file=sys.stderr)

        seq = record.seq
        for i in range(len(seq) + 1 - kmer_len):
            # Slide through the base reference sequence
            kmer = "b'" + str(seq[i: i + kmer_len]).upper() + "'"
            # and append the squiggle value of each window
            forward.append(float(model_kmer_means[kmer]))
        print("Forward ok", file=sys.stderr)

        seq = record.seq.reverse_complement()
        for i in range(len(seq) + 1 - kmer_len):
            # Slide through the base reference sequence
            kmer = "b'" + str(seq[i: i + kmer_len]).upper() + "'"
            # and append the squiggle value of each window
            reverse.append(float(model_kmer_means[kmer]))
        print("Reverse ok", file=sys.stderr)

        kmer_means[record.id]["Fprime"] = skprep.scale(
            forward, axis=0, with_mean=True, with_std=True, copy=True
        )
        kmer_means[record.id]["Rprime"] = skprep.scale(
            reverse, axis=0, with_mean=True, with_std=True, copy=True
        )
        print("===", file=sys.stderr)

    # Iterator for kmer_means
    items = kmer_means.items()
    items_ = map(process_items, items)
    seqIDs, arrays = zip(*items_)
    # print(type(arrays), type(arrays[0]), type(arrays[0][0]), type(arrays[0][0][0]), file=sys.stderr)

    # 3d Array containing [nSeq][nLists][ListValues]
    # nSeq = len(seqIDs)
    # row, col = list(arrays)[0].shape
    # threedarray = mp.Array(ctypes.c_double, nSeq * row * col)
    # threedarrayshared_array = np.ctypeslib.as_array(threedarray.get_obj())

    # for array in arrays:
        # print(len(array), "->", len(array[0]), len(array[1]), type(array[0][0]))
    # threedarrayshared_array = np.array(arrays, dtype=np.float32)
    threedarrayshared_array = arrays
    return seqIDs, threedarrayshared_array


def process_items(d):
    '''
    Description: Helper function for processing kmer_means items
    Args:
        - d: 0 -> ID; 1 -> {Fprime, Rprime}
    Returns:
        - seqid
        - np.array(result): 2D np array of Fprime and Rprime values
    '''
    seqid = d[0]
    result = []
    for _, v in d[1].items():
        result.append(v)
    return seqid, np.array(result)


def squiggle_search(squiggle, seqIDs, threedarray):
    '''
    Description: Perform squiggle search on the reference sequence and decide to go or not
    Args:
        - squiggle: a list of 250 events to be searched
        - channel_id: channel of the nanopore (not used in the offline testing)
        - read_id:
        - args: arguments of the Read Until program (don't need)
        - seqIDs: a list of sequence IDs
        - threedarray: 3D array of the reference squiggle [seqID][F/R][each data]
        - seqLen: kmer string length
    Returns: seqmatchnameR, distanceR, frR, rsR, reR, qsR, qeR
    '''
    result = []
    blocksize = 200000
    overlap = blocksize - 500

    # Search in each seqIDs
    for ref in seqIDs:
        # Separate the Fprime and Rprime
        refID = seqIDs.index(ref)
        Rprime, Fprime = threedarray[refID]

        # Preprocess squiggle
        queryarray = skprep.scale(
            np.array(squiggle, dtype=float),
            axis=0,
            with_mean=True,
            with_std=True,
            copy=True,
        )

        # Start time
        tic = time.time()

        # Search on the forward reference
        refsubset = Fprime
        indexes = np.array(range(len(refsubset)))

        # Reference partitioned into blocks. Stride of ${overlap}
        subrefs = [refsubset[i: i + blocksize] for i in indexes[::overlap]]
        for blockID, ref_ in enumerate(subrefs):
            # Get pre-DTW time
            # tic = time.time()

            # Run DTW for squiggle and forward reference
            # TODO: change dtw_subsequence to a wrapper instead of mlpy
            dist, _, path = dtw.dtw_subsequence(queryarray, ref_)
            # print(dist)
            result.append(
                (
                    dist,   # unnormalised min-distance warp path between sequences
                    ref,    # seqID of current reference
                    "F",    # forward reference
                    (len(Rprime) - (path[1][0] + (blockID * overlap))),
                    (len(Rprime) - (path[1][-1] + (blockID * overlap))),
                    # path[0][0],     # start position for squiggle
                    # path[0][-1],    # end position for squiggle
                )
            )
            print("F", dist, path[1][-1]+ (blockID * overlap), file=sys.stderr)

        # Search on the reverse reference
        refsubset = Rprime
        subrefs = [refsubset[i: i + blocksize] for i in indexes[::overlap]]
        for blockID, ref_ in enumerate(subrefs):
            # Get pre-DTW time
            # tic = time.time()

            # Run DTW for squiggle and reverse reference
            # TODO: change dtw_subsequence to a wrapper instead of mlpy
            dist, cost, path = dtw.dtw_subsequence(queryarray, ref_)
            # Corrected for the fact that this is a reverse complement
            result.append(
                (
                    dist,
                    ref,
                    "R",
                    # start position for subref position
                    path[1][0] + (blockID * overlap),
                    # end position for subref position
                    path[1][-1] + (blockID * overlap),
                    # path[0][0],
                    # path[0][-1],
                )
            )
            print("R", dist, path[1][-1]+ (blockID * overlap), file=sys.stderr)

    total_time = time.time() - tic
    # Note first two elements flipped for return deliberately.
    distance, seqid, direction, refStart, refEnd = min(
        result, key=lambda result: result[0])
    # if (refStart == 8601):
    #     print(queryarray)
    return seqid, distance, direction, refStart, refEnd

    return seqid, distance, direction, refStart, refEnd, max_cost, total_time

def go_or_no(seqid, direction, position, seqlen, args):
    for sequence in args.targets:
        if args.verbose:
            print(sequence, file=sys.stderr)

def go_or_no(seqid, direction, position, args):
    for sequence in args.targets:
        start = int(float(sequence.split(":", 1)[1].split("-", 1)[0]))
        stop = int(float(sequence.split(":", 1)[1].split("-", 1)[1]))
        length = seqlen[seqid]

        if args.verbose:
            print(start, stop, length, file=sys.stderr)
            print(sequence.split(":", 1)[0], file=sys.stderr)
            print(type(seqid), file=sys.stderr)

        # We note that the average template read length is 6kb for the test lambda dataset.
        # Therefore we are interested in reads which start at least 3kb in advance of our position of interest
        balance = args.length / 2
        if seqid.find(sequence.split(":", 1)[0]) >= 0:
            if args.verbose:
                print("Found it", file=sys.stderr)

            if direction == "F":
                if args.verbose:
                    print("Forward Strand", file=sys.stderr)

                if position >= (start - balance) and position <= stop:
                    return "Sequence"
            elif direction == "R":
                if args.verbose:
                    print("Reverse Strand", file=sys.stderr)
                # We assume that coordinates are reported with respect to the forward strand regardless of
                # wether you are matching to forward or reverse.
                if position >= (start - balance) and position <= stop:
                    return "Sequence"
    return "Skip"
