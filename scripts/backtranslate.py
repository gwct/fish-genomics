#!/usr/bin/python3
###########################################################
# For rodent exomes, 2.2021
# Takes an amino acid alignment and the corresponding, unaligned
# CDS sequences and backtranslates the alignment based on 
# string indices.
###########################################################

import sys, os, core, seqparse as seq, argparse

###########################################################

def ntToCodon(nt_seq):
    assert len(nt_seq) % 3 == 0, "\nOUT OF FRAME NUCLEOTIDE SEQUENCE! " + str(len(nt_seq));
    codon_seq = [(nt_seq[i:i+3]) for i in range(0, len(nt_seq), 3)];
    return codon_seq;

###########################################################

parser = argparse.ArgumentParser(description="Exonerate command generator");
parser.add_argument("-aa", dest="aa_dir", help="Directory of input amino acid alignments.", default=False);
parser.add_argument("-nt", dest="cds_dir", help="Directory of unaligned nucleotide sequences corresponding to the AA alignments.", default=False);
parser.add_argument("-o", dest="outdir", help="Desired output directory for backtranslated alignments.", default=False);
parser.add_argument("--overwrite", dest="overwrite", help="If the output directory already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
args = parser.parse_args();
# IO options

if not args.aa_dir or not os.path.isdir(args.aa_dir):
    sys.exit( " * Error 1: An input directory of amino acid alignments must be defined with -aa.");
args.aa_dir = os.path.abspath(args.aa_dir);

if not args.cds_dir or not os.path.isdir(args.cds_dir):
    sys.exit( " * Error 2: An input directory of nucleotide sequences must be defined with -nt.");
args.cds_dir = os.path.abspath(args.cds_dir);

if not args.outdir:
    sys.exit( " * Error 3: An output directory must be defined with -o.");

if os.path.isdir(args.outdir) and not args.overwrite:
    sys.exit( " * Error 4: Output directory (-o) already exists! Explicity specify --overwrite to overwrite it.");

if not os.path.isdir(args.outdir):
    os.system("mkdir " + args.outdir);
args.outdir = os.path.abspath(args.outdir);
# IO option error checking

prequal_dir = False;
# Maybe add functionality for prequal filtered alignments later.

core.PWS("# " + core.getDateTime() + " Starting back translation.");
aa_files = [ f for f in os.listdir(args.aa_dir) if f.endswith(".fas") ];
num_files = len(aa_files);
# Read the AA alignment file names.

counter = 0;
for f in aa_files:
    print(f);
    # if not f.endswith(".fas"):
    #     continue;

    #if counter % 10 == 0:
    #    print(counter, "/", num_files);
    #print(counter)
    counter += 1;

    #pid = f.split("-")[0].replace(".fa", "");
    pid = f.split(".fa")[0];
    #if pid != "ENSMUSP00000021056":
    #    continue;
    # Get the protein id by splitting the file name by - and removing the extension.

    aa_file = os.path.join(args.aa_dir, f);
    cds_file = os.path.join(args.cds_dir, pid + "-cds.fa");

    if not os.path.isfile(cds_file):
        continue;

    outfilename = os.path.join(args.outdir, pid + "-cds.fa");
    # Assign the file names for the input AA, NT, and output NT sequences.

    aa_seqs = seq.fastaGetDict(aa_file);
    cds_seqs = seq.fastaGetDict(cds_file);
    # Read the sequences for the input AA and NT files.

    cds_alns = {};
    # The CDS align dict
    for sample in aa_seqs:
    # Go through every sample in the alignment.
        # while sample == ">mm10" and len(cds_seqs[sample]) % 3 != 0:
        #     print(cds_seqs[sample]);
        #     cds_seqs[sample] = cds_seqs[sample][:-1];
        # # Some of the mouse seqs from Ensembl are off...

        cur_aa_aln = aa_seqs[sample];
        cur_cds_aln = "";
        seq_index = 0;
        # Tracker and sequence vars

        cds_id = [ i for i in cds_seqs if sample.replace("_", "-") in i ];

        if len(cds_id) != 1:
            print(cds_id);
            sys.exit();
        cds_id = cds_id[0];

        codon_seq = ntToCodon(cds_seqs[cds_id]);
        # Convert the NT sequence to a list of codons.

        translated_codons = seq.yabt(codon_seq);

        if translated_codons != cur_aa_aln.replace("-", "").replace("X", ""):


            print(sample);
            print(cur_aa_aln.replace("-", ""));
            print(len(cur_aa_aln.replace("-", "")))
            print(translated_codons);
            print(len(translated_codons))
            print(codon_seq);
            print(len(codon_seq));
            print("--------------------------")

            continue;

        for i in range(len(cur_aa_aln)):
            if cur_aa_aln[i] == "-":
                cur_cds_aln += "---";
            elif cur_aa_aln[i] == "X":
                cur_cds_aln += "NNN";
                if codon_seq[seq_index] == "NNN":
                    seq_index += 1;
            else:
                cur_cds_aln += codon_seq[seq_index];
                seq_index += 1;
        # Build up the CDS alignment, adding a codon of gaps if the AA is a gap.

        cds_alns[cds_id] = cur_cds_aln;
        # Assign the current CDS alignment to the dictionary.

        assert len(cur_aa_aln) == len(cur_cds_aln) / 3, "\nINVALID CDS ALN LENGTH";
        # Make sure the CDS alignment and AA alignment are the same length.

    with open(outfilename, "w") as outfile:
        for sample in cds_alns:
            outfile.write(sample + "\n");
            outfile.write(cds_alns[sample] + "\n");
    # Write out the current alignment.
core.PWS("# " + core.getDateTime() + " Done!");






