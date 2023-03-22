#############################################################################
# Pulls CDS sequences from a ortholog peptide sequences
#
# Authors: Gregg Thomas and Clara Boothby (regex)
#############################################################################

import sys
import os
import re
import core as CORE
import seqparse as SEQ
import gxfparse as GXF
from collections import defaultdict

#############################################################################


specs = ["allis_shad", "atlantic_herring", "cod", "denticle_herring", "milkfish", "shad"];
#specs = ["allis_shad"];
# Species list

seqs = { s : {'pep' : {}, 'cds' : {}, 'ids' : {}} for s in specs }
# A dict to store the sequences and ids for each species in memory

genome_dir = "/n/holylfs05/LABS/informatics/Users/gthomas/fish/genomes/";
# Directory with genomes in sub-folders listed in spec
# Each genome should have a file with CDS sequences and a file with peptide sequences

ortholog_dir = "/n/holylfs05/LABS/informatics/Users/gthomas/fish/Shad-Single_Copy_Orthologue_Sequences-031923";
# The predicted orthogroups from fastortho

outdir = "/n/holylfs05/LABS/informatics/Users/gthomas/fish/seq/";
# The output directory for the run

outfilename = "/n/home07/gthomas/projects/fish-genomics/data/multi-spec-orthogroups.txt";
# A file to track the orthogroups that have at least 4 taxa represented for Guidance

logfilename = "/n/home07/gthomas/projects/fish-genomics/data/get-orths-" + CORE.getLogTime() + ".log";
# A log file for the run

## I/O options
####################

pad = 20;
with open(logfilename, "w") as logfile:
    CORE.runTime("# Get sequences", logfile);
    CORE.PWS(CORE.spacedOut("# Ortholog dir:", pad) + ortholog_dir, logfile);
    CORE.PWS(CORE.spacedOut("# Genome dir:", pad) + genome_dir, logfile);
    CORE.PWS(CORE.spacedOut("# Output dir:", pad) + outdir, logfile);
    CORE.PWS(CORE.spacedOut("# Log file:", pad) + logfilename, logfile);
    CORE.PWS("# ----------------", logfile);
    # Log info

    CORE.PWS("# " + CORE.getDateTime() + " Reading sequences", logfile);

    for spec in specs:
        CORE.PWS("# " + CORE.getDateTime() + " " + spec, logfile);

        #pep_file = [ f for f in os.listdir(os.path.join(genome_dir, spec)) if "protein" in f ][0];
        #pep_file = os.path.join(genome_dir, spec, pep_file);
        # Peptide file

        cds_file = [ f for f in os.listdir(os.path.join(genome_dir, spec)) if "cds_from_genomic" in f ][0];
        cds_file = os.path.join(genome_dir, spec, cds_file);
        # Lookup CDS file

        ## Get files for the current species
        ##########

        #CORE.PWS("# " + CORE.getDateTime() + "    Reading peptides   : " + pep_file, logfile);
        #seqs[spec]['pep'] = SEQ.fastaReadSeqs(pep_file, header_sep=" ");
        #CORE.PWS("# " + CORE.getDateTime() + "    Success            : " + str(len(seqs[spec]['pep'])) + " seqs read", logfile);
        CORE.PWS("# " + CORE.getDateTime() + "    Reading CDS        : " + cds_file, logfile);
        in_seqs = SEQ.fastaReadSeqs(cds_file);
        CORE.PWS("# " + CORE.getDateTime() + "    Success            : " + str(len(in_seqs)) + " seqs read", logfile);
        ## Read sequences for the current species
        ##########

        # CORE.PWS("# " + CORE.getDateTime() + "    Reading annotations: " + gff_file, logfile);
        # cur_exons, comp, tid_to_gid, pid_to_tid = GXF.getExons(gff_file, coding_only=True);
        ## Read annotations for the current species
        ##########

        missing_pid = 0;

        for seq in in_seqs:
            pid_part = re.findall("protein_id=XP_[0-9.]+", seq);

            if pid_part:
                pid = pid_part[0].replace("protein_id=", "");

                seqs[spec]['cds'][pid] = in_seqs[seq];

            else:
                missing_pid += 1;

        CORE.PWS("# " + CORE.getDateTime() + "    Success            : " + str(len(seqs[spec]['cds'])) + " seqs retained, " + str(missing_pid) + " seqs missing protein IDs", logfile);
        # Save the ids for the current species
    ## End spec loop

    ## End sequence read block
    ####################

    CORE.PWS("# " + CORE.getDateTime() + " Reading orthologs and writing sequences", logfile);
    #num_clusters = "41448"; # 16 spec
    #num_clusters = 49561;
    # The total number of clusters from the i3 file (if we want to track progress)

    num_written = 0;
    num_written_d = defaultdict(int);
    # Number of orthogroups written to file

    incorrect_num_seqs = [];
    # Orthogroups with the incorrect number of sequences recovered

    seqs_w_stop = 0;
    seqs_w_stop_d = defaultdict(int);
    # Number of CDS sequences with premature stop codons
    ## Some trackers
    ##########
   

    match_code_counts = defaultdict(int);
    # Counts of each match code case
    # Match codes are for different scenarios of CDS to protein sequence matches or mis-matches (e.g. missing one or more codons)

    with open(outfilename, "w") as outfile:
    # Open the file to track orthogroups with at least 4 taxa

        #i = 1;

        og_4taxa = 0;
        # Orthogroups with at least 4 species

        lt_4taxa_after = 0;
        # Orthogroups with less than 4 species after filtering

        non_div3 = 0;
        non_div3_d = defaultdict(int);
        # CDS Sequences that aren't divisible by 3

        multi_cds = 0;
        # Preteins with multiple CDS recovered

        wrong_translation = 0;
        wrong_translation_d = defaultdict(int);
        # CDS sequences whose translation does not match the peptide sequence

        total_seqs_pre = 0;
        total_seqs_pre_d = defaultdict(int);
        # The total number of sequences before any filtering

        ## Some more trackers
        ##########

        for f in os.listdir(ortholog_dir):
        # Go through every ortholog cluster in the file

            ortho_file = os.path.join(ortholog_dir, f);
            # Get the ortho file name

            cluster_id = f.replace(".fa", "");
            #print(cluster_id)
            # Extract the cluster id

            ortho_seqs = SEQ.fastaReadSeqs(ortho_file);
            # Read the sequences of the current ortholog

            expected_seqs = len(ortho_seqs);
            total_seqs_pre += expected_seqs;
            # Extract the number of sequences in the cluster          

            expected_taxa = len(set( [ s.split("-")[0] for s in ortho_seqs ] ));
            # The number of taxa in the current orthogroup

            if expected_taxa >= 4:
                og_4taxa += 1;
            # Check if the current cluster has at least 4 taxa represented and write it to the file if so

            cur_taxa_written = [];

            #print(str(i) + " / " + num_clusters + " " + cluster_id);
            #i += 1;
            # Print a status update

            #protein_ids = re.sub("\(([^)]+)\)", "", line[1]).split(" ");
            #protein_ids = list(filter(None, protein_ids));
            # Remove the file name from the protein ids

            cds_outfile = os.path.join(outdir, "cds", cluster_id + "-cds.fa");
            # Set the output files

            cur_num_written = 0;
            # Tracking for the current cluster

            with open(cds_outfile, "w") as cdsout:
            # Open the output files

                for seq in ortho_seqs:
                    spec, pid = seq.split("-");
                    pep_seq = ortho_seqs[seq];

                    #print(spec, pid);
                    #sys.exit();

                    #print(spec, pid);

                    total_seqs_pre_d[spec] += 1;

                    if re.search("-PA.[\d]+", pid):
                        pid = pid[:-2];
                    # Remove any version numbers from the protein id (e.g. .1)

                    # pep_seq = seqs[spec]['pep'][orth];
                    # # Lookup the protein sequence

                    # tid = seqs[spec]['ids'][pid];
                    # # Lookup the transcript/CDS id

                    # cds_header = [key for key, val in seqs[spec]['cds'].items() if tid in key];
                    # # Lookup the CDS header by looking for headers that contain the current transcript id

                    # ####################
                    # if len(cds_header) != 1:
                    # # Ideally, there would be only one header per transcript, but that isn't the case for some (8)
                    # # sequences... Skip these for now
                    #     CORE.PWS("# " + CORE.getDateTime() + " ********** ", logfile);
                    #     CORE.PWS("# " + CORE.getDateTime() + " Multiple headers found for transcript ID: " + tid, logfile);
                    #     CORE.PWS("# " + CORE.getDateTime() + " Species: " + spec, logfile);
                    #     CORE.PWS("# " + CORE.getDateTime() + " Associated protein ID: " + pid, logfile);
                    #     CORE.PWS("# " + CORE.getDateTime() + " Ortho cluster: " + cluster_id, logfile);
                    #     CORE.PWS("# " + CORE.getDateTime() + " Headers: " + str(cds_header), logfile);
                    #     CORE.PWS("# " + CORE.getDateTime() + " ********** ", logfile);
                    #     multi_cds += 1;
                    #     continue;
                    # ## End multi-header check
                    ####################

                    # cds_header = cds_header[0];
                    # # Convert the header to a string

                    # if spec in ['abrue', 'bmori', 'crotu', 'cscul', 'dmela', 'hlong', 'iscap', 'lhesp', 'lpoly', 'lrecl', 'mocci', 'nclav', 'ptepi', 'smimo', 'sscab', 'turti', 'vdest']:
                    #     cds_output_header = cds_header.split(" ")[0];
                    # elif spec in ['tanti']:
                    #     cds_output_header = cds_header.split("_")[2];
                    # cds_output_header = cds_output_header.replace("_", "-");
                    # # Get the header for the output CDS file

                    cds_seq = seqs[spec]['cds'][pid];
                    # Lookup the CDS sequence

                    cds_seq_orig = cds_seq;
                    codon_seq_orig = SEQ.ntToCodon(cds_seq);
                    translated_seq_orig = SEQ.yabt(codon_seq_orig);
                    # The original codon and translated seqs for reference later if needed.

                    final_cds, final_pep = "", "";

                    ####################
                    cds_seq_front = "";
                    cds_seq_end = "";
                    # In the case that cds_seq is not divisible by 3, we check the frame in both the front and end of the sequence

                    if len(cds_seq) % 3 != 0:
                    # Check to see if the CDS is divisible by 3

                        non_div3 += 1;

                        cds_seq_front = cds_seq;
                        while len(cds_seq_front) % 3 != 0:
                            cds_seq_front = cds_seq_front[1:];
                        # Remove 1 or 2 bases from the FRONT of the CDS sequence until it is divisible by 3

                        cds_seq_end = cds_seq;
                        while len(cds_seq_end) % 3 != 0:
                            cds_seq_end = cds_seq_end[:-1];
                        # Remove 1 or 2 bases from the END of the CDS sequence until it is divisible by 3

                        cds_seq = "";
                        # Save the original CDS but remove it from the cds_seq var for later
                    ## End CDS div 3 check
                    ####################

                    codon_seq = SEQ.ntToCodon(cds_seq);
                    translated_seq = SEQ.yabt(codon_seq);
                    if translated_seq and translated_seq[-1] == "X":
                        codon_seq = codon_seq[:-1];
                        translated_seq = translated_seq[:-1];
                    # Translate the original CDS sequence if it was divisible by 3
                    # If not, these vars will be empty strings

                    ####################

                    translated_seq_no_x = translated_seq.replace("X", "");
                    codon_seq_no_x = [ codon_seq[i] for i in range(len(codon_seq)) if translated_seq[i] != "X" ];

                    codon_seq_front = SEQ.ntToCodon(cds_seq_front);
                    translated_seq_front = SEQ.yabt(codon_seq_front);
                    if translated_seq_front and translated_seq_front[-1] == "X":
                        codon_seq_front = codon_seq_front[:-1];
                        translated_seq_front = translated_seq_front[:-1];
                    # Translate the CDS sequence with bases removed from the FRONT in the case that it was not divisible by 3
                    # Otherwise, these vars will be empty strings

                    ####################

                    translated_seq_front_no_x = translated_seq_front.replace("X", "");
                    codon_seq_front_no_x = [ codon_seq_front[i] for i in range(len(codon_seq_front)) if translated_seq_front[i] != "X" ];
                    # Remove Xs from the front of the CDS sequence

                    codon_seq_end = SEQ.ntToCodon(cds_seq_end);
                    translated_seq_end = SEQ.yabt(codon_seq_end);
                    if translated_seq_end and translated_seq_end[-1] == "X":
                        codon_seq_end = codon_seq_end[:-1];
                        translated_seq_end = translated_seq_end[:-1];
                    # Translate the CDS sequence with bases removed from the END in the case that it was not divisible by 3
                    # Otherwise, these vars will be empty strings

                    translated_seq_end_no_x = translated_seq_end.replace("X", "");
                    codon_seq_end_no_x = [ codon_seq_end[i] for i in range(len(codon_seq_end)) if translated_seq_end[i] != "X" ];
                    # Remove Xs from the end of the CDS sequence

                    pep_seq_no_x = pep_seq.replace("X", "");
                    # Remove Xs from the peptide sequence since they may not be reflected in the CDS sequence

                    ####################

                    cur_seqs = {    
                                    'front' : {'codons' : codon_seq_front, 'translated' : translated_seq_front, "pep" : pep_seq },
                                    'front.no.x' : {'codons' : codon_seq_front_no_x, 'translated' : translated_seq_front_no_x, "pep" : pep_seq_no_x },
                                    'orig' : {'codons' : codon_seq, 'translated' : translated_seq, "pep" : pep_seq },
                                    'orig.no.x' : {'codons' : codon_seq_no_x, 'translated' : translated_seq_no_x, "pep" : pep_seq_no_x },
                                    'end' : {'codons' : codon_seq_end, 'translated' : translated_seq_end, "pep" : pep_seq },
                                    'end.no.x' : {'codons' : codon_seq_end_no_x, 'translated' : translated_seq_end_no_x, "pep" : pep_seq_no_x },
                                }
                    # Save the possible sequences in this dict to loop over and check for other errors below

                    ####################

                    match_found = False;
                    # Tells us whether we find a match between the CDS and the peptide sequence

                    match_code = "";
                    # A list that stores how the sequences were matched

                    final_cds, final_pep = "", "";

                    for seq_type in cur_seqs:
                    # Loop over each type of CDS sequence: original, bases removed from FRONT, and bases removed from END

                        if cur_seqs[seq_type]['translated'] == "" or match_found:
                            continue;
                        # If the sequence doesn't exist because it wasn't a case for this sample, skip                        

                        if cur_seqs[seq_type]['translated'] == pep_seq:
                            final_cds = "".join(cur_seqs[seq_type]['codons']);
                            final_pep = pep_seq;
                            match_code = seq_type + "-match";
                            match_found = True;
                        # If the translated sequence matches the peptide sequence, save those as final here

                        ####################

                        else:
                        # If the translated sequence and peptide sequence don't match, try a bunch of things to see if we can find a match

                            if not match_found and cur_seqs[seq_type]['translated'] in pep_seq:
                                final_cds = "".join(cur_seqs[seq_type]['codons']);
                                # The final CDS seq is the original CDS seq
                                
                                pep_start = pep_seq.index(cur_seqs[seq_type]['translated']);
                                pep_end = len(cur_seqs[seq_type]['translated']) + pep_start;
                                final_pep = pep_seq[pep_start:pep_end];
                                # For the final peptide seq, we find the start and end coords of the translated seq in the peptide seq

                                match_code = seq_type + "-trans.substr.pep";
                                match_found = True;
                                CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);
                                CORE.PWS("# " + CORE.getDateTime() + " Translated CDS sequence within peptide seq: " + " " + match_code + " " + cluster_id + " " + spec + " " + pid, logfile);
                                CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);
                                # Log info
                            # Check if the translated seq is a subtring of the peptide seq

                            ####################

                            if not match_found and pep_seq in cur_seqs[seq_type]['translated']:
                                codon_start = cur_seqs[seq_type]['translated'].index(pep_seq);
                                codon_end = len(pep_seq) + codon_start;
                                final_cds = "".join(cur_seqs[seq_type]['codons']);
                                # The final CDS seq is the original CDS seq for this sequence type

                                final_pep = pep_seq;
                                # For the final peptide seq, we find the start and end coords of the peptdie seq in the translated seq

                                match_code = seq_type + "-pep.substr.trans"
                                match_found = True;
                                CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);
                                CORE.PWS("# " + CORE.getDateTime() + " Peptide sequence within translated seq: " + " " + match_code + " " + cluster_id + " " + spec + " " + pid, logfile);
                                CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);
                                # Log info
                            ## Check if the peptide sequnce is a substring of the translated CDS seq
                        ## End mis-match cases block
                    ## End CDS seq type loop

                    ####################

                    if not match_found:
                        CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);
                        CORE.PWS("# " + CORE.getDateTime() + " Translated CDS sequence doesn't match peptide sequence " + " " + cluster_id + " " + spec + " " + pid, logfile);
                        CORE.PWS("# " + CORE.getDateTime() + " " + str(len(cds_seq_orig)), logfile);
                        CORE.PWS("# " + CORE.getDateTime() + " " + cds_seq_orig, logfile);
                        CORE.PWS("# " + CORE.getDateTime() + " " + str(len(translated_seq_orig)), logfile);
                        CORE.PWS("# " + CORE.getDateTime() + " " + translated_seq_orig, logfile);
                        CORE.PWS("# " + CORE.getDateTime() + " " + str(len(pep_seq)), logfile);
                        CORE.PWS("# " + CORE.getDateTime() + " " + pep_seq, logfile);
                        CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);
                        wrong_translation += 1;
                        wrong_translation_d[spec] += 1;

                        # for seq_type in cur_seqs:
                        #     print(seq_type);
                        #     for seq in cur_seqs[seq_type]:
                        #         print("\t" + seq + "\t" + str(cur_seqs[seq_type][seq]));
                        # sys.exit();
                    ## Report if even after all match case checks the sequences still don't match

                    elif not final_cds or not final_pep:
                        sys.exit("something went wrong");
                    # If the code thinks a match was found but doesn't have the sequences, exit
                    # This should never happen...

                    else:

                        match_code_counts[match_code] += 1;
                        final_codons = SEQ.ntToCodon(final_cds);
                        final_translated_seq = SEQ.yabt(final_codons);

                        if final_translated_seq != final_pep:
                            CORE.PWS("# " + CORE.getDateTime() + " +++++++++++++++++++++++++++++++++++", logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " Translated CDS sequence doesn't match peptide sequence AFTER CORRECTION " + " " + cluster_id + " " + spec + " " + pid, logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + match_code, logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + str(len(cds_seq_orig)), logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + cds_seq_orig, logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + str(len(translated_seq_orig)), logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + translated_seq_orig, logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + str(len(pep_seq)), logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + pep_seq, logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " ------", logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + str(len(final_cds)), logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + final_cds, logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + str(len(final_translated_seq)), logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + final_translated_seq, logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + str(len(final_pep)), logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + final_pep, logfile);                            
                            CORE.PWS("# " + CORE.getDateTime() + " +++++++++++++++++++++++++++++++++++", logfile);
                            sys.exit();
                        ## One final check to make sure the final translated CDS sequence and final peptide sequence match                         

                        if len(final_cds) % 3 != 0:
                            CORE.PWS("# " + CORE.getDateTime() + " +++++++++++++++++++++++++++++++++++", logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " Final CDS not divisible by 3 AFTER CORRECTION " + " " + cluster_id + " " + spec + " " + pid, logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + match_code, logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + str(len(cds_seq_orig)), logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + cds_seq_orig, logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + str(len(translated_seq_orig)), logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + translated_seq_orig, logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + str(len(pep_seq)), logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + pep_seq, logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " ------", logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + str(len(final_cds)), logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + final_cds, logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + str(len(final_translated_seq)), logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + final_translated_seq, logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + str(len(final_pep)), logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + final_pep, logfile);                            
                            CORE.PWS("# " + CORE.getDateTime() + " +++++++++++++++++++++++++++++++++++", logfile);
                            sys.exit();
                        ## One final check to make sure the final CDS sequence is divisible by 3                                 

                        if any(stop in final_codons for stop in ['TAA', 'TAG', 'TGA']):
                            CORE.PWS("# " + CORE.getDateTime() + " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " Stop codon within final coding sequence: " + " " + cluster_id + " " + spec + " " + pid, logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + str(len(final_cds)), logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + final_cds, logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + str(len(final_translated_seq)), logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " " + final_translated_seq, logfile);
                            CORE.PWS("# " + CORE.getDateTime() + " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", logfile);
                            seqs_w_stop += 1;
                            seqs_w_stop_d[spec] += 1;
                            continue;

                        # pepout.write(">" + orth + "\n");
                        # pepout.write(final_pep + "\n");
                        # Write out the protein sequence

                        cdsout.write(">" + spec + "-" + pid + "\n");
                        cdsout.write(final_cds + "\n");
                        # Write out the CDS sequence

                        cur_num_written += 1;
                        num_written += 1;
                        cur_taxa_written.append(spec);
                        num_written_d[spec] += 1;
                        # Update trackers
                ## End ortho cluster loop
            ## Close output files
            ####################

            if cur_num_written != expected_seqs:
                incorrect_num_seqs.append(cluster_id + ": " + str(cur_num_written) + "/" + str(expected_seqs));
            # Check if the number of sequences written matches the number expected for the current cluster

            if len(set(cur_taxa_written)) < 4:
                CORE.PWS("# " + CORE.getDateTime() + " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^", logfile);
                CORE.PWS("# " + CORE.getDateTime() + " Orthogroup has fewer than 4 taxa after filtering: " + " " + cluster_id, logfile);
                CORE.PWS("# " + CORE.getDateTime() + " Removing CDS sequence file: " + " " + cds_outfile, logfile);
                os.remove(cds_outfile);
                #CORE.PWS("# " + CORE.getDateTime() + " Removing peptide sequence file: " + " " + pep_outfile, logfile);
                #os.remove(pep_outfile);
                CORE.PWS("# " + CORE.getDateTime() + " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^", logfile);
                lt_4taxa_after += 1;
            else:
                outfile.write(cluster_id + "\n");
            # Check if sequences from at least 4 taxa have been written to the output file
        ## End ortholog loop
    ## Close main output file

    ####################

    CORE.PWS("# " + CORE.getDateTime() + " Done!", logfile);
    CORE.PWS("# " + CORE.getDateTime() + " Clusters with at least 4 taxa: " + str(og_4taxa), logfile);
    CORE.PWS("# " + CORE.getDateTime() + " Total sequences in clusters: " + str(total_seqs_pre), logfile);
    for spec in total_seqs_pre_d:
        CORE.PWS("# " + CORE.getDateTime() + " " + spec + " " + str(total_seqs_pre_d[spec]), logfile);
    CORE.PWS("# ----------------", logfile);

    CORE.PWS("# " + CORE.getDateTime() + " Total sequences written: " + str(num_written), logfile);
    for spec in num_written_d:
        CORE.PWS("# " + CORE.getDateTime() + " " + spec + " " + str(num_written_d[spec]), logfile);
    CORE.PWS("# ----------------", logfile);

    if incorrect_num_seqs:
        CORE.PWS("# " + CORE.getDateTime() + " " + str(len(incorrect_num_seqs)) + " clusters didn't have the correct number of sequences written:", logfile);
        for seq in incorrect_num_seqs:
            CORE.PWS("# " + CORE.getDateTime() + "\t" + seq, logfile);
    CORE.PWS("# ----------------", logfile);

    if multi_cds:
        CORE.PWS("# " + CORE.getDateTime() + " " + str(multi_cds) + " sequences had multiple headers. See above.", logfile);
    CORE.PWS("# ----------------", logfile);

    if non_div3:
        CORE.PWS("# " + CORE.getDateTime() + " " + str(non_div3) + " CDS sequences not divisible by 3.", logfile);
        for spec in non_div3_d:
            CORE.PWS("# " + CORE.getDateTime() + " " + spec + " " + str(non_div3_d[spec]), logfile);
    CORE.PWS("# ----------------", logfile);

    if wrong_translation:
        CORE.PWS("# " + CORE.getDateTime() + " " + str(wrong_translation) + " CDS sequences have translations that don't match the peptide sequence. See above.", logfile);
        for spec in wrong_translation_d:
            CORE.PWS("# " + CORE.getDateTime() + " " + spec + " " + str(wrong_translation_d[spec]), logfile);
    CORE.PWS("# ----------------", logfile);

    if seqs_w_stop:
        CORE.PWS("# " + CORE.getDateTime() + " " + str(seqs_w_stop) + " CDS sequences have stop codons. See above.", logfile);
        for spec in seqs_w_stop_d:
            CORE.PWS("# " + CORE.getDateTime() + " " + spec + " " + str(seqs_w_stop_d[spec]), logfile);
    CORE.PWS("# ----------------", logfile);

    CORE.PWS("# " + CORE.getDateTime() + " Match code case counts:", logfile);
    for match_case in match_code_counts:
        CORE.PWS("# " + CORE.getDateTime() + " " + match_case + " " + str(match_code_counts[match_case]), logfile);
    CORE.PWS("# ----------------", logfile);

    if lt_4taxa_after:
        CORE.PWS("# " + CORE.getDateTime() + " " + str(lt_4taxa_after) + " orthogroups with less than 4 genes AFTER correction", logfile);
    CORE.PWS("# ----------------", logfile);
## Close log file


#############################################################################










                    #     if "[partial=3']" in cds_header:
                    #         CORE.PWS("# " + CORE.getDateTime() + " --------------------------------", logfile);
                    #         CORE.PWS("# " + CORE.getDateTime() + " PARTIAL 3' CDS SEQUENCE ENDS OUT OF FRAME. ADJUSTING BY ADDING Ns ", logfile);
                    #         CORE.PWS("# " + CORE.getDateTime() + " " + cluster_id + " " + spec + " " + cds_header + " " + pid, logfile);
                    #         CORE.PWS("# " + CORE.getDateTime() + " " + str(len(cds_seq)), logfile);

                    #         while len(cds_seq % 3) != 0:
                    #             cds_seq += "N";

                    #         CORE.PWS("# " + CORE.getDateTime() + " " + cds_seq, logfile);


                    #     elif "[partial=5']" in cds_header:
                    #         CORE.PWS("# " + CORE.getDateTime() + " --------------------------------", logfile);
                    #         CORE.PWS("# " + CORE.getDateTime() + " PARTIAL 5' CDS SEQUENCE STARTS OUT OF FRAME. ADJUSTING BY ADDING Ns ", logfile);
                    #         CORE.PWS("# " + CORE.getDateTime() + " " + cluster_id + " " + spec + " " + cds_header + " " + pid, logfile);
                    #         CORE.PWS("# " + CORE.getDateTime() + " " + str(len(cds_seq)), logfile);

                    #         while len(cds_seq % 3) != 0:
                    #             cds_seq = "N" + cds_seq;

                    #         CORE.PWS("# " + CORE.getDateTime() + " " + cds_seq, logfile);

                    #     elif "[partial=5',3']" in cds_header:
                    #         pass;

                        # corrected = False;
                        # cds_seq_orig = cds_seq;
                        # bases_rm_end, bases_rm_front = 0, 0;
                        # while len(cds_seq) % 3 != 0:
                        #     cds_seq = cds_seq[:-1];
                        #     bases_rm_end += 1;

                        # codon_seq = SEQ.ntToCodon(cds_seq);
                        # translated_seq = SEQ.yabt(codon_seq);

                        # if translated_seq == pep_seq:
                        #     CORE.PWS("# " + CORE.getDateTime() + " --------------------------------", logfile);
                        #     CORE.PWS("# " + CORE.getDateTime() + " NON-DIVISIBLE BY 3 CDS CORRECTED BY REMOVING " + str(bases_rm_end) + " BASES FROM END OF SEQ: " + " " + cluster_id + " " + spec + " " + cds_header + " " + pid, logfile);
                        #     CORE.PWS("# " + CORE.getDateTime() + " --------------------------------", logfile);    
                        #     corrected = True;
                        # elif translated_seq in pep_seq:
                        #     cds_start = pep_seq.index(translated_seq);
                        #     cds_end = len(translated_seq) + cds_start;


                            
                        #     CORE.PWS("# " + CORE.getDateTime() + " --------------------------------", logfile);
                        #     CORE.PWS("# " + CORE.getDateTime() + " TRANSLATED SEQ IN PEPTIDE SEQ AFTER NON-DIVISIBLE BY 3 CDS CORRECTED BY REMOVING " + str(bases_rm_end) + " BASES FROM END OF SEQ: " + " " + cluster_id + " " + spec + " " + cds_header + " " + pid, logfile);
                        #     CORE.PWS("# " + CORE.getDateTime() + " --------------------------------", logfile);    
                        #     corrected = True;

                        # elif pep_seq in translated_seq:
                        #     CORE.PWS("# " + CORE.getDateTime() + " --------------------------------", logfile);
                        #     CORE.PWS("# " + CORE.getDateTime() + " PEPTIDE SEQ IN TRANSLATED SEQ AFTER NON-DIVISIBLE BY 3 CDS CORRECTED BY REMOVING " + str(bases_rm_end) + " BASES FROM END OF SEQ: " + " " + cluster_id + " " + spec + " " + cds_header + " " + pid, logfile);
                        #     CORE.PWS("# " + CORE.getDateTime() + " --------------------------------", logfile);   

                        # else:
                        #     cds_seq = cds_seq_orig;
                        #     while len(cds_seq) % 3 != 0:
                        #         cds_seq = cds_seq[1:];
                        #         bases_rm_front += 1;

                        #     codon_seq = SEQ.ntToCodon(cds_seq);
                        #     translated_seq = SEQ.yabt(codon_seq);
                            
                        #     if translated_seq == pep_seq or translated_seq in pep_seq or pep_seq in translated_seq:
                        #         CORE.PWS("# " + CORE.getDateTime() + " --------------------------------", logfile);
                        #         CORE.PWS("# " + CORE.getDateTime() + " NON-DIVISIBLE BY 3 CDS CORRECTED BY REMOVING " + str(bases_rm_front) + " BASES FROM FRONT OF SEQ: " + " " + cluster_id + " " + spec + " " + cds_header + " " + pid, logfile);
                        #         CORE.PWS("# " + CORE.getDateTime() + " --------------------------------", logfile);    
                        #         corrected = True;

                        # if not corrected:
                        #     CORE.PWS("# " + CORE.getDateTime() + " --------------------------------", logfile);
                        #     CORE.PWS("# " + CORE.getDateTime() + " NON-DIVISIBLE BY 3 CDS: " + " " + cluster_id + " " + spec + " " + cds_header + " " + pid, logfile);
                        #     CORE.PWS("# " + CORE.getDateTime() + " " + str(len(cds_seq_orig)), logfile);
                        #     CORE.PWS("# " + CORE.getDateTime() + " " + cds_seq_orig, logfile);
                        #     CORE.PWS("# " + CORE.getDateTime() + " " + str(len(pep_seq)), logfile);
                        #     CORE.PWS("# " + CORE.getDateTime() + " " + pep_seq, logfile);
                        #     CORE.PWS("# " + CORE.getDateTime() + " --------------------------------", logfile);
                        #     non_div3 += 1;
                        #     non_div3_d[spec] += 1;
                        #     continue;





                    # if translated_seq != pep_seq:

                    #     corrected = False;
                    #     codon_seq_orig = codon_seq;
                    #     translated_seq_orig = translated_seq

                    #     if not corrected and translated_seq.replace("X", "") == pep_seq.replace("X", ""):
                    #         CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);
                    #         CORE.PWS("# " + CORE.getDateTime() + " Translated CDS sequence corrected by removing X's from seqs " + " " + cluster_id + " " + spec + " " + cds_header + " " + pid, logfile);
                    #         CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);   
                    #         corrected = True;


                        # if not corrected:
                        #     codon_rm_end = 0;
                        #     while codon_seq:
                        #         codon_seq = codon_seq[:-1];
                        #         translated_seq = SEQ.yabt(codon_seq);
                        #         codon_rm_end += 1;
                        #         if translated_seq == pep_seq:
                        #             CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);
                        #             CORE.PWS("# " + CORE.getDateTime() + " Translated CDS sequence corrected by removing " + str(codon_rm_end) + " codons from end of seq " + " " + cluster_id + " " + spec + " " + cds_header + " " + pid, logfile);
                        #             CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);   
                        #             corrected = True;

                        # if not corrected:
                        #     codon_rm_front = 0;
                        #     codon_seq = codon_seq_orig;
                        #     while codon_seq:
                        #         codon_seq = codon_seq[1:];
                        #         translated_seq = SEQ.yabt(codon_seq);
                        #         codon_rm_front += 1;
                        #         if translated_seq == pep_seq:
                        #             CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);
                        #             CORE.PWS("# " + CORE.getDateTime() + " Translated CDS sequence corrected by removing " + str(codon_rm_front) + " codons from front of seq " + " " + cluster_id + " " + spec + " " + cds_header + " " + pid, logfile);
                        #             CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);   
                        #             corrected = True;                                                   

                        # if not corrected:
                            # if translated_seq_orig in pep_seq:
                            #     CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);
                            #     CORE.PWS("# " + CORE.getDateTime() + " Translated CDS sequence within peptide seq: " + " " + cluster_id + " " + spec + " " + cds_header + " " + pid, logfile);
                            #     CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);   
                            #     corrected = True;
                            # elif translated_seq_orig.replace("X", "") in pep_seq.replace("X", ""):                     
                            #     CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);
                            #     CORE.PWS("# " + CORE.getDateTime() + " Translated CDS sequence within peptide seq AFTER removing X's: " + " " + cluster_id + " " + spec + " " + cds_header + " " + pid, logfile);
                            #     CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);   
                            #     corrected = True;
                            # elif pep_seq in translated_seq_orig:
                            #     CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);
                            #     CORE.PWS("# " + CORE.getDateTime() + " Peptide sequence within translated seq: " + " " + cluster_id + " " + spec + " " + cds_header + " " + pid, logfile);
                            #     CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);   
                            #     corrected = True;     
                            # elif pep_seq.replace("X", "") in translated_seq_orig.replace("X", ""):
                            #     CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);
                            #     CORE.PWS("# " + CORE.getDateTime() + " Peptide sequence within translated seq AFTER removing X's: " + " " + cluster_id + " " + spec + " " + cds_header + " " + pid, logfile);
                            #     CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);   
                            #     corrected = True;                            


                            # if not match_found and cur_seqs[seq_type]['translated'].replace("X", "") == pep_seq.replace("X", ""):

                            #     final_codons = [];
                            #     for i in range(len(cur_seqs[seq_type]['translated'])):
                            #         if cur_seqs[seq_type]['translated'][i] != "X":
                            #             final_codons.append(cur_seqs[seq_type]['codons'][i]);
                            #     final_cds = "".join(final_codons);
                            #     # For the final CDS, we have to join the codons that don't result in X's in the translated seq

                            #     final_pep = pep_seq.replace("X", "");
                            #     # For the final peptide seq, just replace the X's

                            #     CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);
                            #     CORE.PWS("# " + CORE.getDateTime() + " Translated CDS sequence corrected by removing X's from seqs " + " " + cluster_id + " " + spec + " " + cds_header + " " + pid, logfile);
                            #     CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);
                            #     match_code = [seq_type, "match.after.Xs"];
                            #     match_found = True;
                            #     # Log info
                            # ## Check if the seqs match after removing all X's

                            ####################

                            # if not match_found and cur_seqs[seq_type]['translated'].replace("X", "") in pep_seq.replace("X", ""):
                            #     translated_seq = cur_seqs[seq_type]['translated'].replace("X", "");
                            #     # Replace the X's in the translated sequence

                            #     final_codons = [];
                            #     for i in range(len(cur_seqs[seq_type]['translated'])):
                            #         if cur_seqs[seq_type]['translated'][i] != "X":
                            #             final_codons.append(cur_seqs[seq_type]['codons'][i]);
                            #     # Remove the codons that result in X's from the codon sequence                                

                            #     final_cds = "".join(final_codons);
                            #     # For the final CDS, we find the start and end coords of the translated seq in the peptide seq

                            #     pep_start = pep_seq.index(translated_seq);
                            #     pep_end = len(translated_seq) + pep_start;
                            #     final_pep = pep_seq.replace("X", "")[pep_start:pep_end];
                            #     # For the final peptide seq, just replace the X's

                            #     CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);
                            #     CORE.PWS("# " + CORE.getDateTime() + " Translated CDS sequence within peptide seq AFTER removing X's: " + " " + cluster_id + " " + spec + " " + cds_header + " " + pid, logfile);
                            #     CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);
                            #     match_code = [seq_type, "trans.substr.pep.after.Xs"];
                            #     match_found = True;
                            #     # Log info
                            # ## Check if the translated seq is a substring of the petide seq after removing all X's

                            # ####################

                            # if not match_found and pep_seq.replace("X", "") in cur_seqs[seq_type]['translated'].replace("X", ""):
                            #     translated_seq = cur_seqs[seq_type]['translated'].replace("X", "");
                            #     # Remove all the X's from the translated sequence

                            #     final_codons = [];
                            #     for i in range(len(cur_seqs[seq_type]['translated'])):
                            #         if cur_seqs[seq_type]['translated'][i] != "X":
                            #             final_codons.append(cur_seqs[seq_type]['codons'][i]);
                            #     final_cds = "".join(final_codons);
                            #     # Remove all the codons that result in X's and join them as the final CDS sequence

                            #     final_pep = pep_seq.replace("X", "");
                            #     # Remove all X's from the peptide sequence

                            #     pep_start = translated_seq.index(final_pep);
                            #     pep_end = len(final_pep) + pep_start;
                            #     final_pep = final_pep[pep_start:pep_end];
                            #     # For the final peptide seq, we find the start and end coords of the peptdie seq in the translated seq after removing X's from both

                            #     CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);
                            #     CORE.PWS("# " + CORE.getDateTime() + " Peptide sequence within translated seq AFTER removing X's: " + " " + cluster_id + " " + spec + " " + cds_header + " " + pid, logfile);
                            #     CORE.PWS("# " + CORE.getDateTime() + " ===================================", logfile);
                            #     match_code = [seq_type, "pep.substr.trans.after.Xs"];  
                            #     match_found = True;
                            #     # Log info
                            # ## Check if the peptide sequnce is a substring of the translated CDS seq after removing all X's