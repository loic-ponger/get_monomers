#!/usr/bin/env python3


import os
import time
import argparse
import get_monomers_lib_bed as gm
import Bio


##################################################################################	
###################################################################################					
###################################################################################					
###################################################################################					
###################################################################################					
	
#~ def are_hits_overlapping(hit1, hit2):
	#~ if hit1['
###################################################################################					
parser = argparse.ArgumentParser(description='Search for tandem repeats ...', epilog='')
parser.add_argument('-a', dest="action", required=True, choices=['search', 'extract', 'search_and_extract'])
## for search action
parser.add_argument('-s', dest="seq_file",   required=True, default=None)
parser.add_argument('-r', dest="ref_file",   required=False, default=None)
# parser.add_argument('-out_tab_file',                       dest="out_tab_file",               required=False, default=None)
# parser.add_argument('-out_block_seq_index',                dest="out_block_seq_index_file",   required=False, default=None)
# parser.add_argument('-out_block_seq_index_after_merge',    dest="out_block_seq_index_file_after_merge",   required=False, default=None)
# parser.add_argument('-out_block_seq',                      dest="out_block_seq_file",         required=False, default=None)
parser.add_argument('-m',                                  dest="minimal_score",   type=int, default=0)
parser.add_argument('-hmmer1',      dest="hmmer_options1",             default="")
parser.add_argument('-hmmer2',      dest="hmmer_options2",             default="")
parser.add_argument('-extrapol_length',         dest="extrapol_length",         type=int, default=5)

## for extract action
parser.add_argument('-i',                       dest="index_file", required=False, default=None)
parser.add_argument('-p',                       dest="phase",           type=int, default=None)
parser.add_argument('-x',                       dest="suffix",   required=True,default="")

#parser.add_argument('-out_monomer_seq_index',   dest="out_monomer_seq_index_file", required=False, default=None)
#parser.add_argument('-out_monomer_seq',         dest="out_monomer_seq_file",       required=False, default=None)
#parser.add_argument('-out_monomer_gtf',         dest="out_monomer_gtf_file",       required=False, default=None)
## for all actions
parser.add_argument('-v',  dest="verbose",         type=int, default=0)
#parser.add_argument('-k',  dest="keep_temporary_files",         type=boolean, default=0)

##
args = parser.parse_args()
args = gm.parse_options(args)

###################################################################################						
#################################################################
if args.action == 'search' or args.action == 'search_and_extract':
	if args.verbose > 0:
		print("Getting monomer length ...")
		t1=time.time()
	monomer_length = gm.get_monomer_ali_length(args.ref_file)
	if args.verbose > 2:
		print("     monomer length %d (%f s)" % (monomer_length, time.time()-t1))

#################################################################
if args.action == 'search' or args.action == 'search_and_extract':
	if args.verbose > 0:
		print("Formating the hmmer db ...")
		t1=time.time()
	monomer_hmm_file = gm.format_hmmdb(args.ref_file)
	if args.verbose > 2:
		print("     hmm db in %s (%f s)" % (monomer_hmm_file, time.time()-t1))

#################################################################
if args.action == 'search' or args.action == 'search_and_extract':
	if args.verbose > 0:
		print("Doubling the monomer sequence ...")
		t1=time.time()
	double_ref_file = gm.get_double_ref_ali_file(args.ref_file)
	if args.verbose > 2:
		print("     double monomers in %s (%f s)" % (double_ref_file, time.time()-t1))

#################################################################
if args.action == 'search' or args.action == 'search_and_extract':
	if args.verbose > 0:
		print("Formating the hmmer db (double seq) ...")
		t1=time.time()
	double_monomer_hmm_file = gm.format_hmmdb(double_ref_file)
	if args.verbose > 2:
		print("     hmm db in %s (%f s)" % (double_monomer_hmm_file, time.time()-t1))

#################################################################
if args.action == 'search' or args.action == 'search_and_extract':
	if args.verbose > 0:
		print("Searching the sequences ...")
		t1=time.time()
	hmmer_ouput_file = gm.search_hmm(args.seq_file, monomer_hmm_file, options=args.hmmer_options1)
	if args.verbose > 2:
		print("     nhmmscan output in %s (%f s)" % (hmmer_ouput_file, time.time()-t1))

#################################################################
#if args.verbose > 0:
	#print("Searching the sequences (double seq)...")
	#t1=time.time()
#double_hmmer_ouput_file = gm.search_hmm(args.seq_file, double_monomer_hmm_file, options=args.hmmer_options2)
#if args.verbose > 2:
	#print("     nhmmscan output in %s (%f s)" % (double_hmmer_ouput_file, time.time()-t1))

#################################################################
if args.action == 'search' or args.action == 'search_and_extract':
	if args.verbose > 0:
		print("Parsing the hmmer output ...")
		t1=time.time()
	hmm_hits = gm.parse_hmm_ouput(hmmer_ouput_file, args.minimal_score)
	if args.verbose > 2:
		print("     number of sequences: %d (%f s)" % (len(hmm_hits), time.time()-t1))

#################################################################	
if args.action == 'search' or args.action == 'search_and_extract':
	# if args.out_tab_file != None:
		# if args.verbose > 0:
			# print("Writing monomer output file ...")
			# t1=time.time()
		# gm.write_monomer_info(hmm_hits, args.out_tab_file)
		# if args.verbose > 2:
			# print("     output in %s (%f s)" % (args.out_tab_file, time.time()-t1))
	hmm_hits_copy = dict(hmm_hits)

#################################################################	
if args.action == 'search' or args.action == 'search_and_extract':
	if args.out_block_seq_index_file != None:
		if args.verbose > 2:
			print("Writing block sequence/index output file ...")
			t1=time.time()
		gm.write_block_seq_index(hmm_hits, args.seq_file, args.out_block_seq_index_file)
		if args.verbose > 2:
			print("     output in %s (%f s)" % (args.out_block_seq_index_file, time.time()-t1))

#################################################################	
if args.action == 'search' or args.action == 'search_and_extract':
	if args.out_block_seq_file != None:
		if args.verbose > 0:
			print("Writing block sequence output file ...")
			t1=time.time()
		gm.write_block_seq(hmm_hits, args.seq_file, args.out_block_seq_file)
		if args.verbose > 2:
			print("     output in %s (%f s)" % (args.out_block_seq_file, time.time()-t1))

#################################################################	
if args.action == 'search' or args.action == 'search_and_extract':
	if args.verbose > 0:
		print("Extrapolating 5' and 3' index ...")
		t1=time.time()
	hmm_hits = gm.extrapol_index(hmm_hits, monomer_length, args.extrapol_length)
	if args.verbose > 2:
		print("     number of sequences: %d (%f s)" % (len(hmm_hits), time.time()-t1))
		
#################################################################	
nb_junctions = 0
if args.action == 'search' or args.action == 'search_and_extract':
	if args.verbose > 0:
		print("Exporting monomer junctions ...")
		t1=time.time()
	junction_seq_file, nb_junctions = gm.export_junctions(hmm_hits, args.seq_file)
	if args.verbose > 2:
		print("     output in %s (%d junctions, %f s)" % (junction_seq_file, nb_junctions, time.time()-t1))

	
#################################################################
if args.action == 'search' or args.action == 'search_and_extract':
	if nb_junctions > 0:
		if args.verbose > 0:
			print("Searching the sequences (double seq / junctions)...")
			t1=time.time()
		junction_hmmer_ouput_file = gm.search_hmm(junction_seq_file, double_monomer_hmm_file)
		if args.verbose > 2:
			print("     nhmmscan output in %s (%f s)" % (junction_hmmer_ouput_file, time.time()-t1))
		
	    #################################################################
		if args.verbose > 0:
			print("Parsing the hmmer output for junctions...")
			t1=time.time()
		junction_hmm_hits = gm.parse_hmm_ouput(junction_hmmer_ouput_file)
		if args.verbose > 2:
			print("     number of hits: %d (%f s)" % (len(junction_hmm_hits), time.time()-t1))
			
		#################################################################
		if args.verbose > 0:
			print("Merging initial and junction index ...")
			t1=time.time()
		checked_hmm_hits = gm.merging_index(hmm_hits, junction_hmm_hits, monomer_length, args.verbose)
		if args.verbose > 2:
			print("     number of hits: %d (%f s)" % (len(checked_hmm_hits), time.time()-t1))
		#################################################################	
		if args.verbose > 0:
			print("Writing block sequence/index output file (after merging) ...")
			t1=time.time()
		os.rename(args.out_block_seq_index_file, args.out_block_seq_index_file_before_merge)
		if args.verbose > 2:
			print("     output in %s backup to %s" % (args.out_block_seq_index_file, args.out_block_seq_index_file_before_merge))
		gm.write_block_seq_index(checked_hmm_hits, args.seq_file, args.out_block_seq_index_file)
		if args.verbose > 2:
			print("     output in %s (%f s)" % (args.out_block_seq_index_file, time.time()-t1))
		#################################################################	
		if args.verbose > 0:
			print("Checking stuffsss [TO DO: validate this step]...")
			t1=time.time()
			junction_seq_file_before_after, nb_junctions_before_after = gm.export_junctions_before_after(hmm_hits_copy, checked_hmm_hits, args.seq_file)
		if args.verbose > 2:
			print("     output in %s (%d junctions, %f s)" % (junction_seq_file_before_after, nb_junctions_before_after, time.time()-t1))
	else:
		if args.verbose > 0:
			print("     No ambigous junctions ...")
#################################################################	






if args.action == 'search_and_extract':
	if args.verbose > 0:
		t1=time.time()
		print("Reading sequence/index in %s ..."% (args.out_block_seq_index_file))
	seq_index = gm.read_sequence_index_data(args.out_block_seq_index_file)
	if args.verbose > 2:
		print("     %d sequence/index (%f s)" % (len(seq_index), time.time()-t1))
		
#################################################################	
if args.action == 'extract':
	if args.verbose > 0:
		print("Reading sequence/index ...")
		t1=time.time()
	seq_index = gm.read_sequence_index_data(args.out_block_seq_index_file)
	if args.verbose > 2:
		print("     %d sequence/index (%f s)" % (len(seq_index), time.time()-t1))
		
#################################################################	
if args.action == 'search_and_extract' or args.action == 'extract':
	if args.verbose > 0:
		print("Spliting monomers (using phase %d) ..." % (args.phase))
		t1=time.time()
	monomers_seq_index, monomers_descriptions = gm.split_monomers(seq_index, args.phase)
	if args.verbose > 2:
		print("     %d monomers (%f s)" % (len(monomers_seq_index), time.time()-t1))
		
#################################################################	
if args.action == 'search_and_extract' or args.action == 'extract':
	if args.out_monomer_seq_file != None:	
		if args.verbose > 0:
			print("Writing monomer sequences ...")
			t1=time.time()
		gm.writing_monomers_sequence(monomers_seq_index, args.out_monomer_seq_file)
		if args.verbose > 2:
			print("     %d monomers (%f s)" % (len(monomers_seq_index), time.time()-t1))
			print("     output in %s" % (args.out_monomer_seq_file))
#################################################################	
if args.action == 'search_and_extract' or args.action == 'extract':
	if args.out_monomer_gtf_file != None:	
		if args.verbose > 0:
			print("Writing monomer description (GTF format) ...")
			t1=time.time()
		gm.writing_monomers_descriptions(monomers_descriptions, args.out_monomer_gtf_file)
		if args.verbose > 2:
			print("     %d descriptions (%f s)" % (len(monomers_descriptions), time.time()-t1))
			print("     output in %s" % (args.out_monomer_gtf_file))

#################################################################
if args.action == 'search_and_extract' or args.action == 'extract':
	if args.out_monomer_seq_index_file != None:	
		if args.verbose > 0:
			print("Writing monomers sequences/index ...")
			t1=time.time()
		gm.writing_monomers_sequence_index(monomers_seq_index, args.out_monomer_seq_index_file)
		if args.verbose > 2:
			print("     %d monomers (%f s)" % (len(monomers_seq_index), time.time()-t1))
			print("     output in %s" % (args.out_monomer_seq_index_file))
		if args.verbose > 2:
			if nb_junctions > 0:
				print("     nhmmscan output in %s (%f s)" % (junction_hmmer_ouput_file, time.time()-t1))
				print("     output in %s (%d junctions, %f s)" % (junction_seq_file, nb_junctions, time.time()-t1))
			else:
				print("     no junction files (%f s)" % (time.time()-t1))
#################################################################
if args.verbose > 0:
	print("Bye bye!!!")
quit()



