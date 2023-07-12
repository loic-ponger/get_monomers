import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import tempfile
from itertools import compress
import time
import shutil as sh
from operator import itemgetter
#from werkzeug import secure_filename



######################################
### to get sequence length
### first seq is used for multifasta file 
######################################
def get_monomer_length(monomer_file):
	monomer_seq = list(SeqIO.parse(monomer_file, "fasta"))[0]
	return(len(monomer_seq.seq))

######################################
### to get sequence length
### check that all sequences have the same length
######################################

def get_monomer_ali_length(monomer_file):
	all_seqs = list(SeqIO.parse(monomer_file, "fasta"))
	length = len(all_seqs[0].seq)
	for ss in all_seqs:
		if len(ss.seq) != length:
			print("Error: sequences in %s exhibit different lengths (%d vs. %d). Sequences should be aligned." % (monomer_file, length, len(ss.seq)))
			print("Bye bye !!")
			quit()
	return(length)
	
######################################
### write genomic sequences of monomers
######################################
	
def writing_monomers_sequence(monomer_sequence_index, output_file):
	f = open(output_file, 'w')
	for name in monomer_sequence_index:
			f.write(">%s\n" % (name))
			f.write("%s\n" % (monomer_sequence_index[name]['sequence']))
	f.close()		
	return(0)
######################################
### write descriptions of monomers (in GTF format)
######################################
	
def writing_monomers_descriptions(monomer_descriptions, output_file):
	f = open(output_file, 'w')
	for name in monomer_descriptions:
			f.write("%s\t" % (monomer_descriptions[name]['sequence_name']))
			f.write("%s\t" % ("monosearch"))
			f.write("%s\t" % ("monomer"))
			f.write("%d\t" % (monomer_descriptions[name]['sequence_begin']))
			f.write("%d\t" % (monomer_descriptions[name]['sequence_end']))
			f.write("%d\t" % (monomer_descriptions[name]['length']))
			f.write("%s\t" % (monomer_descriptions[name]['strand']))
			f.write("%d\t" % (monomer_descriptions[name]['phase']))
			f.write("%s=%d;" % ("block_begin", monomer_descriptions[name]['block_begin']))
			f.write("%s=%d;" % ("block_end", monomer_descriptions[name]['block_end']))
			f.write("%s=%d;" % ("begin_in_block", monomer_descriptions[name]['begin_in_block']))
			f.write("%s=%d;" % ("end_in_block", monomer_descriptions[name]['end_in_block']))
			f.write("%s=%d;" % ("initial_phase", monomer_descriptions[name]['initial_phase']))
			f.write("%s=%d;" % ("num", monomer_descriptions[name]['nb']))
			f.write("%s=%d;" % ("final_phase", monomer_descriptions[name]['final_phase']))
			f.write("%s=%s;" %   ("monomer_name", name.replace(">","")))
			f.write("%s=%s;" % ("block_name", monomer_descriptions[name]['block_name']))

			
			
			f.write("\n")
	f.close()		
	return(0)######################################
### write genomic sequences + indexes of monomers
######################################

def writing_monomers_sequence_index(monomer_sequence_index, output_file):
	f = open(output_file, 'w')
	for name in monomer_sequence_index:
			f.write(">%s\n" % (name))
			f.write("%s\n" % (monomer_sequence_index[name]['sequence']))
			f.write(">%s\n" % (name))
			f.write("%s\n" % (monomer_sequence_index[name]['index']))
	f.close()		
	return(0)

######################################
### split blocks into monomers
### according to the given phase
######################################
def split_monomers(sequence_index, phase=1):
	monomers = {}
	descriptions = {}
	for name in sequence_index:
		ii = 0
		nb=0
		ll = len(sequence_index[name]['sequence'])
		# to escape the first 0s and -1s
		while ii < ll and sequence_index[name]['index'][ii] < 1:
			ii = ii + 1				
		if ii < ll:
			begin = ii
			initial_phase = sequence_index[name]['index'][ii]
			end = ii
			final_phase = sequence_index[name]['index'][ii]
			prev_index = ii
			ii = ii + 1
			while ii < ll:
				if sequence_index[name]['index'][ii] > 0:
					if sequence_index[name]['index'][ii] == phase or (sequence_index[name]['index'][ii] > phase and sequence_index[name]['index'][prev_index] < phase) or (sequence_index[name]['index'][prev_index] > sequence_index[name]['index'][ii] > phase and sequence_index[name]['index'][ii] > phase):
						monomer_name = name + "_" + str(begin) + "_" + str(end)
						monomers[monomer_name] = {'sequence': None, 'index': None}
						monomers[monomer_name]['sequence'] = sequence_index[name]['sequence'][begin:(end+1)]
						monomers[monomer_name]['index']    = sequence_index[name]['index'][begin:(end+1)]
						descriptions[monomer_name] = {}
						tmp=name.split("_")
						descriptions[monomer_name]['sequence_name']='_'.join(tmp[0:-5]).replace(">","")
						descriptions[monomer_name]['block_name']=name.replace(">","")
						descriptions[monomer_name]['nb']=nb
						descriptions[monomer_name]['strand']=tmp[-1]
						descriptions[monomer_name]['block_begin']=int(tmp[-3])
						descriptions[monomer_name]['block_end']=int(tmp[-2])
						descriptions[monomer_name]['begin_in_block']=begin
						descriptions[monomer_name]['end_in_block']=end
						if descriptions[monomer_name]['strand'] == 'plus':
							descriptions[monomer_name]['sequence_begin']=descriptions[monomer_name]['block_begin']+descriptions[monomer_name]['begin_in_block']
							descriptions[monomer_name]['sequence_end']=descriptions[monomer_name]['block_begin']+descriptions[monomer_name]['end_in_block']
						else:
							descriptions[monomer_name]['sequence_begin']=descriptions[monomer_name]['block_end']-descriptions[monomer_name]['begin_in_block']+1
							descriptions[monomer_name]['sequence_end']=descriptions[monomer_name]['block_end']-descriptions[monomer_name]['end_in_block']+1
						descriptions[monomer_name]['length']=end-begin+1
						descriptions[monomer_name]['initial_phase']=initial_phase
						descriptions[monomer_name]['final_phase']=final_phase
						descriptions[monomer_name]['phase']=phase
						begin = ii
						nb = nb + 1
						initial_phase = sequence_index[name]['index'][ii]
					else:
						end = ii
						final_phase = sequence_index[name]['index'][ii]			
						
					prev_index = ii
				ii = ii + 1				
			monomer_name = name + "_" + str(begin) + "_" + str(end)
			monomers[monomer_name] = {'sequence': None, 'index': None}
			monomers[monomer_name]['sequence'] = sequence_index[name]['sequence'][begin:(end+1)]
			monomers[monomer_name]['index']    = sequence_index[name]['index'][begin:(end+1)]
			descriptions[monomer_name] = {}
			tmp=name.split("_")
			descriptions[monomer_name]['sequence_name']='_'.join(tmp[0:-5]).replace(">","")
			descriptions[monomer_name]['block_name']=name.replace(">","")
			descriptions[monomer_name]['nb']=nb
			descriptions[monomer_name]['strand']=tmp[-1]
			descriptions[monomer_name]['block_begin']=int(tmp[-3])
			descriptions[monomer_name]['block_end']=int(tmp[-2])
			descriptions[monomer_name]['begin_in_block']=begin
			descriptions[monomer_name]['end_in_block']=end
			if descriptions[monomer_name]['strand'] == 'plus':
				descriptions[monomer_name]['sequence_begin']=descriptions[monomer_name]['block_begin']+descriptions[monomer_name]['begin_in_block']
				descriptions[monomer_name]['sequence_end']=descriptions[monomer_name]['block_begin']+descriptions[monomer_name]['end_in_block']
			else:
				descriptions[monomer_name]['sequence_begin']=descriptions[monomer_name]['block_end']-descriptions[monomer_name]['begin_in_block']+1
				descriptions[monomer_name]['sequence_end']=descriptions[monomer_name]['block_end']-descriptions[monomer_name]['end_in_block']+1
			descriptions[monomer_name]['length']=end-begin+1
			descriptions[monomer_name]['initial_phase']=initial_phase
			descriptions[monomer_name]['final_phase']=final_phase
			descriptions[monomer_name]['phase']=phase

		
	return(monomers, descriptions)
	

######################################
### read DNA sequence and index from pseudo fasta file
### >seq1
### ACCGGTGCGT
### >seq1
### [0, 0, 0, ...0]
######################################
def read_sequence_index_data(input_file):
	rc={}
	titleline = re.compile(r"^>.*$")
	mode = None
	with open(input_file, "rU") as f:
		for line in f:
			if re.search("^>.*_seq$", line):
				mode ='sequence'				
				nm=line.split()[0]
				nm = nm[:-4]
				if nm not in rc:
					rc[nm]={'sequence': None, 'index': None}
				if rc[nm]['sequence'] != None:
					print("Error in sequence/index file (%s)." %(input_file))
					print("Sequence should be unique")
					print("Not the case for %s" %(nm))
					quit(-1)
					
			elif re.search("^>.*_index$", line):
				mode ='index'
				nm=line.split()[0]
				nm = nm[:-6]
				if nm not in rc:
					rc[nm]={'sequence': None, 'index': None}
				if rc[nm]['index'] != None:
					print("Error in sequence/index file (%s)." %(input_file))
					print("Index should be unique")
					print("Not the case for %s" %(nm))
					quit(-1)
			elif not re.search("^>.*$", line):
				if re.search("^[a-zA-Z]*$", line) and mode == 'sequence':
					rc[nm]['sequence']=line.rstrip()
				elif mode == 'index':
					rc[nm]['index']=list(map (int, line.split(", ")[1:]))	
				else:
					print("Error in sequence/index file (%s)." %(input_file))
					print("Problem of format")
					quit(-1)
	check_sequence_index(rc, stop_on_error=True)
	return(rc)
	
	
######################################
### check if index/sequences are okay
######################################
def check_sequence_index(sequence_index, stop_on_error=True):
	nb_error = 0
	for kk in sequence_index:
		if sequence_index[kk]['sequence'] == None:
			print("Error: %s as no DNA sequence" %(kk))
			nb_error = nb_error + 1
			if stop_on_error == True:
				quit(-2)
		if sequence_index[kk]['index'] == None:
			print("Error: %s as no index" %(kk))
			nb_error = nb_error + 1
			if stop_on_error == True:
				quit(-2)
		lseq = len(sequence_index[kk]['sequence'])
		lidx = len(sequence_index[kk]['index'])
		if sequence_index[kk]['sequence'] != None and sequence_index[kk]['index'] != None and lidx!= lseq:
			print("Error: %s sequence/index with different length (%d bp, %d integers)" %(kk, lseq, lidx))
			nb_error = nb_error + 1
			print(sequence_index[kk]['sequence'])
			print(sequence_index[kk]['sequence'][-2:])
			print(sequence_index[kk]['index'])
			if stop_on_error == True:
				quit(-2)
	return(nb_error)

##############################################################
### merge index of initial search with index of search
### using a double monomer against the ambigous junctions
###############################################################


def merging_index(inital_hmmer_hits, junction_hmmer_hits, length_monomer, verbose=1):
	for name in junction_hmmer_hits:
		#print(name)
		string=name.split("_")[-1]
		end=int(name.split("_")[-2])
		begin=int(name.split("_")[-3])
		block_end=int(name.split("_")[-4])
		block_begin=int(name.split("_")[-5])
		original_name=re.sub(r"_[^_]*_[^_]*_[^_]*_[^_]*_[^_]*$",'', name)
		# print(original_name)
		# print(string)
		# print("----")
		# print(begin)
		# print(end)
		# print("----")		
		# print(block_begin)
		# print(block_end)
		# print("----")
		if string == 'plus':
			gap=begin-block_begin
		else:
			gap=block_end-end
			#to check
		# print("str: %s gap: %d" %(string, gap))
		
		# print(junction_hmmer_hits[name])
		fd=-1
		ii=0
		while fd == -1 and ii < len(inital_hmmer_hits[original_name]):
			# print("______________")
			# print(fd)
			# print(ii)
			# print(inital_hmmer_hits[original_name][ii]['begin'])
			# print(block_begin)
			# print(inital_hmmer_hits[original_name][ii]['end'])
			# print(block_end)
			if inital_hmmer_hits[original_name][ii]['begin'] == block_begin and inital_hmmer_hits[original_name][ii]['end'] == block_end:
				fd=ii
			ii=ii+1
		if fd == -1:
			print("ERROR: check the junction file. Sequence %s not found." %(name))
			quit()
		if len(junction_hmmer_hits[name]) > 1:
			print("ERROR: check the junction file. More than 1 junction for %s  (%d junctions)." %(name, len(junction_hmmer_hits[name])))
			continue
			#quit()
			
			
		for ii in range(1,len(junction_hmmer_hits[name][0]['index'])):			
			if junction_hmmer_hits[name][0]['index'][ii] > length_monomer:
				junction_hmmer_hits[name][0]['index'][ii]=junction_hmmer_hits[name][0]['index'][ii]-length_monomer
			
		offset = get_offset(inital_hmmer_hits[original_name][fd]['index'], junction_hmmer_hits[name][0]['index'], gap) 	
		
		for ii in range(1,len(junction_hmmer_hits[name][0]['index'])):
			#inital_hmmer_hits[original_name][fd]['index'][begin+ii-1]
			#print(begin)
			idx = junction_hmmer_hits[name][0]['index'][ii]
			#ind = inital_hmmer_hits[original_name][fd]['index'][ii+gap+offset]
			#print("%s -> ii: %d idx: %d  ||| sq: %s bl: %d blbg: %d idx: %d %d" %(name, ii, idx, original_name, fd, block_begin, inital_hmmer_hits[original_name][fd]['index'][ii+gap+offset]  , offset      ))
			#
			#
			# to check
			#
			
			if idx != -1 and idx != 0:
				#print("##############")
				#print("seq: %s     fd: %d     %d+%d+%d  -> %d  / %d" %   (original_name,fd,ii,gap,offset, idx, inital_hmmer_hits[original_name][fd]['index'][ii+gap+offset]))
				#print (len(inital_hmmer_hits[original_name][fd]['index']))
				inital_hmmer_hits[original_name][fd]['index'][ii+gap+offset] = idx 
				#print("okkkk")
		#print(inital_hmmer_hits[original_name][fd]['index'])
		#print(junction_hmmer_hits[name][0]['index'])
	return(inital_hmmer_hits)
##############################################################
### 
###
###############################################################
def get_offset(index1, index2, gap, offset_max=5):
	offset=-offset_max
	ll=len(index2)
	lll=len(index1)
	#print("offset:----------------------")
	mintmp=10000000000000000000000
	minnbtmp=10000000000000000000000
	for offs in range(-offset_max, offset_max+1):
		tmp=0
		tmpnb=0
		for ii in range(1, ll):
			if gap+ii+offs < lll and index1[gap+ii+offs] > 0 and index2[ii] > 0:
				if index1[gap+ii+offs] != index2[ii]:
					tmpnb=tmpnb+1
				tmp=tmp+abs(index1[gap+ii+offs]-index2[ii])
		if tmpnb < minnbtmp:
			minnbtmp = tmpnb
			offset = offs
		#print("offset: %d  %d  %d  %d" % (offs, tmpnb, minnbtmp, gap))
	return(offset)

##############################################################
### merge index of initial search with index of search
### using a double monomer against the ambigous junctions
###############################################################
	
def merging_index_old(inital_hmmer_hits, junction_hmmer_hits, length_monomer, verbose=1):
	for target_name in inital_hmmer_hits.keys():
		for hits in inital_hmmer_hits[target_name]:
			minus_previous_is_gap_or_overlap = False
			for hit in hits['hits']:	
				#if hits['hits'][0]['strand'] == 'plus' and (hit['junction'] == 'gap' or hit['junction'] == 'overlap' or hit['junction'] == 'perfect_ali-env'):
				if hits['hits'][0]['strand'] == 'plus' and (hit['junction'] == 'gap' or hit['junction'] == 'overlap'):
					#print("############")
					mkeyyyy = '%s_%d_%s' % (target_name, previous_alifrom, hits['hits'][0]['strand'])
					if mkeyyyy not in junction_hmmer_hits.keys():
						if verbose > -1:
							print ("WARNING: The key (sequence name) %s not found in junction hmmer output. Check the junction fasta file and the hmmer output file" % (mkeyyyy))
							quit()
					else:
						for junction_hits in junction_hmmer_hits[mkeyyyy]:
							#print(junction_hits)
							begin = previous_alifrom + junction_hits['begin'] - 1 - hits['begin']
							for ii in range(1,len(junction_hits['index'])):
								ind = junction_hits['index'][ii] 
								if ind > length_monomer:
									#print("ind %d mono %d" %(ind, length_monomer))
									ind = ind - length_monomer
								if hits['index'][begin+ii] != ind:
									if hits['index'][begin+ii] < 1 and (ind > 0 or ind == -2): 
										if verbose > 9:
											print ("NORMAL: replacement %d  by %d at position %d" % (hits['index'][begin+ii], ind, ii) )
										hits['index'][begin+ii] = ind
									else:
										if verbose > 9:
											print ("WARNING: replacement %d  by %d at position %d" % (hits['index'][begin+ii], ind, ii) )
				
				print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")				
				print(hit)

				if hits['hits'][0]['strand'] == 'minus' and minus_previous_is_gap_or_overlap == True:
					mkeyyyy = '%s_%d_%s' % (target_name, previous_alifrom, hits['hits'][0]['strand'])
					if mkeyyyy not in junction_hmmer_hits.keys():
						if verbose > -1:
							print(hits)
							print(junction_hmmer_hits.keys())
							print ("ERROR MINUS: The key (sequence name) %s not found in junction hmmer output. Check the junction fasta file and the hmmer output file" % (mkeyyyy))
							quit()
					else:
						
						for junction_hits in junction_hmmer_hits[mkeyyyy]:
							print("======##################################============")
							print(hit)
							print(junction_hits)
							begin = previous_alifrom + junction_hits['begin'] - 1 - hits['begin']
							print("%d %d %d" % (previous_alifrom, junction_hits['begin'], hits['begin']))
							for ii in range(1,len(junction_hits['index'])):
								ind = junction_hits['index'][ii] 
								print("=0 ind: %d %d"  %(ind, length_monomer))
								if ind > length_monomer:
									ind = ind - length_monomer
									print("===> ind: %d %d"  %(ind, length_monomer))
								if hits['index'][begin+ii] != ind:
									if hits['index'][begin+ii] < 1 and (ind > 0 or ind == -2): 
										if verbose > 9:
											print ("NORMAL MINUSSSSSSSSSSSSSSSS: replacement %d  by %d at position %d" % (hits['index'][begin+ii], ind, ii) )
										print("begin: %d, begin+ii: %d, ii: %d, ind: %d old: %d" %(begin,begin+ii,ii, ind, hits['index'][begin+ii]))
										hits['index'][begin+ii] = ind
										print("##################################")
										print(junction_hits)
										print("##################################")
										#else:
										#if verbose > 9:
											#print ("WARNING MINUS: replacement %d  by %d at position %d" % (hits['index'][begin+ii], ind, ii) )



				minus_previous_is_gap_or_overlap = False		
				#if hits['hits'][0]['strand'] == 'minus' and (hit['junction'] == 'gap' or hit['junction'] == 'overlap' or hit['junction'] == 'perfect_ali-env'):
				if hits['hits'][0]['strand'] == 'minus' and (hit['junction'] == 'gap' or hit['junction'] == 'overlap' ):
					minus_previous_is_gap_or_overlap = True
				previous_alifrom = hit['envfrom']

	return(inital_hmmer_hits)




##############################################################
### doubling monomer sequence
###############################################################
def get_double_ref_file(monomer_file):
	double_monomer_file = tempfile.NamedTemporaryFile().name
	monomer_seq = list(SeqIO.parse(monomer_file, "fasta"))[0]
	monomer_seq.seq=monomer_seq.seq + monomer_seq.seq
	SeqIO.write(monomer_seq, double_monomer_file, "fasta") 
	return(double_monomer_file)
	
##############################################################
### doubling monomer sequence in alignement
###############################################################
def get_double_ref_ali_file(monomer_file):
	double_monomer_file = tempfile.NamedTemporaryFile().name
	all_double_seqs = []
	for ss in list(SeqIO.parse(monomer_file, "fasta")):
		monomer_seq = ss
		monomer_seq.seq = monomer_seq.seq + monomer_seq.seq
		all_double_seqs.append(monomer_seq)
		
	SeqIO.write(all_double_seqs, double_monomer_file, "fasta") 
	return(double_monomer_file)

##############################################################
### 
###############################################################

def write_monomer_info(hmm_hits, outfile):
	f = open(outfile, 'w')
	f.write('sequence length hits_nb seq_begin seq_end -> hit_junction hit_hmm_from hit_hmm_to hit_ali_from hit_alito hit_env_from hit_env_to hit_strand hit_evalue hit_score monomer_seq target_seq\n')
	for target_name in hmm_hits.keys():
		alito=1
		envto=1
		for hits in hmm_hits[target_name]:
			for hit in hits['hits']:	
	#			print("%s %d  %d   %d -- %d  %d -- %d (%s) %f" % (target_name, hit['alifrom']-alito-1,  hit['envfrom']-envto-1, hit['alifrom'], hit['alito'], hit['envfrom'], hit['envto'], hit["strand"]))
				f.write("%s %d %d %d - %d -> %s %d - %d %d - %d  %d - %d %s %f %f %s %s\n" % (target_name,  hits['length'], hits['nb_hits'], hits['begin'], hits['end'],hit['junction'], hit['hmmfrom'], hit['hmmto'], hit['alifrom'], hit['alito'], hit['envfrom'], hit['envto'], hit["strand"], hit['evalue'], hit['score'], hit['monomer_seq'], hit['target_seq']))
				alito=hit['alito']
				envto=hit['envto']	
			#f.write("%s\n" % (hits['index']))
	f.close()
	
def write_monomer_gtf(descriptions, outfile):
	f = open(outfile, 'w')
	f.write('sequence length hits_nb seq_begin seq_end -> hit_junction hit_hmm_from hit_hmm_to hit_ali_from hit_alito hit_env_from hit_env_to hit_strand hit_evalue hit_score monomer_seq target_seq\n')
	for target_name in hmm_hits.keys():
		alito=1
		envto=1
		for hits in hmm_hits[target_name]:
			for hit in hits['hits']:	
	#			print("%s %d  %d   %d -- %d  %d -- %d (%s) %f" % (target_name, hit['alifrom']-alito-1,  hit['envfrom']-envto-1, hit['alifrom'], hit['alito'], hit['envfrom'], hit['envto'], hit["strand"]))
				f.write("%s %d %d %d - %d -> %s %d - %d %d - %d  %d - %d %s %f %f %s %s\n" % (target_name,  hits['length'], hits['nb_hits'], hits['begin'], hits['end'],hit['junction'], hit['hmmfrom'], hit['hmmto'], hit['alifrom'], hit['alito'], hit['envfrom'], hit['envto'], hit["strand"], hit['evalue'], hit['score'], hit['monomer_seq'], hit['target_seq']))
				alito=hit['alito']
				envto=hit['envto']	
			#f.write("%s\n" % (hits['index']))
	f.close()
	
	
	
##############################################################
### write block sequences (fasta)
### write block index (pseudo fasta)
###############################################################

def write_block_seq_index_put_all_seq_into_memory(hmm_hits, seq_file, outfile):
	record_dict = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))
	f = open(outfile, 'w')
	all_blocks =  []
	for target_name in hmm_hits.keys():
		for hits in hmm_hits[target_name]:
			#print(record_dict[target_name])
			block_frag = record_dict[target_name].seq[(hits['begin']-1):(hits['end']-1+1)]	
			if hits['hits'][0]['strand'] == 'minus':
				block_frag = SeqRecord(block_frag,'toto', '', '').reverse_complement().seq
	
			record = ('%s_%d_%d_%d_%d_%s' % (target_name, hits['length'], hits['nb_hits'], hits['begin'], hits['end'], hits['hits'][0]['strand']),  block_frag, hits['index'])
			all_blocks.append(record)	
			f.write(">%s_seq\n" % (record[0]))
			f.write("%s\n" % (record[1]))
			f.write(">%s_index\n" % (record[0]))
			a=", "
			f.write("%s\n" % (a.join(list(map (str, record[2])))))		
			#print("###############")
			#print(len(record[2]))
			#print(len(record[1]))
			#print("###############")
	f.close()		

###############################################################
def write_block_seq_index(hmm_hits, seq_file, outfile):
	f = open(outfile, 'w')
	all_blocks =  []
	
	with open(seq_file, "rU") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			target_name = record.id
			if target_name in hmm_hits.keys():
				for hits in hmm_hits[target_name]:
					#print(record_dict[target_name])
					block_frag = record.seq[(hits['begin']-1):(hits['end']-1+1)]	
					if hits['hits'][0]['strand'] == 'minus':
						block_frag = SeqRecord(block_frag,'toto', '', '').reverse_complement().seq
			
					data = ('%s_%d_%d_%d_%d_%s' % (target_name, hits['length'], hits['nb_hits'], hits['begin'], hits['end'], hits['hits'][0]['strand']),  block_frag, hits['index'])
					all_blocks.append(data)	
					f.write(">%s_seq\n" % (data[0]))
					f.write("%s\n" % (data[1]))
					f.write(">%s_index\n" % (data[0]))
					a=", "
					f.write("%s\n" % (a.join(list(map (str, data[2])))))		
			#print("###############")
			#print(len(record[2]))
			#print(len(record[1]))
			#print("###############")
	handle.close()
	f.close()		
##############################################################
### write block sequences only
#################################

def write_block_seq(hmm_hits, seq_file, outfile):
	f = open(outfile, 'w')
	all_blocks =  []
	with open(seq_file, "rU") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			target_name = record.id
			if target_name in hmm_hits.keys():
				for hits in hmm_hits[target_name]:
					#print(record_dict[target_name])
					block_frag = record.seq[(hits['begin']-1):(hits['end']-1+1)]	
					if hits['hits'][0]['strand'] == 'minus':
						block_frag = SeqRecord(block_frag,'toto', '', '').reverse_complement().seq
			
					data = ('%s_%d_%d_%d_%d_%s' % (target_name, hits['length'], hits['nb_hits'], hits['begin'], hits['end'], hits['hits'][0]['strand']),  block_frag, hits['index'])
					all_blocks.append(data)	
					f.write(">%s\n" % (data[0]))
					f.write("%s\n" % (data[1]))
	handle.close()
	f.close()		

def write_block_seq_put_all_seq_into_memory(hmm_hits, seq_file, outfile):
	record_dict = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))
	f = open(outfile, 'w')
	all_blocks =  []
	for target_name in hmm_hits.keys():
		for hits in hmm_hits[target_name]:
			#print(record_dict[target_name])
			block_frag = record_dict[target_name].seq[(hits['begin']-1):(hits['end']-1+1)]	
			if hits['hits'][0]['strand'] == 'minus':
				block_frag = SeqRecord(block_frag,'toto', '', '').reverse_complement().seq
	
			record = ('%s_%d_%d_%d_%d_%s' % (target_name, hits['length'], hits['nb_hits'], hits['begin'], hits['end'], hits['hits'][0]['strand']),  block_frag, hits['index'])
			all_blocks.append(record)	
			f.write(">%s\n" % (record[0]))
			f.write("%s\n" % (record[1]))
	f.close()		





##############################################################
### create temporary directory
###############################################################
def create_temp_folder(directory, verbose=0):
	fullname_dir=tempfile.mkdtemp(dir=directory)
	shortname_dir=os.path.basename(fullname_dir)
	if verbose >  5:
		print("#####")
		print("Directory full name: %s" % fullname_dir)
		print("Directory short name: %s" % shortname_dir)	
	return((shortname_dir, fullname_dir))


##############################################################
### get parameter from form
###############################################################

def get_parameters(request_form):
	parameters={}
	parameters['index'] = int(request_form['index'])
	return(parameters)


##############################################################
### get genome file
###############################################################
def get_genome_file(request_files, directory, verbose=0):
	### parsing file name from the form
	genome_file=request_files['genome_file']
	genome_filename = secure_filename(genome_file.filename)	
	if genome_filename != "":
		genome_file.save(os.path.join(directory, 'genome.fst'))
	else:
		sh.copyfile("static/genome_example2.fst", os.path.join(directory, "genome.fst"))
	return("genome.fst")


##############################################################
### get reference file
###############################################################
def get_reference_file(request_files, directory, verbose=0):
	### parsing file name from the form
	reference_file=request_files['reference_file']
	reference_filename = secure_filename(reference_file.filename)	
	if reference_filename != "":
		reference_file.save(os.path.join(directory, 'reference.fst'))
	else:
		sh.copyfile("static/reference_double_example.fst", os.path.join(directory, "reference.fst"))
	return("reference.fst")

##############################################################
### delete old files
###########################################################
def delete_old_files(path, n_days, verbose=9):	
	nb = 0	
	now = time.time()
	for dd in os.listdir(path):
		if os.stat(os.path.join(path, dd)).st_mtime > (now - n_days * 24 * 60 *60):
			if os.path.isdir(os.path.join(path, dd)):
				if verbose > 8:
					print("     Removing %s ..." % os.path.join(path, dd))
				sh.rmtree(os.path.join(path, dd))
				nb = nb + 1	
	return(nb)

##############################################################
### format hmmer db 
##############################################################
def format_hmmdb(seq_file, db_name=None, hmmpress=True):
	name = "monomer"  # ne pas modifier
	if db_name == None:
		db_name = tempfile.NamedTemporaryFile().name

	# create command line
	command = "hmmbuild --dna -n " + name + " " + db_name + " " + seq_file
	# running the formatting
	result = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)

	# hmmpress if needed 
	if hmmpress == True:
		command = "hmmpress " + db_name
		result = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)
	return(db_name)
	
##############################################################
### search a sequence file with a hmmer db 
##############################################################
def search_hmm(seq_file, monomer_hmm_file, output_file=None, cut_ga=False, options="", prog="nhmmer"):
	if output_file == None:
		output_file = tempfile.NamedTemporaryFile().name
	# create command line
	command = "nhmmer -o " + output_file
	if cut_ga == True:
		command = command + " --cut_ga "
	command = command + " " + options + " "
	command = command + " " + monomer_hmm_file + " " + seq_file
	# running the search
	result = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)
	return(output_file)

##############################################################
### extract regions 
##############################################################
def extract_region(regions, seq_file, output_file=None):
	if output_file == None:
		output_file = tempfile.NamedTemporaryFile().name

	fasta_sequences = SeqIO.parse(open(seq_file),'fasta')
	with open(output_file) as out_file:
		for fasta in fasta_sequences:
			if fasta.id in regions:
				for ii in range(0, len(regions[fasta.id])):
					print (regions[fasta.id][ii]['begin'])
					print (regions[fasta.id][ii]['end'])
					print (regions[fasta.id][ii]['strand'])
					
	return(output_file)




##############################################################
### parse a hmmer output
##############################################################
def parse_hmm_ouput(hmmer_ouput_file, minimal_score=0):
	output = {}
	fffile = open(hmmer_ouput_file,'r')
	hit = None
	for line in fffile:
		if line.startswith("Query:"):
			model_name, model_length = re.split("[ \n]+", line)[1:3]	

		if line.startswith("//"):  # to get the last hit
				if hit != None and hit['score'] >= minimal_score:
					if target_name not in output:
						output[target_name] = []
					output[target_name].append(hit)
				
		if line.startswith(">>"):
				if hit != None and hit['score'] >= minimal_score:
					if target_name not in output:
						output[target_name] = []
					output[target_name].append(hit)
			
				target_name = re.split("[ \n]+", line)[1]
				#if target_name not in output:
					 #output[target_name] = []
				fffile.readline()
				fffile.readline()
				line = fffile.readline()
				hit = get_hit(line)
				hit['monomer_seq']=""
				hit['target_seq']=""
		#if re.match(r' *[^ ]+ +- -----+ +-', line):
		if 	re.match(r' *[^ ]+ +[0-9]+ .+ [0-9]+', line) or re.match(r' *[^ ]+ +- -----+ +-', line): 
				if re.split("[ ]+", line)[1] == "monomer":
					monomer_seq_ali = re.split(" +", line)[3]
					hit['monomer_seq']=hit['monomer_seq'] + monomer_seq_ali
					#print("monomer: %s" % monomer_seq_ali)
				else:
					target_seq_ali  = re.split(" +", line)[3]
					hit['target_seq']=hit['target_seq'] + target_seq_ali
					#print("target: %s" % target_seq_ali)				
	for target_name in 	output.keys():
			output[target_name] = sorted(output[target_name], key=lambda k: k['alifrom'])					
			output[target_name] = remove_overlaping_hits_plusminus(output[target_name])
			output[target_name] = split_blocks(output[target_name], target_name)
			output[target_name] = reverse_all_minus_blocks(output[target_name])
			output[target_name] = junction_all_blocks(output[target_name])
			output[target_name] = index_all_blocks(output[target_name])					
			#output[target_name] = sorted(output[target_name], key=lambda k: k['alifrom'])					
			for hits in output[target_name]:
				hits['hits'] = sorted(hits['hits'], key=itemgetter('alito')) 
			
			#print_all_blocks(output[target_name])
	#print(output)
	#quit()
	return(output)
	
##############################################################
### check and adjust monomer junctions
### all blocks
##############################################################
def junction_all_blocks(blocks):
	for ii in range(0, len(blocks)):
			blocks[ii]['index']=junction_one_blocks(blocks[ii])
	return(blocks)
##############################################################
### check and adjust monomer junctions
### one block
##############################################################
def junction_one_blocks(block):
		block['hits'][0]['junction']='first'
		for ii in range(1, block['nb_hits']):
			if block['hits'][ii]['envfromB'] == (block['hits'][ii-1]['envtoB'] + 1):
				block['hits'][ii]['junction'] = 'perfect_env-env'
			else:
				if block['hits'][ii]['envfromB'] == (block['hits'][ii-1]['alitoB'] + 1):
					block['hits'][ii]['junction'] = 'perfect_ali-env'
				else:
					if block['hits'][ii]['alifromB'] == (block['hits'][ii-1]['envtoB'] + 1):
						block['hits'][ii]['junction'] = 'perfect_env-ali'
					else:
						if block['hits'][ii]['alifromB'] == (block['hits'][ii-1]['alitoB'] + 1):
							block['hits'][ii]['junction'] = 'perfect_ali-ali'
						else:
							if block['hits'][ii]['envfromB'] > (block['hits'][ii-1]['envtoB'] + 1):
								block['hits'][ii]['junction'] = 'gap'
							else:
								block['hits'][ii]['junction'] = 'overlap'
		return(block)
##############################################################
### reverse all minus blocks
##############################################################
def reverse_all_minus_blocks(blocks):
	for ii in range(0, len(blocks)):
			blocks[ii]['index']=reverse_one_minus_block(blocks[ii])
	return(blocks)
	
##############################################################
### reverse one minus block
##############################################################
def reverse_one_minus_block(block):
	if block['hits'][0]['strand'] == 'plus':
		return(block)
	for ii in range(0, block['nb_hits']):
		#print(block['length'])
		tmp                           = block['length'] - block['hits'][ii]['alitoB'] + 1
		block['hits'][ii]['alitoB']   = block['length'] - block['hits'][ii]['alifromB'] + 1
		block['hits'][ii]['alifromB'] = tmp
		
		tmp                           = block['length'] - block['hits'][ii]['envtoB'] + 1
		block['hits'][ii]['envtoB']   = block['length'] - block['hits'][ii]['envfromB'] + 1
		block['hits'][ii]['envfromB'] = tmp
		
		#tmp = block['hits'][ii]['alifrom']
		#block['hits'][ii]['alifrom'] = block['hits'][ii]['alito']
		#block['hits'][ii]['alito'] =  tmp
		
		#tmp = block['hits'][ii]['envfrom']
		#block['hits'][ii]['envfrom'] = block['hits'][ii]['envto']
		#block['hits'][ii]['envto'] = tmp
		
		#print ("#####")
		#print(block['end'])
		#print(block['begin'])
	
		#print(block['hits'][ii]['envfrom'])
		#print(block['hits'][ii]['envto'])
		#print ("#####")
		#quit()
		
		
		
		
		#block['hits'][ii]['envfrom'] = block['end'] - block['hits'][ii]['envfrom'] + 1
		#block['hits'][ii]['envto']   = block['end'] - block['hits'][ii]['envto'] + 1
		#block['hits'][ii]['alifrom'] = block['end'] - block['hits'][ii]['alifrom'] + 1
		#block['hits'][ii]['alito']   = block['end'] - block['hits'][ii]['alito'] + 1
	block['hits']=list(reversed(block['hits']))
	return(block)
	
##############################################################
### index all blocks
##############################################################
def index_all_blocks(blocks):
	for ii in range(0, len(blocks)):
			# print ("##################")
			# print(ii)
			# print(blocks[ii])
			blocks[ii]['index']=index_one_block(blocks[ii])
	return(blocks)
	
##############################################################
### index one block
##############################################################
def index_one_block(block):
	#print(block['hits'][0]['alifromB'])
	
	output = [0] * (block['length']+1)
	for ii in range(1, block['hits'][0]['alifromB']):
		output[ii] = -1
	for kk in range(0, block['nb_hits']):
		ii=block['hits'][kk]['alifromB']
		#print("zobby: %d" %(block['hits'][kk]['alifromB']))
		index = block['hits'][kk]['hmmfrom']

#		for jj in range (1, len(block['hits'][kk]['monomer_seq'])):
		for jj in range (0, len(block['hits'][kk]['monomer_seq'])):
			if block['hits'][kk]['monomer_seq'][jj] == '.':
				output[ii] = -2
				ii = ii + 1
			else:
				#print("%d" % (jj))
				if block['hits'][kk]['target_seq'][jj] != '-':
					# print("--------")
					# print(output)
					# print(ii)
					# print(len(output))
					# print(block)
					# print(block['length'])
					output[ii] = index	
					# print ("xx")
					#print("un %d  %c" %(index, block['hits'][kk]['target_seq'][jj]))
					ii = ii + 1
				index=index+1
				
		if kk < (block['nb_hits']-1):
			if (block['hits'][kk+1]['envfromB'] - block['hits'][kk]['envtoB']) == 1:
				# envfromB - envtoB link
				while ii <= block['hits'][kk]['envtoB']:
					output[ii] = index
					#print("deux %d %c" %(index, block['hits'][kk]['target_seq'][jj]))
					index = index + 1
					ii = ii + 1
				index = 1
				while ii < block['hits'][kk+1]['alifromB']:
					output[ii] = index
					#print("trois %d %c" %(index, block['hits'][kk]['target_seq'][jj]))
					index = index + 1
					ii = ii + 1
	for ii in range(block['hits'][block['nb_hits']-1]['alitoB']+1, block['length']+1):
		output[ii] = -1
	return(output)
##############################################################
### print all blocks
##############################################################
def print_all_blocks(blocks):
	for ii in range(0, len(blocks)):
		print(blocks[ii]['nb_hits'])
		for jj in 	range(0, blocks[ii]['nb_hits']):
			print ("block: %d / %d " % (ii+1,  len(blocks)))
			print ("block: %d - %d (%d bp)  " % (blocks[ii]['begin'],  blocks[ii]['end'], blocks[ii]['length']))
			print ("hit: %d / %d"   % (jj+1, blocks[ii]['nb_hits']))
			print ("ini ->   coord: %d - %d - %s" % (blocks[ii]['hits'][jj]['envfrom'],  blocks[ii]['hits'][jj]['envto'],  blocks[ii]['hits'][jj]["strand"]))
			print ("env ->   coord: %d - %d     " % (blocks[ii]['hits'][jj]['envfromB'], blocks[ii]['hits'][jj]['envtoB']))
			print ("ali ->   coord: %d - %d     " % (blocks[ii]['hits'][jj]['alifromB'], blocks[ii]['hits'][jj]['alitoB']))
			print ("%s" % (blocks[ii]['hits'][jj]['monomer_seq']))
			print ("%s" % (blocks[ii]['hits'][jj]['target_seq']))
		print(blocs[ii]['index'])
		#~ alito=hit['alito']print_hit(blocks|ii][jj]
##############################################################
### print all blocks
##############################################################
def print_all_blocks3(blocks):
	for ii in range(0, len(blocks)):
		print(blocks[ii]['nb_hits'])
		for jj in 	range(0, blocks[ii]['nb_hits']):
			print ("%d / %d      %d - %d (%d bp)      %d / %d     %d - %d . %s %d - %d %d - %d  %s  %s" % (ii+1,  len(blocks), blocks[ii]['begin'],  blocks[ii]['end'], blocks[ii]['length'], jj+1, blocks[ii]['nb_hits'], blocks[ii]['hits'][jj]['envfrom'],  blocks[ii]['hits'][jj]['envto'],  blocks[ii]['hits'][jj]["strand"], blocks[ii]['hits'][jj]['envfromB'], blocks[ii]['hits'][jj]['envtoB'], blocks[ii]['hits'][jj]['alifromB'], blocks[ii]['hits'][jj]['alitoB'], blocks[ii]['hits'][jj]['monomer_seq'], blocks[ii]['hits'][jj]['target_seq']))
		#~ alito=hit['alito']print_hit(blocks|ii][jj]

##############################################################
### write all blocks into a file
##############################################################
def write_all_blocks(blocks, outfile, seqname=None):
	exist=False
	if  os.path.isfile(outfile): 
		exist=True
	handle=open(outfile,'a')
	if not exist:
		if seqname != None:
			handle.write("sequence ")
		handle.write("block block_nb block_begin block_end block_length")
		handle.write(" hit hit_nb envfrom envto strand envfromB envtoB alifromB alitoB junction")
		handle.write(" reference_seq monomer_seq\n")
	
	for ii in range(0, len(blocks)):
		for jj in 	range(0, blocks[ii]['nb_hits']):
			if seqname != None:
				handle.write(seqname+" ")			
			handle.write("%d %d %d %d %d %d %d %d %d %s %d %d %d %d %s %s %s\n" % (ii+1,  len(blocks), blocks[ii]['begin'],  blocks[ii]['end'], blocks[ii]['length'], jj+1, blocks[ii]['nb_hits'], blocks[ii]['hits'][jj]['envfrom'],  blocks[ii]['hits'][jj]['envto'],  blocks[ii]['hits'][jj]["strand"], blocks[ii]['hits'][jj]['envfromB'], blocks[ii]['hits'][jj]['envtoB'], blocks[ii]['hits'][jj]['alifromB'], blocks[ii]['hits'][jj]['alitoB'],  blocks[ii]['hits'][jj]['junction'], blocks[ii]['hits'][jj]['monomer_seq'], blocks[ii]['hits'][jj]['target_seq']))
	handle.close()
	return(0)	

##############################################################
### write all blocks for all sequences into a file
##############################################################
def write_all_blocks_all_seq(output, outfile):
	for target_name in 	output.keys():
		write_all_blocks(output[target_name], outfile, target_name)
	return(0)
	
##############################################################
### print all blocks
##############################################################
def print_all_blocks2(blocks):
	for ii in range(0, len(blocks)):
		print_one_block(blocks[ii])
		
##############################################################
### print one block
##############################################################
def print_one_block(block):
		print("########################")
		print("Nb of hits: %d" % block['nb_hits'])
		print("From-To: %d-%d" % (block['begin'],block['end']))
		print("Length: %d"     % block['length'])
		for jj in 	range(0, block['nb_hits']):
			print("#####")
			print ("hit: %d / %d      coord: %d - %d - %s" % (jj+1, block['nb_hits'], block['hits'][jj]['envfrom'], block['hits'][jj]['envto'], block['hits'][jj]["strand"]))
			print ("                  coord: %d - %d"      % (block['hits'][jj]['envfromB'], block['hits'][jj]['envtoB']))
			print ("                  coord: %d - %d"      % (block['hits'][jj]['alifromB'], block['hits'][jj]['alitoB']))
			print ("       %s" % (block['hits'][jj]['monomer_seq']))
			print ("       %s" % (block['hits'][jj]['target_seq']))
		#~ alito=hit['alito']print_hit(blocks|ii][jj]

##############################################################
### read hit list from tabular file
##############################################################
def read_hit_list(infile, sep=" "):
	out=[]
	with open(infile, "r") as f: 
		nb=0
		for line in f.readlines():
			li = line.lstrip()
			if not li.startswith("#"):
				nb=nb+1
				word = line.split(sep)
				if nb == 1:
					nbc=0
					name={}
					while nbc < len(word) and word[nbc] != '1':
						name[nbc]=word[nbc]
						nbc=nbc+1
				if nb > 1:
					tmp={}
					for ii in range(nbc):
						tmp[name[ii]]=word[ii]
					out.append(tmp)
	#print(out)
	return(out)


##############################################################
### get hit from hmmer output
##############################################################
def get_hit(line):				
	parsed_line = re.split("[ \n]+", line)
	nb=len(parsed_line)
	hit = {}
	hit["score"] = float(parsed_line[2])
	hit["evalue"] = float(parsed_line[4])
	hit["strand"] = "plus"
	hit["hmmfrom"] = int(parsed_line[5])
	hit["hmmto"] = int(parsed_line[6])
	hit["strand"] = "plus"
	hit["envfrom"] = int(parsed_line[11])
	hit["envto"] = int(parsed_line[12])
	hit["alifrom"] = int(parsed_line[8])
	hit["alito"] = int(parsed_line[9])
	if hit["alifrom"] > hit["alito"]:
		hit["strand"] = "minus"
		hit["alifrom"] = int(parsed_line[9])
		hit["alito"] = int(parsed_line[8])
		hit["envfrom"] = int(parsed_line[12])
		hit["envto"] = int(parsed_line[11])
	return(hit)		
	
				
##############################################################
### get hit from hmmer output
##############################################################
def get_hit2(file_iterator):
	file_iterator.newline()
	file_iterator.newline()
	line = file_iterator.newline()
	parsed_line = re.split("[ \n]+", line)
	nb=len(parsed_line)
	hit = {}
	hit["score"] = float(parsed_line[2])
	hit["evalue"] = float(parsed_line[4])
	hit["strand"] = "plus"
	hit["hmmfrom"] = int(parsed_line[5])
	hit["hmmto"] = int(parsed_line[6])
	hit["strand"] = "plus"
	hit["envfrom"] = int(parsed_line[11])
	hit["envto"] = int(parsed_line[12])
	hit["alifrom"] = int(parsed_line[8])
	hit["alito"] = int(parsed_line[9])
	if hit["alifrom"] > hit["alito"]:
		hit["strand"] = "minus"
		hit["alifrom"] = int(parsed_line[9])
		hit["alito"] = int(parsed_line[8])
		hit["envfrom"] = int(parsed_line[12])
		hit["envto"] = int(parsed_line[11])
	return(hit)					
	
		
##############################################################
### remove overlaping hits on different strands
##############################################################
def remove_overlaping_hits_plusminus(hits):
	nb_hits = len(hits)
	as_to_be_conserved = [True] * nb_hits
	for ii in range(0, nb_hits): 
		#as_to_be_conserved[ii]= is_best_overlapping_hit(hits, ii)
		#as_to_be_conserved[ii]= is_best_overlapping_hit_nostrandtest(hits, ii)
		as_to_be_conserved[ii]= is_best_overlapping_hit_new(hits, ii)
	#print(as_to_be_conserved)
	output = list(compress(hits, as_to_be_conserved))
	return(output)
##############################################################
### compare and return the best hit
##############################################################
def is_best_overlapping_hit(hits, hit_index, overlap=20):
	nb_hits = len(hits)
	ii = hit_index+1
	output = True
	while output == True and ii < nb_hits  and (hits[ii]['envfrom']-overlap)  <= hits[hit_index]['envto']:
		if hits[ii]['strand'] != hits[hit_index]['strand'] and hits[ii]['score'] > hits[hit_index]['score']:
			output = False
		ii = ii+1
	return(output)
##############################################################
### compare and return the best hit
##############################################################
def is_best_overlapping_hit_new(hits, hit_index, overlap=20):
	nb_hits = len(hits)
	ii = hit_index+1
	output = True
	#print("debug: %d %d %d" %(ii, hit_index, len(hits)))
	while output == True and ii < nb_hits  and (hits[ii]['envfrom']+overlap)  <= hits[hit_index]['envto']:
		if hits[ii]['score'] > hits[hit_index]['score']:
			output = False
		#print("xxxxx: %d %d %d" %(ii, hit_index, len(hits)))
		ii = ii+1
	ii = hit_index-1
	#print("debug: %d %d %d" %(ii, hit_index, len(hits)))
	while output == True and ii >= 0  and (hits[ii]['envto']-overlap)  >= hits[hit_index]['envfrom']:
		if hits[ii]['score'] > hits[hit_index]['score']:
			output = False
		#print("yyyyy: %d %d %d" %(ii, hit_index, len(hits)))
		ii = ii-1
	return(output)
##############################################################
### compare and return the best hit
##############################################################
def is_best_overlapping_hit_nostrandtest(hits, hit_index, overlap=20, verbose=0):
	#print("in")
	nb_hits = len(hits)
	#ii = hit_index+1
	ii = 0
	output = True
	while output == True and ii < nb_hits: 
#		if (hits[ii]['envfrom']-overlap)  <= hits[hit_index]['envto']:
		if ii != hit_index and ((hits[ii]['envfrom']+overlap)  <= hits[hit_index]['envto'] and (hits[ii]['envto']-overlap)  >= hits[hit_index]['envfrom']) :
			if verbose > 8:
				print("################")
				print(" removing one of the overlapping hits")
				print(hits[ii])
				print(hits[hit_index])
				print(ii)
				print(hit_index)
				print("################")				
			if hits[ii]['score'] > hits[hit_index]['score']:
				output = False
		ii = ii+1
	#print("out")
	return(output)
	
##############################################################
### split hmmer output into blocks
##############################################################
def split_blocks(hits, seqname,distance=50):
	nb_hits = len(hits)
	output = []

	block = {}
	block['sequence' ]= seqname
	block['hits']     = []
	block['hits'].append(hits[0])
	block['begin']   = hits[0]['envfrom']
	block['end']     = hits[0]['envto']
	block['nb_hits'] = 1
#	print(block)
	for ii in range(1, nb_hits): 
		# not on the same strand
		# or
		# separated by more than distance
		if hits[ii-1]['strand'] != hits[ii]['strand']  or  hits[ii-1]['envto'] < (hits[ii]['envfrom'] - distance):
			block['length'] = block['end'] - block['begin'] + 1
			block = get_hit_positions_along_block(block)
			output.append(block)
			block = {}
			block['hits'] = []
			block['hits'].append(hits[ii])
			block['begin'] = hits[ii]['envfrom']
			block['end'] = hits[ii]['envto']
			block['nb_hits'] = 1
		else:
			block['hits'].append(hits[ii])
			block['end'] = hits[ii]['envto']
			block['nb_hits'] = block['nb_hits'] + 1
	block['length'] = block['end'] - block['begin'] + 1
	block = get_hit_positions_along_block(block)
	output.append(block)
	return(output)
	
##############################################################
### calculate hit positions along the blocks instead of 
###     positions along the sequence
##############################################################
def get_hit_positions_along_block(block):
	output=block
	for ii in range(0, block['nb_hits']):
		block['hits'][ii]['alifromB'] = block['hits'][ii]['alifrom'] - block['begin'] + 1
		block['hits'][ii]['alitoB']   = block['hits'][ii]['alito']   - block['begin'] + 1
		block['hits'][ii]['envfromB'] = block['hits'][ii]['envfrom'] - block['begin'] + 1
		block['hits'][ii]['envtoB']   = block['hits'][ii]['envto']   - block['begin'] + 1	
	return(block)

##############################################################
### export junction sequences
##############################################################
def export_junctions_put_all_seq_into_memory(hmm_hits, seq_file):
	junction_file = tempfile.NamedTemporaryFile().name
	record_dict = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))
	all_junctions =  []
	nb_junctions = 0
	for target_name in hmm_hits.keys():
		for hits in hmm_hits[target_name]:
			#hits['hits'] = sorted(hits['hits'], key=itemgetter('alito')) 
			hits['hits'] = sorted(hits['hits'], key=lambda k: k['alito']) 
			minus_previous_is_gap_or_overlap = False
			for hit in hits['hits']:
				if hits['hits'][0]['strand'] == 'plus' and (hit['junction'] == 'gap' or hit['junction'] == 'overlap'):
				#if hits['hits'][0]['strand'] == 'plus' and (hit['junction'] == 'gap' or hit['junction'] == 'overlap' or hit['junction'] == 'perfect_ali-env'):
					junction_frag = record_dict[target_name].seq[(previous_alifrom - 1):(hit['envto'] - 1 + 1)]		
					record = SeqRecord(junction_frag, '%s_%d_%d_%d_%s' % (target_name, hits['begin'], hits['end'], previous_alifrom, hits['hits'][0]['strand']))
					nb_junctions = nb_junctions +1
					all_junctions.append(record)
					
				if hits['hits'][0]['strand'] == 'minus' and minus_previous_is_gap_or_overlap == True:
					junction_frag = record_dict[target_name].seq[(previous_alifrom - 1):(hit['envto'] - 1 + 1)]		
					record = SeqRecord(junction_frag, '%s_%d_%d_%d_%s' % (target_name, hits['begin'], hits['end'], previous_alifrom, hits['hits'][0]['strand']))
					nb_junctions = nb_junctions +1
					all_junctions.append(record)
				minus_previous_is_gap_or_overlap = False
					
				if hits['hits'][0]['strand'] == 'minus' and (hit['junction'] == 'gap' or hit['junction'] == 'overlap'):
				#if hits['hits'][0]['strand'] == 'minus' and (hit['junction'] == 'gap' or hit['junction'] == 'overlap' or hit['junction'] == 'perfect_ali-env'):
					minus_previous_is_gap_or_overlap = True
				
				previous_alifrom = hit['envfrom']
				
	SeqIO.write(all_junctions, junction_file, "fasta") 
	return(junction_file, nb_junctions)
	
	
def export_junctions(hmm_hits, seq_file):
	junction_file = tempfile.NamedTemporaryFile().name
	all_junctions =  []
	nb_junctions = 0
	with open(seq_file, "rU") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			target_name = record.id
			if target_name in hmm_hits.keys():
				for hits in hmm_hits[target_name]:
					#hits['hits'] = sorted(hits['hits'], key=itemgetter('alito')) 
					hits['hits'] = sorted(hits['hits'], key=lambda k: k['alito']) 
					minus_previous_is_gap_or_overlap = False
					for hit in hits['hits']:
						if hits['hits'][0]['strand'] == 'plus' and (hit['junction'] == 'gap' or hit['junction'] == 'overlap'):
						#if hits['hits'][0]['strand'] == 'plus' and (hit['junction'] == 'gap' or hit['junction'] == 'overlap' or hit['junction'] == 'perfect_ali-env'):
# envto    replaced by alito							
							junction_frag = record.seq[(previous_alifrom - 1):(hit['alito'] - 1 + 1)]		
							data = SeqRecord(junction_frag, '%s_%d_%d_%d_%d_%s' % (target_name,  hits['begin'], hits['end'], previous_alifrom, hit['alito'] - 1 + 1, hits['hits'][0]['strand']))
							nb_junctions = nb_junctions +1
							all_junctions.append(data)
 # envto   replaced by alito							

						if hits['hits'][0]['strand'] == 'minus' and minus_previous_is_gap_or_overlap == True:
							junction_frag = record.seq[(previous_alifrom - 1):(hit['alito'] - 1 + 1)].reverse_complement()	
							#print("zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz")
							data = SeqRecord(junction_frag, '%s_%d_%d_%d_%d_%s' % (target_name,  hits['begin'], hits['end'], previous_alifrom, hit['alito'] - 1 + 1, hits['hits'][0]['strand']))
							nb_junctions = nb_junctions +1
							all_junctions.append(data)
						minus_previous_is_gap_or_overlap = False
							
						if hits['hits'][0]['strand'] == 'minus' and (hit['junction'] == 'gap' or hit['junction'] == 'overlap'):
						#if hits['hits'][0]['strand'] == 'minus' and (hit['junction'] == 'gap' or hit['junction'] == 'overlap' or hit['junction'] == 'perfect_ali-env'):
							minus_previous_is_gap_or_overlap = True
						
 # envfrom replaced by alito							
						previous_alifrom = hit['alifrom']
				
	SeqIO.write(all_junctions, junction_file, "fasta") 
	return(junction_file, nb_junctions)	
#############################################################
### export junction sequences
##############################################################
def export_junctions_before_after_put_all_seq_into_memory(hmm_hits_before, hmm_hits_after, seq_file):
	junction_file = tempfile.NamedTemporaryFile().name
	f = open(junction_file, "w")
	record_dict = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))
	all_junctions =  []
	nb_junctions = 0
	for target_name in hmm_hits_before.keys():
		for hits in hmm_hits_before[target_name]:
			#hits['hits'] = sorted(hits['hits'], key=itemgetter('alito')) 
			hits['hits'] = sorted(hits['hits'], key=lambda k: k['alito']) 
			minus_previous_is_gap_or_overlap = False
			index_after = get_index_after(hmm_hits_after, target_name, hits['begin'])
			for hit in hits['hits']:
				if hits['hits'][0]['strand'] == 'plus' and (hit['junction'] == 'gap' or hit['junction'] == 'overlap'):
				#if hits['hits'][0]['strand'] == 'plus' and (hit['junction'] == 'gap' or hit['junction'] == 'overlap' or hit['junction'] == 'perfect_ali-env'):
					#junction_frag = record_dict[target_name].seq[(previous_alifrom - 1):(hit['envto'] - 1 + 1)]		
					nb_junctions = nb_junctions +1
					before = ','.join(map(str,hits['index'][(previous_alifrom-hits['begin']):(hit['envto']-hits['begin'])]))
					after = ','.join(map(str,  index_after[(previous_alifrom-hits['begin']):(hit['envto']-hits['begin'])]))
					f.write('->'+before+'\n')
					f.write('=>'+after+'\n')
					f.write('++> '+before+' '+after+'\n')
					
				if hits['hits'][0]['strand'] == 'minus' and minus_previous_is_gap_or_overlap == True:
					junction_frag = record_dict[target_name].seq[(previous_alifrom - 1):(hit['envto'] - 1 + 1)]		
					#record = SeqRecord(junction_frag, '%s_%d_%s' % (target_name, previous_alifrom, hits['hits'][0]['strand']))
					nb_junctions = nb_junctions +1
					f.write('%s_%d_%s\n' % (target_name, previous_alifrom, hits['hits'][0]['strand']))
					
					before = ','.join(map(str,hits['index'][(previous_alifrom-hits['begin']):(hit['envto']-hits['begin'])]))
					after = ','.join(map(str,  index_after[(previous_alifrom-hits['begin']):(hit['envto']-hits['begin'])]))
					f.write('->'+before+'\n')
					f.write('=>'+after+'\n')
					f.write('++> '+before+' '+after+'\n')
					
				minus_previous_is_gap_or_overlap = False
					
				if hits['hits'][0]['strand'] == 'minus' and (hit['junction'] == 'gap' or hit['junction'] == 'overlap'):
				#if hits['hits'][0]['strand'] == 'minus' and (hit['junction'] == 'gap' or hit['junction'] == 'overlap' or hit['junction'] == 'perfect_ali-env'):
					minus_previous_is_gap_or_overlap = True
				
				previous_alifrom = hit['envfrom']
	f.close()			
	SeqIO.write(all_junctions, junction_file, "fasta") 
	return(junction_file, nb_junctions)
	
def export_junctions_before_after(hmm_hits_before, hmm_hits_after, seq_file):
	junction_file = tempfile.NamedTemporaryFile().name
	f = open(junction_file, "w")
	all_junctions =  []
	nb_junctions = 0
	with open(seq_file, "rU") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			target_name = record.id
			if target_name in hmm_hits_before.keys():
				for hits in hmm_hits_before[target_name]:
					#hits['hits'] = sorted(hits['hits'], key=itemgetter('alito')) 
					hits['hits'] = sorted(hits['hits'], key=lambda k: k['alito']) 
					minus_previous_is_gap_or_overlap = False
					index_after = get_index_after(hmm_hits_after, target_name, hits['begin'])
					for hit in hits['hits']:
						if hits['hits'][0]['strand'] == 'plus' and (hit['junction'] == 'gap' or hit['junction'] == 'overlap'):
						#if hits['hits'][0]['strand'] == 'plus' and (hit['junction'] == 'gap' or hit['junction'] == 'overlap' or hit['junction'] == 'perfect_ali-env'):
							#junction_frag = record_dict[target_name].seq[(previous_alifrom - 1):(hit['envto'] - 1 + 1)]		
							nb_junctions = nb_junctions +1
							before = ','.join(map(str,hits['index'][(previous_alifrom-hits['begin']):(hit['envto']-hits['begin'])]))
							after = ','.join(map(str,  index_after[(previous_alifrom-hits['begin']):(hit['envto']-hits['begin'])]))
							f.write('->'+before+'\n')
							f.write('=>'+after+'\n')
							f.write('++> '+before+' '+after+'\n')
							
						if hits['hits'][0]['strand'] == 'minus' and minus_previous_is_gap_or_overlap == True:
							junction_frag = record.seq[(previous_alifrom - 1):(hit['envto'] - 1 + 1)]		
							#record = SeqRecord(junction_frag, '%s_%d_%s' % (target_name, previous_alifrom, hits['hits'][0]['strand']))
							nb_junctions = nb_junctions +1
							f.write('%s_%d_%s\n' % (target_name, previous_alifrom, hits['hits'][0]['strand']))
							
							before = ','.join(map(str,hits['index'][(previous_alifrom-hits['begin']):(hit['envto']-hits['begin'])]))
							after = ','.join(map(str,  index_after[(previous_alifrom-hits['begin']):(hit['envto']-hits['begin'])]))
							f.write('->'+before+'\n')
							f.write('=>'+after+'\n')
							f.write('++> '+before+' '+after+'\n')
							
						minus_previous_is_gap_or_overlap = False
							
						if hits['hits'][0]['strand'] == 'minus' and (hit['junction'] == 'gap' or hit['junction'] == 'overlap'):
						#if hits['hits'][0]['strand'] == 'minus' and (hit['junction'] == 'gap' or hit['junction'] == 'overlap' or hit['junction'] == 'perfect_ali-env'):
							minus_previous_is_gap_or_overlap = True
						
						previous_alifrom = hit['envfrom']
	handle.close()
	f.close()			
	#SeqIO.write(all_junctions, junction_file, "fasta") 
	return(junction_file, nb_junctions)
#############################################################
### get index after ...
##############################################################

	
def get_index_after(hmm_hits_after, target_name, begin):
		for hitsss in hmm_hits_after[target_name]:
			if hitsss['begin'] == begin:
				return(hitsss['index'])

#############################################################
### parse options ...
##############################################################
def parse_options(args):
	args_error = 0
	##############################
	if args.action == 'search':
		if args.seq_file == None:
			args_error = 1
			print('error: the following argument is required if action is search: -s')
			quit(args_error)
		else:	
			if not os.path.exists(args.seq_file):
				args_error = 2
				print("error: the following sequence file doesn't exist: "+args.seq_file)
				quit(args_error)
		if args.ref_file == None:
			args_error = 3
			print('error: the following argument is required if action is search: -r')
			quit(args_error)
		else:	
			if not os.path.exists(args.ref_file):
				args_error = 4
				print("error: the following monomer file doesn't exist: "+args.ref_file)
				quit(args_error)
		if args.minimal_score < 0:
			args_error = 8
			print("the minimal score should be greater or equal to 0 : "+str(args.minimal_score))	
		if args.phase != None:
			print('warning: the following argument will be ignored if action is search: -p')
		args.out_tab_file = args.seq_file+'_'+args.suffix+'.tab.dat'
		args.out_block_seq_index_file = args.seq_file+'_'+args.suffix+'.block_seq_index.fst'
		args.out_block_seq_index_file_before_merge = args.seq_file+'_'+args.suffix+'.block_seq_index_before_merge.fst'
		args.out_block_seq_file = args.seq_file+'_'+args.suffix+'.block_seq.fst'
	
	##############################
	if args.action == 'search_and_extract':
		if args.seq_file == None:
			args_error = 1
			print('error: the following argument is required if action is search_and_extract: -s')
			quit(args_error)
		else:	
			if not os.path.exists(args.seq_file):
				args_error = 2
				print("error: the following sequence file doesn't exist: "+args.seq_file)
				quit(args_error)
		if args.ref_file == None:
			args_error = 3
			print('error: the following argument is required if action is search_and_extract: -r')
			quit(args_error)
		else:	
			if not os.path.exists(args.ref_file):
				args_error = 4
				print("error: the following monomer file doesn't exist: "+args.ref_file)
				quit(args_error)
		if args.phase == None:
			args_error = 7
			print('error: the following argument is required if action is search_and_extract: -p')
			quit(args_error)
		else:	
			if args.phase <= 0:
				args_error = 8
				print("error: the phase value should be greater than 0 (and lower than monomer length, not tested): "+str(args.phase))
				quit(args_error)
		if args.minimal_score < 0:
			args_error = 8
			print("error: the minimal score should be greater or equal to 0 : "+str(args.minimal_score))	
			quit(args_error)
		if args.index_file != None:
			print('warning: the following argument will be ignored if action is search_and_extract: -i')			
		args.out_tab_file = args.seq_file+'_'+args.suffix+'.tab.dat'
		args.out_block_seq_index_file = args.seq_file+'_'+args.suffix+'.block_seq_index.fst'
		args.out_block_seq_index_file_before_merge = args.seq_file+'_'+args.suffix+'.block_seq_index_before_merge.fst'
		args.out_block_seq_file = args.seq_file+'_'+args.suffix+'.block_seq.fst'
		args.out_monomer_seq_index_file = args.out_block_seq_index_file+'.monomer_seq_index_'+str(args.phase)+'.fst'
		args.out_monomer_seq_file = args.out_block_seq_index_file+'.monomer_seq_'+str(args.phase)+'.fst'
		args.out_monomer_gtf_file = args.out_block_seq_index_file+'.monomer_seq_'+str(args.phase)+'.gtf'
		args.out_block_seq_index_file_before_merge = args.seq_file+'_'+args.suffix+'.block_seq_index_before_merge.fst'
	
	##############################
	if args.action == 'extract':
		if args.seq_file == None:
			args_error = 1
			print('error: the following argument is required if action is extract: -s')
			quit(args_error)
		else:	
			args.out_block_seq_index_file = args.seq_file+'_'+args.suffix+'.block_seq_index.fst'
			if not os.path.exists(args.out_block_seq_index_file):
				args_error = 6
				print("error: the following index file doesn't exist: "+args.out_block_seq_index_file)		
				quit(args_error)
		if args.phase == None:
			args_error = 7
			print('error: the following argument is required if action is extract: -p')
			quit(args_error)
		else:	
			if args.phase <= 0:
				args_error = 8
				print("error: the phase value should be greater than 0 (and lower than monomer length, not tested): "+str(args.phase))
				quit(args_error)
		if args.ref_file != None:
			print('warning: the following argument will be ignored if action is extract: -r')	
		if args.minimal_score != 0:
			print('warning: the following argument will be ignored if action is search: -m')
		args.out_monomer_seq_index_file = args.out_block_seq_index_file+'.monomer_seq_index_'+str(args.phase)+'.fst'
		args.out_monomer_seq_file = args.out_block_seq_index_file+'.monomer_seq_'+str(args.phase)+'.fst'
		args.out_monomer_gtf_file = args.out_block_seq_index_file+'.monomer_seq_'+str(args.phase)+'.gtf'
	return(args)
############################################################
### extrapolation of index at 5' and 3' extremities ...
##############################################################
def extrapol_index(hmm_hits, monomer_length, extrapol_length=5):
		for target_name in hmm_hits.keys():
			for hits in hmm_hits[target_name]:
					ok = 0
					old_index =  hits['index'].copy()
					
					# 5' extremity
					ii = 1
					while ii < len(hits['index']) and hits['index'][ii] == -1:
						ii = ii + 1
					beg_ii = ii
					val =  hits['index'][ii]

					ii  = ii -1
					while ii > 0 and ii > (beg_ii - extrapol_length) and hits['index'][ii] == -1:
						ok = 1
						val = val - 1
						if val > 0:
							hits['index'][ii] = val 
						else: 
							val = monomer_length
							hits['index'][ii] = val 
						ii = ii - 1


					# 3' extremity
					ii = beg_ii
					while ii <  len(hits['index']) and hits['index'][ii] != -1:
						val = hits['index'][ii]
						ii = ii + 1
					beg_ii = ii
					
					while ii < len(hits['index']) and ii < (beg_ii + extrapol_length) and hits['index'][ii] == -1:
						ok = 1
						val = val + 1
						if val <= monomer_length:
							hits['index'][ii] = val 
						else:
							val = 1
							hits['index'][ii] = val 
						ii = ii + 1
						
					# if ok == 1:
						# print(target_name)
						# print (old_index)
						# print (hits['index'])
						# print ("##########################")
		return(hmm_hits)
