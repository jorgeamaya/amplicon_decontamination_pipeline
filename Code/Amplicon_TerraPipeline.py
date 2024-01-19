#!/usr/bin/env python

import os
import sys
import argparse
import json
import subprocess
#import multiprocessing

import amplicon_decontamination as ad
import asv_to_cigar as ac

def main():
	"""
	Implementation of the amplicon decontamination pipeline for use and 
	distribution in TERRA v=1.0

	Usage: python Code/Amplicon_TerraPipeline.py --config config.json ...
	
	Returns:
	None. However, see the functions documentation.
	"""

	### LOAD ARGUMENTS

	#Parse command line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('--config', help="Path to config.json file.", required =True)
	parser.add_argument('--overlap_reads', action="store_true", help="Specify whether the data contains paired reads that overlap on their targets. For example, MiSeq paired reads of 250bp for 500bp or shorter targets. The pipeline will perform the contamination detection steps with bbmerge or dada2, depending on the choice made be the user.")
	parser.add_argument('--mixed_reads', action="store_true", help="Specify whether the data contains a mix of reads that overlap and do not overlap on their targets. For example, iSeq paired reads of 150bp for targets that are shorter than 300bp long and targets that are longer than 300bp.")
	parser.add_argument('--terra', action="store_true", help="Specify whether the pipeline is being run in the Terra platform.")
	parser.add_argument('--meta', action="store_true", help="Specify if metadata file must be created. This flag runs the pipeline from the beginning.")
	parser.add_argument('--repo', action="store_true", help="Specify if the reports must be created.")
	parser.add_argument('--bbmerge', action="store_true", help="Specify if preprocess merge with bbmerge is performed.")
	parser.add_argument('--bbmerge_report', action="store_true", help="Aggregate bbmerge reports and generate a visual summary. This flag will only work if the --merge and --overlap_reads flags are active.")
	parser.add_argument('--adaptor_removal', action="store_true", help="Specify if adaptor removal needed.")
	parser.add_argument('--primer_removal', action="store_true", help="Specify if primer removal needed.")
	parser.add_argument('--dada2', action="store_true", help="Specifiy if standard preprocess merge with DADA2 is performed.")
	#parser.add_argument('--dada2_contamination', action="store_true", help="Specifiy if contamination detection is performed with DADA2. Unlike in the standard DADA2 preproprocess, dada2_contamination merges the reads after adaptor removal, not after primer removal. Primer removal will cut the barcodes, making it impossible to perform the matching procedure for contamination detection.")
	parser.add_argument('--postproc_dada2', action="store_true", help="Specifiy if postProcess of DADA2 results is perfomed.")
	parser.add_argument('--asv_to_cigar', action="store_true", help="Specifiy if the ASV to CIGAR transformation is perfomed.")

	args = parser.parse_args()

	#Check minimum arguments and contracdicting flags
	if not any([args.mixed_reads, args.overlap_reads]):
		sys.exit('Pre-process halted: User must declare the mixed_reads or overlap_reads flags.')

	if args.mixed_reads and args.bbmerge and args.bbmerge_report:
		sys.exit('Merging of non-overlapping reads with bbmerge unsupported. Use the dada2 and dada2_contamination flags instead. Read documentation for limitations of this method.')

	if args.terra:
		print("Pipeline is running in Terra. Adjusted paths will be used.")
	else:
		print("Pipeline not running in Terra. Default paths will be used.")
				
	#Configuration aguments will be parsed from config.json
	with open(args.config, 'r') as config_file:
		config_inputs = json.load(config_file)
		path_to_fq = config_inputs['path_to_fq']
		path_to_flist = config_inputs['path_to_flist']
		if 'pattern_fw' in config_inputs.keys(): pattern_fw = config_inputs['pattern_fw']
		if 'pattern_rv' in config_inputs.keys(): pattern_rv = config_inputs['pattern_rv']
		if 'read_maxlength' in config_inputs.keys(): read_maxlength = config_inputs['read_maxlength']
		if 'pairread_minlength' in config_inputs.keys(): pairread_minlength = config_inputs['pairread_minlength']
		if 'merge_minlength' in config_inputs.keys(): merge_minlength = config_inputs['merge_minlength']
		if 'pr1' in config_inputs.keys(): pr1 = config_inputs['pr1']
		if 'pr2' in config_inputs.keys(): pr2 = config_inputs['pr2']
		if 'Class' in config_inputs.keys(): Class = config_inputs['Class']
		if 'maxEE' in config_inputs.keys(): maxEE = config_inputs['maxEE']
		if 'trimRight' in config_inputs.keys(): trimRight = config_inputs['trimRight']
		if 'minLen' in config_inputs.keys(): minLen = config_inputs['minLen']
		if 'truncQ' in config_inputs.keys(): truncQ = config_inputs['truncQ']
		if 'matchIDs' in config_inputs.keys(): matchIDs = config_inputs['matchIDs']
		if 'max_consist' in config_inputs.keys(): max_consist = config_inputs['max_consist']
		if 'omegaA' in config_inputs.keys(): omegaA = config_inputs['omegaA']
		if 'saveRdata' in config_inputs.keys(): saveRdata = config_inputs['saveRdata']
		if 'justConcatenate' in config_inputs.keys(): justConcatenate = config_inputs['justConcatenate']
		if 'maxMismatch' in config_inputs.keys(): maxMismatch = config_inputs['maxMismatch']
		if 'overlap_pr1' in config_inputs.keys(): overlap_pr1 = config_inputs['overlap_pr1']
		if 'overlap_pr2' in config_inputs.keys(): overlap_pr2 = config_inputs['overlap_pr2']
		if 'reference' in config_inputs.keys(): reference = config_inputs['reference']
		if 'adjust_mode' in config_inputs.keys(): adjust_mode = config_inputs['adjust_mode']
		if 'path_to_snv' in config_inputs.keys(): path_to_snv = config_inputs['path_to_snv']
		if 'no_ref' in config_inputs.keys(): no_ref = config_inputs['no_ref']
		if 'reference2' in config_inputs.keys(): reference2 = config_inputs['reference2']
		if 'strain' in config_inputs.keys(): strain = config_inputs['strain']
		if 'strain2' in config_inputs.keys(): strain2 = config_inputs['strain2']
		if 'polyN' in config_inputs.keys(): polyN = int(config_inputs['polyN'])
		if 'min_reads' in config_inputs.keys(): min_reads = int(config_inputs['min_reads'])
		if 'min_samples' in config_inputs.keys(): min_samples = int(config_inputs['min_samples'])
		if 'max_snv_dist' in config_inputs.keys(): max_snv_dist = int(config_inputs['max_snv_dist'])
		if 'max_indel_dist' in config_inputs.keys(): max_indel_dist = int(config_inputs['max_indel_dist'])
		if 'include_failed' in config_inputs.keys(): include_failed = eval(config_inputs['include_failed'])
		if 'exclude_bimeras' in config_inputs.keys(): exclude_bimeras = eval(config_inputs['exclude_bimeras'])
		if 'amp_mask' in config_inputs.keys(): amp_mask = config_inputs['amp_mask']
		if 'verbose' in config_inputs.keys(): verbose = eval(config_inputs['verbose'])
		#barcodes_file = config_inputs['barcodes_file']

	### PREPARE OUTPUT DIRECTORIES

	#Generate the path to the Results directory.
	global res_dir
	global rep_dir
	res_dir = os.path.abspath(os.path.join("Results"))
	rep_dir = os.path.abspath(os.path.join("Report"))

	#Restart Results and Report directory if the workflow runs form the start
	if args.meta:
		os.system("rm -rf " + res_dir)
		os.mkdir(res_dir)

	if args.repo:
		os.system("rm -rf " + rep_dir)
		os.mkdir(rep_dir)
	
	#Create metadata files
	if args.meta:
		ad.flush_dir(res_dir, "Fq_metadata")
		ad.create_meta(path_to_fq, res_dir, "Fq_metadata", "rawfilelist.tsv", pattern_fw, pattern_rv)

		#List missing file from the complete list of files
		with open(path_to_flist, 'r') as input_file:
			next(input_file)
			samples = [line.split(',')[0] for line in input_file]

		with open('Results/Fq_metadata/rawfilelist.tsv', 'r') as raw_files:
			raw_file_samples = [line.split('\t')[0] for line in raw_files]

		missing_samples = [sample for sample in samples if sample not in raw_file_samples]

		with open('Results/missing_files.tsv', 'w') as output_file:
			output_file.write('\n'.join(missing_samples))

	### EXECUTE

	#Remove adaptors
	#Most sequences must be adaptor free; just in case, run this step to eliminate any lingering adaptors.
	if args.adaptor_removal:
		ad.flush_dir(res_dir, "AdaptorRem")
		meta = open(os.path.join(res_dir, "Fq_metadata", "rawfilelist.tsv"), 'r')
		samples = meta.readlines()
		
		for sample in samples:
			slist = sample.split()
			ad.adaptor_rem(slist[0], slist[1], slist[2], res_dir, "AdaptorRem")
	
		ad.create_meta(os.path.join(res_dir, "AdaptorRem"), res_dir, "AdaptorRem", "adaptorrem_meta.tsv",
			pattern_fw="*_val_1.fq.gz", pattern_rv="*_val_2.fq.gz")

	#Merge forward and reverse reads with bbmerge. Only for reads that overlap.
	if args.bbmerge and args.overlap_reads:
		ad.flush_dir(res_dir, "Merge")
		meta = open(os.path.join(res_dir, "AdaptorRem", "adaptorrem_meta.tsv") , 'r')
		samples = meta.readlines()
		for sample in samples:
			slist = sample.split()
			ad.mergereads(slist[0], slist[1], slist[2], res_dir, "Merge", read_maxlength, pairread_minlength, merge_minlength)

	#Process merge report and generate RStudio plots
	if args.bbmerge_report and args.overlap_reads:
		ad.flush_dir(rep_dir, "Merge")
		with open(os.path.join(rep_dir, "Merge", "bbmergefields.tsv"), 'a') as f:
			header = ['SampleID', 'Pairs', 'Joined', 'JoinedP', 'Ambiguous', 'AmbiguousP', 'No_Solution', 'No_SolutionP', 
				'Too_Short', 'Too_ShortP', 'Avg_Insert', 'Standard_Deviation', 'Mode', 'Insert_range_low', 
				'Insert_range_up', '90th_pc', '75th', '50th_pc', '25th_pc', '10th_pc']	
			f.write('\t'.join(header) + '\n')
		meta = open(os.path.join(res_dir, "Merge", "merge_meta.tsv"), 'r')
		samples = meta.readlines()
		
		for sample in samples:
			slist = sample.split()
			ad.extract_bbmergefields(slist[0], slist[1], slist[3], path_to_flist, res_dir, rep_dir, "Merge", args.terra)

	#Remove primers
	#For a set where all reads have overlap
	if args.overlap_reads and args.primer_removal:
		ad.flush_dir(res_dir, "PrimerRem")
		meta = open(os.path.join(res_dir, "AdaptorRem", "adaptorrem_meta.tsv"), 'r')
		samples = meta.readlines()
		
		for sample in samples:
			slist = sample.split()
			ad.trim_primer(slist[0], slist[1], slist[2], res_dir, "PrimerRem", pr1, pr2, "prim")

		ad.create_meta(os.path.join(res_dir,"PrimerRem"), res_dir, "PrimerRem", "primrem_meta.tsv",
			pattern_fw="*_prim_1.fq.gz", pattern_rv="*_prim_2.fq.gz")

	#For a set that mixes reads with and without overlap
	if args.mixed_reads and args.primer_removal:
		ad.flush_dir(res_dir, "PrimerRem")
		meta = open(os.path.join(res_dir, "AdaptorRem", "adaptorrem_meta.tsv"), 'r')

		#Trim primers off Overlapping short targets and demux them to a different file
		samples = meta.readlines()
		for sample in samples:
			slist = sample.split()
			ad.trim_primer(slist[0], slist[1], slist[2], res_dir, "PrimerRem", overlap_pr1, overlap_pr2, "mixed_op", True)

		#Metafile for trimmed overlapping target reads
		ad.create_meta(os.path.join(res_dir, "PrimerRem"), res_dir, "PrimerRem", "mixed_op_prim_meta.tsv",
			pattern_fw="*_mixed_op_1.fq.gz", pattern_rv="*_mixed_op_2.fq.gz")

		#Metafile for un-trimmed non-op target reads
		ad.create_meta(os.path.join(res_dir, "PrimerRem"), res_dir, "PrimerRem", "mixed_temp_meta.tsv",
			pattern_fw="*_temp_1.fq.gz", pattern_rv="*_temp_2.fq.gz")
		temp_meta = open(os.path.join(res_dir, "PrimerRem", "mixed_temp_meta.tsv"), 'r')

		#Trim primers off second subset of non-op long targets 
		samples = temp_meta.readlines()	
		for sample in samples:
			slist = sample.split()
			ad.trim_primer(slist[0], slist[1], slist[2], res_dir, "PrimerRem", pr1, pr2, "mixed_nop")

		#Metafile for trimmed non-op target reads
		ad.create_meta(os.path.join(res_dir, "PrimerRem"), res_dir, "PrimerRem", "mixed_nop_prim_meta.tsv", 
			pattern_fw="*_mixed_nop_1.fq.gz", pattern_rv="*_mixed_nop_2.fq.gz")
		temp_meta.close()

	#Perform denoising with DADA2
	#For a set where all reads have overlap
	if args.overlap_reads and args.dada2:
		ad.flush_dir(res_dir, "DADA2", "QProfile")
		path_to_meta = os.path.join(res_dir, "PrimerRem", "primrem_meta.tsv")
		ad.run_dada2(path_to_meta, path_to_fq, path_to_flist, Class, maxEE, trimRight, minLen, truncQ, matchIDs, max_consist, omegaA, justConcatenate, maxMismatch, saveRdata, res_dir, "DADA2", args.terra)
		cmd = ['cp', os.path.join(res_dir, 'DADA2', 'seqtab.tsv'), 
			os.path.join(res_dir, '.'), 
			'\\', 
			'cp', os.path.join(res_dir, 'DADA2', 'ASVBimeras.txt'),
                        os.path.join(res_dir, '.')
			]

		proccp = subprocess.Popen(cmd)
		proccp.wait()

	#For a set that mixes reads with and without overlap
	if args.mixed_reads and args.dada2:
		#Run DADA2 on op targets
		ad.flush_dir(res_dir, "DADA2_OP", "QProfile")
		path_to_meta = os.path.join(res_dir, "PrimerRem", "mixed_op_prim_meta.tsv")
		justConcatenate=0	
		ad.run_dada2(path_to_meta, path_to_fq, path_to_flist, Class, maxEE, trimRight, minLen, truncQ, matchIDs, max_consist, omegaA, justConcatenate, maxMismatch,saveRdata, res_dir, "DADA2_OP", args.terra)
		seqtab_op = os.path.join(res_dir, 'DADA2_OP', 'seqtab.tsv')
		bimera_op = os.path.join(res_dir, 'DADA2_OP', 'ASVBimeras.txt')

		#Run DADA2 on non-op targets
		ad.flush_dir(res_dir, "DADA2_NOP", "QProfile")
		path_to_meta = os.path.join(res_dir, "PrimerRem", "mixed_nop_prim_meta.tsv")
		justConcatenate=1	
		ad.run_dada2(path_to_meta, path_to_fq, path_to_flist, Class, maxEE, trimRight, minLen, truncQ, matchIDs, max_consist, omegaA, justConcatenate, maxMismatch,saveRdata, res_dir, "DADA2_NOP", args.terra)
		seqtab_nop = os.path.join(res_dir, 'DADA2_NOP', 'seqtab.tsv')
		bimera_nop = os.path.join(res_dir, 'DADA2_NOP', 'ASVBimeras.txt')

		#ASV modification block for non-op targets and merge two ASV tables
		if reference is not None:
			if args.terra:
				path_to_program = os.path.join("/", "Code/adjustASV.R")
			else:
				path_to_program = os.path.join("Code/adjustASV.R")
			adjASV = ['Rscript', path_to_program, '-s', seqtab_nop, '-ref', str(reference),
			'-dist', adjust_mode,
			'-o', os.path.join(res_dir, 'DADA2_NOP', 'correctedASV.txt')]
			print(adjASV)
			procASV = subprocess.Popen(adjASV)
			procASV.wait()
			seqtab_corrected = os.path.join(res_dir, 'DADA2_NOP', 'seqtab_corrected.tsv')
			seqtab = ad.merge_seqtab(seqtab_op, seqtab_corrected)
			bimera = ad.merge_bimeras(bimera_op, bimera_nop)
		else:
			print('--reference file not found. skipping ASV correction..')
			seqtab = ad.merge_seqtab(seqtab_op, seqtab_nop)
			bimera = ad.merge_bimeras(bimera_op, bimera_nop)

		seqtab.to_csv(os.path.join(res_dir, 'seqtab.tsv'), sep = "\t", index=False)
		bimera.to_csv(os.path.join(res_dir, 'ASVBimeras.txt'), sep = "\t", index=False)

	if args.postproc_dada2:		
		ad.flush_dir(res_dir, "PostProc_DADA2")
		
		path_to_seqtab = os.path.join(res_dir, 'seqtab.tsv')

		if args.terra:
			path_to_program = os.path.join("/", "Code/postProc_dada2.R")
		else:
			path_to_program = os.path.join("Code/postProc_dada2.R")

		postProc = ['Rscript', path_to_program, 
				'-s', path_to_seqtab, 
				'-b', os.path.join(res_dir, 'ASVBimeras.txt'),
				'-snv', os.path.join(path_to_snv),
				'--indel_filter', '0.895',
				'-o', os.path.join(res_dir, 'PostProc_DADA2', 'ASVTable.txt'),
				'--fasta']

		if no_ref == 'True':
			postProc.extend(['-no_ref'])
		else:
			postProc.extend(['--reference', str(reference), '--strain', strain])
			if reference2 != "":
				postProc.extend(['--reference2', str(reference2), '--strain2', strain2])

		print(postProc)
		procASV = subprocess.Popen(postProc)
		procASV.wait()

	#ASV to CIGAR
	#Convert ASVs from DADA2 pipeline to pseudo-CIGAR strings.
	if args.asv_to_cigar:		
		ad.flush_dir(res_dir, "ASV_to_CIGAR", "alingments")

		path_to_seqtab = os.path.join(res_dir, 'seqtab.tsv')
		path_to_fasta = os.path.join(res_dir, "PostProc_DADA2", "ASVSeqs.fasta") #Fasta file of ASV sequences from DADA2 pipeline"
		path_to_table = os.path.join(res_dir, "PostProc_DADA2", "ASVTable.txt") #ASV table from DADA2 pipeline
		path_to_out = os.path.join(res_dir, "CIGARVariants_Bfilter.out.tsv") #Output seqtab tsv file with amplicon/variant counts
		path_asv_to_cigar = os.path.join(res_dir, "ASV_to_CIGAR", "ASV_to_CIGAR.out.txt") #Output file for ASV -> CIGAR string table 
		path_to_amp_db = reference #Amplicon sequence fasta file
		path_to_alignments = os.path.join(res_dir, "ASV_to_CIGAR", "alingments") #Directory to store ASV alignment files

		print(f"INFO: Loading {path_to_amp_db}")
		amplicons = ac.parse_amp_db(path_to_amp_db)
		if not amplicons:
			print(f"ERROR: No amplicons in {path_to_amp_db}")
			sys.exit(1)

		#Disabled. Possibly deprecated
		if os.path.exists(amp_mask):
			print(f"INFO: Loading {amp_mask}")
			mask = ac.parse_dustmasker(amp_mask)
		else:
			print(f"INFO: No mask data specified.")
			mask = {}

		print(f"INFO: Loading {path_to_fasta}")
		asvs = ac.get_asv_seqs(path_to_fasta)
		if not asvs:
			print(f"ERROR: No ASV sequences in {path_to_fasta}")
			sys.exit(1)

		print(f"INFO: Parsing {path_to_table} with total reads >= {min_reads}, samples >= {min_samples}, snv_dist <= {max_snv_dist}, indel_dist <= {max_indel_dist}")

		if include_failed:
			print("WARNING: Including ASVs that failed post-DADA2 filters! This is not recommended.")
		else:
			print("INFO: Excluding ASVs that failed post-DADA2 filters.")

		if exclude_bimeras:
			print("INFO: Excluding ASVs that DADA2 marked as bimeras.")

		bins = ac.parse_asv_table(path_to_table, min_reads=min_reads, min_samples=min_samples, max_snv_dist=max_snv_dist, max_indel_dist=max_indel_dist, include_failed=include_failed, exclude_bimeras=exclude_bimeras) #This function only matches to the first strain.
		if not bins:
			print(f"ERROR: No useable data in {path_to_table}")
			sys.exit(1)

		print(f"INFO: Writing amplicon fasta files to {path_to_alignments}")
		ac.write_amplicon_fastas(asvs, bins, amplicons, outdir=path_to_alignments)

		print("INFO: Running MUSCLE aligner on amplicon fasta files. Please wait...")
		ac.run_muscle(bins, outdir=path_to_alignments)

		print("INFO: Parsing alignments to CIGAR strings")
		cigars = ac.parse_alignments(bins, mask=mask, min_homopolymer_length=polyN, outdir=path_to_alignments, verbose=False)
		if not cigars:
			print("ERROR: could not determine CIGAR strings")
			sys.exit(1)

		if path_asv_to_cigar:
			ac.write_cigar_strings(cigars, path_asv_to_cigar)
			print(f"INFO: Wrote ASV->CIGAR table to {path_asv_to_cigar}")

		print(f"INFO: Converting DADA2 seqtab file {path_to_seqtab} to {path_to_out}")
		if ac.convert_seqtab(path_to_seqtab, cigars, path_to_out):
			print("INFO: Completed successfully!")

if __name__ == "__main__":
	main()
