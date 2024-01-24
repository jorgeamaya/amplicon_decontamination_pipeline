version 1.0

workflow mixed_reads_ampseq {
	input {	
		Boolean CI
		Boolean AMPSEQ

		#General commands
		String type_of_reads
		String path_to_fq 
		File path_to_flist
		String pattern_fw = "*_L001_R1_001.fastq.gz"
		String pattern_rv = "*_L001_R2_001.fastq.gz"

		#Commands for AmpSeq
		File pr1
		File pr2
		String Class = "parasite"
		String maxEE = "5,5"
		String trimRight = "0,0"
		Int minLen = 30
		String truncQ = "5,5"
		String matchIDs = "0"
		Int max_consist = 10
		Float omegaA = 0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
		String saveRdata = ""
		Int justConcatenate = 0
		Int maxMismatch = 0
		File overlap_pr1
		File overlap_pr2
		File path_to_snv
		String no_ref = 'False'
		File reference
		String adjust_mode = "absolute"
		File reference2
		String strain = "3D7"
		String strain2 = "DD2"
		String polyN = "5"
		String min_reads = "0"
		String min_samples = "0"
		String max_snv_dist = "-1"
		String max_indel_dist = "-1"
		String include_failed = "False"
		String exclude_bimeras = "False"
		String amp_mask = "None"

		#Command for the decontamination pipeline
		Int read_maxlength = 200
		Int pairread_minlength = 100
		Int merge_minlength = 100
		Int joined_threshold = 1000
		Float contamination_threshold = 0.5
		String verbose = "False"		
	}

	if (CI) {
		call inline_barcodes_process {
			input: 	
				type_of_reads = type_of_reads,
				path_to_fq = path_to_fq,
				path_to_flist = path_to_flist,
				pattern_fw = pattern_fw,
				pattern_rv = pattern_rv,
				pr1 = pr1,
				pr2 = pr2,
				Class = Class,
				maxEE = maxEE,
				trimRight = trimRight,
				minLen = minLen,
				truncQ = truncQ,
				matchIDs = matchIDs,
				max_consist = max_consist,
				omegaA = omegaA,
				saveRdata = saveRdata,
				justConcatenate = justConcatenate,
				maxMismatch = maxMismatch,
				overlap_pr1 = overlap_pr1,
				overlap_pr2 = overlap_pr2,
				read_maxlength = read_maxlength,
				pairread_minlength = pairread_minlength,
				merge_minlength = merge_minlength,
				joined_threshold = joined_threshold,
				contamination_threshold = contamination_threshold,
				verbose = verbose
		}
	}

	if (AMPSEQ) {
		call mixed_reads_ampseq_process {
			input:
				type_of_reads = type_of_reads,
				path_to_fq = path_to_fq,
				path_to_flist = path_to_flist,
				pattern_fw = pattern_fw,
				pattern_rv = pattern_rv,
				pr1 = pr1,
				pr2 = pr2,
				Class = Class,
				maxEE = maxEE,
				trimRight = trimRight,
				minLen = minLen,
				truncQ = truncQ,
				matchIDs = matchIDs,
				max_consist = max_consist,
				omegaA = omegaA,
				saveRdata = saveRdata,
				justConcatenate = justConcatenate,
				maxMismatch = maxMismatch,
				overlap_pr1 = overlap_pr1,
				overlap_pr2 = overlap_pr2,
				path_to_snv = path_to_snv,
				no_ref = no_ref,
				reference = reference,
				adjust_mode = adjust_mode,
				reference2 = reference2,
				strain = strain,
				strain2 = strain2,
				polyN = polyN,
				min_reads = min_reads,
				min_samples = min_samples,
				max_snv_dist = max_snv_dist,
				max_indel_dist = max_indel_dist,
				include_failed = include_failed,
				exclude_bimeras = exclude_bimeras,
				amp_mask = amp_mask
		}
	}

	output {
		File? ASVBimeras_f = mixed_reads_ampseq_process.ASVBimeras
		File? CIGARVariants_Bfilter_f = mixed_reads_ampseq_process.CIGARVariants_Bfilter
		File? ASV_to_CIGAR_f = mixed_reads_ampseq_process.ASV_to_CIGAR
		File? seqtab_f = mixed_reads_ampseq_process.seqtab
		File? ASVTable_f = mixed_reads_ampseq_process.ASVTable
		File? ASVSeqs_f = mixed_reads_ampseq_process.ASVSeqs
		File? rawfilelist_f = inline_barcodes_process.rawfilelist
		File? missing_files_f = inline_barcodes_process.missing_files
		File? merge_tar_f = inline_barcodes_process.merge_tar
		File? dada2_ci_tar_f = inline_barcodes_process.dada2_ci_tar
	}
}

task mixed_reads_ampseq_process {
	input {
		String type_of_reads
		String path_to_fq 
		File path_to_flist
		String pattern_fw = "*_L001_R1_001.fastq.gz"
		String pattern_rv = "*_L001_R2_001.fastq.gz"
		File pr1
		File pr2
		String Class = "parasite"
		String maxEE = "5,5"
		String trimRight = "0,0"
		Int minLen = 30
		String truncQ = "5,5"
		String matchIDs = "0"
		Int max_consist = 10
		Float omegaA = 0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
		String saveRdata = ""
		Int justConcatenate = 0
		Int maxMismatch = 0
		File overlap_pr1
		File overlap_pr2
		File path_to_snv
		String no_ref = 'False'
		File reference
		String adjust_mode = "absolute"
		File reference2
		String strain = "3D7"
		String strain2 = "DD2"
		String polyN = "5"
		String min_reads = "0"
		String min_samples = "0"
		String max_snv_dist = "-1"
		String max_indel_dist = "-1"
		String include_failed = "False"
		String exclude_bimeras = "False"
		String amp_mask = "None"
	}

	Map[String, String] in_map = {
		"path_to_fq": "fq_dir",
		"path_to_flist": sub(path_to_flist, "gs://", "/cromwell_root/"),
		"pattern_fw": pattern_fw,
		"pattern_rv": pattern_rv,
		"pr1": sub(pr1, "gs://", "/cromwell_root/"),
		"pr2": sub(pr2, "gs://", "/cromwell_root/"),
		"Class": Class,
		"maxEE": maxEE,
		"trimRight": trimRight,
		"minLen": minLen,
		"truncQ": truncQ,
		"matchIDs": matchIDs,
		"max_consist": max_consist,
		"omegaA": omegaA,
		"saveRdata": saveRdata,
		"justConcatenate": justConcatenate,
		"maxMismatch": maxMismatch,
		"overlap_pr1" : sub(overlap_pr1, "gs://", "/cromwell_root/"),
		"overlap_pr2" : sub(overlap_pr2, "gs://", "/cromwell_root/"),
		"path_to_snv": sub(path_to_snv, "gs://", "/cromwell_root/"),
		"no_ref": no_ref,
		"reference": sub(reference, "gs://", "/cromwell_root/"),
		"adjust_mode": adjust_mode,
		"reference2": sub(reference2, "gs://", "/cromwell_root/"),
		"strain": strain,
		"strain2": strain2,
		"polyN": polyN,
		"min_reads": min_reads,
		"min_samples": min_samples,
		"max_snv_dist": max_snv_dist,
		"max_indel_dist": max_indel_dist,
		"include_failed": include_failed,
		"exclude_bimeras": exclude_bimeras,
		"amp_mask": amp_mask
	}
	File config_json = write_json(in_map)

	command <<<
	set -euxo pipefail
	#set -x
	mkdir fq_dir

	gsutil ls ~{path_to_fq}
	gsutil -m cp -r ~{path_to_fq}* fq_dir/
	find . -type f
	python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --~{type_of_reads} --terra --meta --repo --adaptor_removal --primer_removal --dada2 --postproc_dada2 --asv_to_cigar
	find . -type f
	>>>
	output {
		File? ASVBimeras = "Results/ASVBimeras.txt"
		File? CIGARVariants_Bfilter = "Results/CIGARVariants_Bfilter.out.tsv"
		File? ASV_to_CIGAR = "Results/ASV_to_CIGAR/ASV_to_CIGAR.out.txt"
		File? seqtab = "Results/seqtab.tsv"
		File? ASVTable = "Results/PostProc_DADA2/ASVTable.txt"
		File? ASVSeqs = "Results/PostProc_DADA2/ASVSeqs.fasta"
	}
	runtime {
		cpu: 1
		memory: "15 GiB"
		disks: "local-disk 10 HDD"
		bootDiskSizeGb: 10
		preemptible: 3
		maxRetries: 1
		docker: 'jorgeamaya/mixed_reads_ampseq'
	}
}

task inline_barcodes_process {
	input {
		String type_of_reads
		String path_to_fq 
		File path_to_flist
		String pattern_fw = "*_L001_R1_001.fastq.gz"
		String pattern_rv = "*_L001_R2_001.fastq.gz"
		File pr1
		File pr2
		String Class = "parasite"
		String maxEE = "5,5"
		String trimRight = "0,0"
		Int minLen = 30
		String truncQ = "5,5"
		String matchIDs = "0"
		Int max_consist = 10
		Float omegaA = 0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
		String saveRdata = ""
		Int justConcatenate = 0
		Int maxMismatch = 0
		File overlap_pr1
		File overlap_pr2
		Int read_maxlength = 200
		Int pairread_minlength = 100
		Int merge_minlength = 100
		Int joined_threshold = 1000
		Float contamination_threshold = 0.5
		String verbose = "False"
	}

	Map[String, String] in_map = {
		"path_to_fq": "fq_dir",
		"path_to_flist": sub(path_to_flist, "gs://", "/cromwell_root/"),
		"pattern_fw": pattern_fw,
		"pattern_rv": pattern_rv,
		"pr1": sub(pr1, "gs://", "/cromwell_root/"),
		"pr2": sub(pr2, "gs://", "/cromwell_root/"),
		"Class": Class,
		"maxEE": maxEE,
		"trimRight": trimRight,
		"minLen": minLen,
		"truncQ": truncQ,
		"matchIDs": matchIDs,
		"max_consist": max_consist,
		"omegaA": omegaA,
		"saveRdata": saveRdata,
		"justConcatenate": justConcatenate,
		"maxMismatch": maxMismatch,
		"overlap_pr1" : sub(overlap_pr1, "gs://", "/cromwell_root/"),
		"overlap_pr2" : sub(overlap_pr2, "gs://", "/cromwell_root/"),
		"read_maxlength": read_maxlength,
		"pairread_minlength": pairread_minlength,
		"merge_minlength": merge_minlength,
		"joined_threshold": joined_threshold,
		"contamination_threshold": contamination_threshold,
		"verbose": verbose
	}
	File config_json = write_json(in_map)
	command <<<
	set -euxo pipefail
	#set -x
	mkdir fq_dir

	gsutil ls ~{path_to_fq}
	gsutil -m cp -r ~{path_to_fq}* fq_dir/

	if [ "~{type_of_reads}" == "overlap_reads" ]; then

		python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --~{type_of_reads} --terra --meta --repo --adaptor_removal --bbmerge --bbmerge_report
		#Rscript /Code/BBMerge.R Report/Merge/ Report/

		ls Report/Merge/
		#Rscript /Code/Contamination.R Report/Merge/ Report/ ~{path_to_flist} ~{joined_threshold} ~{contamination_threshold}
		tar -czvf Merge.tar.gz Results/Merge
		find . -type f	
	else
		python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --~{type_of_reads} --terra --meta --repo --adaptor_removal --dada2_contamination
		ls Report/DADA2_Contamination
		tar -cvzf DADA2_Contamination.tar.gz Report/DADA2_Contamination 
		find . -type f
	fi
	>>>
	output {
		File? rawfilelist = "Results/Fq_metadata/rawfilelist.tsv"
		File? missing_files = "Results/missing_files.tsv" 
		File? merge_tar = "Merge.tar.gz"
		File? dada2_ci_tar = "DADA2_Contamination.tar.gz"
		#File bbmergefields = "Report/Merge/bbmergefields.tsv"
		#File BBmerge_performance_absolute_report = "Report/BBmerge_performance_absolute_report.svg"
		#File BBmerge_performance_percentage_report = "Report/BBmerge_performance_percentage_report.svg"
		#File BBmerge_performace_absolute_discarded = "Report/BBmerge_performace_absolute_discarded.svg"	
		#File Barcode_report_abs = "Report/Barcode_report_abs.svg"
		#File Barcode_report_per = "Report/Barcode_report_per.svg"
		#File Insert_size = "Report/Insert_size.png"
		#File Match_report_abs = "Report/Match_report_abs.svg"
		#File Match_report_per = "Report/Match_report_per.svg"
		#File barcodes_report_bbmerge = "Report/barcodes_report_bbmerge.tsv"
		#File hamming_distances_forward = "Report/hamming_forward.tsv"
		#File hamming_distances_reverse = "Report/hamming_reverse.tsv"	
	}
	runtime {
		cpu: 1
		memory: "15 GiB"
		disks: "local-disk 10 HDD"
		bootDiskSizeGb: 10
		preemptible: 3
		maxRetries: 1
		docker: 'jorgeamaya/mixed_reads_ampseq'
	}
}
