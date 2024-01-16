version 1.0

workflow mixed_reads_ampseq {
	input {	
		#General commands
		String type_of_reads
		String path_to_fq 
		File path_to_flist
		String pattern_fw = "*_L001_R1_001.fastq.gz"
		String pattern_rv = "*_L001_R2_001.fastq.gz"

		#Command for the decontamination pipeline
		Int read_maxlength = 200
		Int pairread_minlength = 100
		Int merge_minlength = 100
		Int joined_threshold = 1000
		Float contamination_threshold = 0.5
		String verbose = "False"
	}

	call inline_barcodes_process {
		input: 	
			type_of_reads = type_of_reads,
			path_to_fq = path_to_fq,
			path_to_flist = path_to_flist,
			pattern_fw = pattern_fw,
			pattern_rv = pattern_rv,
			read_maxlength = read_maxlength,
			pairread_minlength = pairread_minlength,
			merge_minlength = merge_minlength,
			joined_threshold = joined_threshold,
			contamination_threshold = contamination_threshold,
			verbose = verbose
	}

	output {
		File rawfilelist_f = inline_barcodes_process.rawfilelist
		File missing_files_f = inline_barcodes_process.missing_files
		File merge_tar_f = inline_barcodes_process.merge_tar
		File bbmergefields_f = inline_barcodes_process.bbmergefields
		File BBmerge_performance_absolute_report_f = inline_barcodes_process.BBmerge_performance_absolute_report
		File BBmerge_performance_percentage_report_f = inline_barcodes_process.BBmerge_performance_percentage_report
		File BBmerge_performace_absolute_discarded_f = inline_barcodes_process.BBmerge_performace_absolute_discarded
		File Barcode_report_abs_f = inline_barcodes_process.Barcode_report_abs
		File Barcode_report_per_f = inline_barcodes_process.Barcode_report_per
		File Insert_size_f = inline_barcodes_process.Insert_size
		File Match_report_abs_f = inline_barcodes_process.Match_report_abs
		File Match_report_per_f = inline_barcodes_process.Match_report_per
		File barcodes_report_bbmerge_f = inline_barcodes_process.barcodes_report_bbmerge
		File hamming_distances_forward_f = inline_barcodes_process.hamming_distances_forward
		File hamming_distances_reverse_f = inline_barcodes_process.hamming_distances_reverse
	}
}

task inline_barcodes_process {
	input {
		String type_of_reads
		String path_to_fq 
		File path_to_flist
		String pattern_fw = "*_L001_R1_001.fastq.gz"
		String pattern_rv = "*_L001_R2_001.fastq.gz"
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

	python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --~{type_of_reads} --meta --repo --adaptor_removal --merge --bbmerge_report
	Rscript /Code/BBMerge.R Report/Merge/ Report/

	ls Report/Merge/
	Rscript /Code/Contamination.R Report/Merge/ Report/ ~{path_to_flist} ~{joined_threshold} ~{contamination_threshold}
	tar -czvf Merge.tar.gz Results/Merge
	find . -type f
	>>>
	output {
		File rawfilelist = "Results/Fq_metadata/rawfilelist.tsv"
		File missing_files = "Results/missing_files.tsv" 
		File merge_tar = "Merge.tar.gz"
		File bbmergefields = "Report/Merge/bbmergefields.tsv"
		File BBmerge_performance_absolute_report = "Report/BBmerge_performance_absolute_report.svg"
		File BBmerge_performance_percentage_report = "Report/BBmerge_performance_percentage_report.svg"
		File BBmerge_performace_absolute_discarded = "Report/BBmerge_performace_absolute_discarded.svg"	
		File Barcode_report_abs = "Report/Barcode_report_abs.svg"
		File Barcode_report_per = "Report/Barcode_report_per.svg"
		File Insert_size = "Report/Insert_size.png"
		File Match_report_abs = "Report/Match_report_abs.svg"
		File Match_report_per = "Report/Match_report_per.svg"
		File barcodes_report_bbmerge = "Report/barcodes_report_bbmerge.tsv"
		File hamming_distances_forward = "Report/hamming_forward.tsv"
		File hamming_distances_reverse = "Report/hamming_reverse.tsv"	
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

