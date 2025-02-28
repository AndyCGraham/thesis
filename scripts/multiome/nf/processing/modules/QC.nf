//Library type
params.library = ""

process FastQC {
	cpus = params.threads_med
	//Some basic qc for the GEX fastq files - GC%, read #, duplicated sequences etc.
	
	container "staphb/fastqc:0.12.1"
	
	tag "FASTQC on $sample_id"
	publishDir "$params.outdir/fastqc/$params.library/", mode: 'copy'

	input:
    		tuple val(sample_id), file(readPairs)
    	
	output:
        	tuple val(sample_id), path("*.{html,zip}"), emit: Results

	script:
	"""
    
	#Rename files with library type if specified
	for f in * ; do mv -- "\$f" "${params.library}_\$f" ; done
	
	find "./" -name "*.fastq.gz" | xargs -n 1 fastqc -t ${params.threads_med} 
	
    
	"""
}


process MultiQC {
	//Concatenates fastqc results for each GEX libraries
	
	container "staphb/multiqc:1.8"

	tag "MultiQC"
	publishDir "$params.outdir/fastqc/", mode: 'copy'
	
	input:
		path "*"
	
	output:
        	file "*multiqc_report.html"
    		file "*_data"

	script:
	"""
	
	##concatenate with multiqc
	multiqc .
	
	"""
}
