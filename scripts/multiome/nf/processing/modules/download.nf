process downloadSRA {
	cpus = params.threads_med
	//Some basic qc for the GEX fastq files - GC%, read #, duplicated sequences etc.
	
	container "nfcore/cellranger:7.1.0"
	
	tag "Download SRA"

	input:
    		val(accession)
    		val(ID)
    	
	output:
        	tuple env(ID), path("fastqs/"), emit: fastqs

	script:
	"""
	#Remove brackets nextflow adds to values
	accession=\$(echo $accession | sed 's/[][]//g')
	ID=\$(echo $ID | sed 's/[][]//g')
	
	#Get file from one of two possible links
	wget https://sra-pub-src-2.s3.amazonaws.com/\${accession}/\${ID}.bam.1 || wget https://sra-pub-src-1.s3.amazonaws.com/\${accession}/\${ID}.bam.1 
	
	#Convert BAM back to fastqs
	cellranger bamtofastq --nthreads=${params.threads_high} \${ID}.bam.1 fastqs/
	
	#Remove BAM
	rm \${ID}.bam.1
	
	#Rename fastqs to ID name
	for f in \$(ls fastqs/*/*.fastq.gz); do mv -v "\$f" "\${f/bamtofastq/\${ID}}"; done;
    
	"""
}

process fasterqDump {
	cpus = params.threads_high
	//Some basic qc for the GEX fastq files - GC%, read #, duplicated sequences etc.
	
	container "jjjermiah/sra-tools:0.1"
	
	tag "fasterq Dump on $ID"

	input:
    		tuple val(ID), val(accession), val(lane)
    		
    	
	output:
        	tuple val(ID), path("*.fastq.gz"), emit: fastqs

	script:
	"""
	
	#Get files
	fasterq-dump --split-files --threads ${params.threads_high} --include-technical ${accession}
	
	###Compress and format header to format for CellRanger
	gzip *
	mv ${accession}_1.fastq.gz ${accession}_S4_L00${lane}_I1_001.fastq.gz
	mv ${accession}_2.fastq.gz ${accession}_S4_L00${lane}_R1_001.fastq.gz
	mv ${accession}_3.fastq.gz ${accession}_S4_L00${lane}_R2_001.fastq.gz
    
	"""
}
