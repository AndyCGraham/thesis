process CellRangerARC {
	cpus = params.threads_high
	// Automatically make Library file per sample to run cellranger arc, then
	// Aligns + quantifies RNA/ATAC counts with cellranger-ARC
	
	container "andycgraham/cellrangerarc:2.0.2"
        //containerOptions ' --mount type=bind,source=${GEXpath} --mount type=bind,source=${ATACpath} '
	
	tag "CellRangerARC on $sample_id"
	publishDir "$params.outdir/CellRangerARC/", mode: 'copy'
	
	input:
		tuple val(sample_id), path(ATACdir, stageAs: "ATACdir"), path(GEXdir, stageAs: "GEXdir")
		tuple val(sample_id), val(ATACpath), val(GEXpath)
		path(Ref)
    	
    	output:
		tuple val(sample_id), file("${sample_id}/outs/raw_feature_bc_matrix.h5"), emit: CellRangerFM
		file("*")

	script:
	"""
	
	echo "fastqs,sample,library_type" > ${sample_id}_libraries.csv
	echo "${GEXpath},${sample_id},Gene Expression" >> ${sample_id}_libraries.csv
	echo "${ATACpath},${sample_id},Chromatin Accessibility" >> ${sample_id}_libraries.csv
	
	##Run cellranger arc
	cellranger-arc count --id=${sample_id} \
                       --reference=${Ref} \
                       --libraries=${sample_id}_libraries.csv \
                       --localcores=${params.threads_high}
	
	"""
}

process CellRanger {
	// Aligns + quantifies RNA counts
	
	tag "CellRanger on $sample_id"
	publishDir "$params.outdir/CellRanger/", mode: 'copy'
	
	container "nfcore/cellranger:7.1.0"
	
	input:
		tuple val(sample_id), file(GEXreads)
		path(Ref)
    	
    	output:
    	        tuple val(sample_id), file("${sample_id}/outs/raw_feature_bc_matrix.h5"), emit: CellRangerFM
		file("*")

	script:
	"""
	
	#Identify runs from the same sample
	samples=\$(ls *.fastq.gz | cut -d_ -f1 | awk '{for (i=1;i<=NF;i++) if (!a[\$i]++) printf("%s,",\$i,FS)}') #Get stirng of fastq prefixes seperated by commas
        samples=\${samples::-1} #Remove last comma
	
	##Run cellranger
	cellranger count --id ${sample_id} \
	--sample \${samples} \
	--transcriptome ${Ref} \
	--fastqs ./ \
	--localcores ${params.threads_high} \
	--no-bam 
	
	"""
}


process CellBender {
	accelerator = params.GPUs
	// This process uses RNAs found in empty droplets to remove counts for RNAs found in the cell/nuclei solution, 
	// from lysed cells etc., which get into droplets and contaminate cells' RNA counts
	
	// Lets do this in cellbender v 0.3.0 docker container
	container "us.gcr.io/broad-dsde-methods/cellbender:0.3.0"

	tag "CellBender on ${sampleID}"
	publishDir "$params.outdir/cellBenderOutput/${sampleID}/", mode: 'copy'
	
	input:
		tuple val(sampleID), path(CellRangerFM)
	
	output:
        	path("CellbenderOutput*")

	script:
	"""
	
	if [ ${params.GPUenabled} == true ]; then
		#Use cellbender to correct for ambient RNA
	  	cellbender remove-background \
		         --input ${CellRangerFM} \
		         --output CellbenderOutput.h5 \
		         --fpr 0.01 \
		         --epochs 150 \
		         --cuda \
		         --exclude-feature-types Peaks
	else
		cellbender remove-background \
		         --input ${CellRangerFM} \
		         --output CellbenderOutput.h5 \
		         --fpr 0.01 \
		         --epochs 150 \
		         --cpu-threads ${params.threads_high} \
		         --exclude-feature-types Peaks
	fi
                 
	"""
}


