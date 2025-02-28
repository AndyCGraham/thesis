#!/usr/bin/env nextflow 

workflow {

	// Show help message
	if (params.help) {
	    helpMessage()
	    exit 0
	}
	
	CellSubsets = Channel.fromPath( params.subsetDir ,  type: 'dir' )
	
	expMat = Channel.value( file( params.expMat ) )
	
	bamDir = Channel.value( file( params.bamDir ) )

	genome = Channel.value( file( params.genome ) )
	
	motifs = Channel.value( file( params.motifPath ) )
	
	geneBed = Channel.value( file( params.geneBed ) )
	
	subsetRNA( CellSubsets, expMat)
	combineSubsetBams( CellSubsets , bamDir )
	CopyData( CellSubsets )
	
	peaks = Channel.fromPath( params.peaksDir ,  type: 'dir' )
	.map {
	[ it.baseName , it ]
    	}
    	.combine(combineSubsetBams.out.SubsetBam, by: 0)
    	.view()
	
	Footprinting( peaks ) 
	MotifSearch( Footprinting.out.Footprints.combine( subsetRNA.out.SubsetRNA, by: 0 ), genome , motifs)
	ScoreTFBinding( MotifSearch.out.Motifs )
	tssDist( MotifSearch.out.Motifs.combine( subsetRNA.out.SubsetRNA, by: 0 ) , geneBed) 
	LinkTFtoGenes( tssDist.out.TSSdist.combine( ScoreTFBinding.out.BindingScore, by: 0 ) )
	DetermineLinkedGenes( LinkTFtoGenes.out.TFtoGenes )
	DetermineLinkedGenes.out.TFtoGenes.combine( subsetRNA.out.SubsetRNA, by: 0 ).view()
	ConstructGRN( DetermineLinkedGenes.out.TFtoGenes.combine( subsetRNA.out.SubsetRNA, by: 0 ) )
	RefineGRN( ConstructGRN.out.GRN )
	IndirectNetwork( ConstructGRN.out.GRN )
	
	tmpDir = Channel.fromPath( params.tmpDir ,  type: 'dir' )
	dataDir = Channel.fromPath( params.dataDir ,  type: 'dir' )
	NetworktoFile( tmpDir, dataDir, RefineGRN.out.Ready.collect(), IndirectNetwork.out.Ready.collect(), CopyData.out.Ready.collect() )
  
}


process subsetRNA {

	//Subsets the RNA data separately for cells of each context/cell type
	
	label "with_cpus"
	tag "subsetRNA on ${subset.baseName}"
	publishDir "$params.tmpDir/${subset.baseName}", mode: 'copy'

	input:
    		file(subset)
    		file(expMat)
    	
	output:
		tuple val(subset.baseName), path("exp.tsv.gz"), emit: SubsetRNA

	script:
	"""
    
	python3 -m dictys preproc selects_rna ${expMat} $subset/names_rna.txt exp.tsv.gz
    
	"""
}

process CopyData {

	//Copies some data to another folder

	label "with_cpus"
	tag "CopyData on ${subset.baseName}"
	publishDir "$params.tmpDir/${subset.baseName}", mode: 'copy'
	
	input:
		file(subset)
		
	output:
		path("${subset.baseName}_subsets.txt"), emit: Ready
		path("names_rna.txt")
		path("names_atac.txt")

	script:
	"""
	
	cp -i $subset/* ./ #Move some files for network aggregation
	
	echo "${subset.baseName}" >> ${subset.baseName}_subsets.txt
	
	"""
}


process combineSubsetBams {
	
	//Combines per-cell bam files to a single bam file for each cell subset
	
	label "with_cpus"
	tag "combineSubsetBams on ${subset.baseName}"
	publishDir "$params.tmpDir/${subset.baseName}", mode: 'copy'

	input:
    		file(subset)
    		val(bamDir)
    	
	output:
        	tuple val(subset.baseName), file("reads.bam"), file("reads.bam.bai"), emit: SubsetBam

	script:
	"""
	
	#Create bam file for custom cells list, filter out chrM/chrUn/chrRandom, sort and index
	awk '{printf("%s\\n","'${bamDir}/'"\$1)}' $subset/names_atac.txt > "00-cells.txt"
	sed -e 's/\$/.bam/' -i "00-cells.txt" #add bam suffix
	( samtools view -h -@ "$params.threads" "\$(head -n 1 "00-cells.txt" )" | grep -v '^@HD' | grep -v '^@PG' ; tail -n +2 "00-cells.txt" | while read l; do samtools view -@ "$params.threads" \$l; done ) | awk '\$3!="chrM"' |  grep -v chrUn_ | grep -v GL00 | grep -v -e "random" | samtools view -1 -@ "$params.threads" -o "02-filtered.bam" -

	#filter, sort and index bam file.
	samtools sort -o "reads.bam" -@ "$params.threads" -l 1 02-filtered.bam
	#rm 02-filtered.bam
	samtools index "reads.bam" "reads.bam.bai"
    
	"""
}

process Footprinting {

	//Finds transcription factor footprints with wellington/pyDNase

	label "with_cpus"
	tag "Footprinting on $subset"
	publishDir "$params.tmpDir/$subset", mode: 'copy'
	
	input:
		tuple val(subset), file(peaks), file(BAM), file(BAI)
		
	
	output:
        	tuple val(subset), file("footprints.bed"), emit: Footprints

	script:
	"""
	
	python3 -m dictys chromatin wellington --nth $params.threads $BAM $BAI $peaks/peaks.bed footprints.bed
	
	"""
}

process MotifSearch {

	//Scans for motif occurrences within provided regions (here open chromatin peaks or footprints) with homer

	label "with_cpus"
	tag "MotifSearch on $subset"
	publishDir "$params.tmpDir/$subset", mode: 'copy'
	
	input:
		tuple val(subset), file(footprints), file(ExpressionMatrix)
		path(genome)
		path(motifs)
	
	output:
        	tuple val(subset), file("motifs.bed"), file("wellington.tsv.gz"), file("homer.tsv.gz"), emit: Motifs

	script:
	"""
	
	python3 -m dictys chromatin homer --nth $params.threads footprints.bed ${motifs} "${genome}" $ExpressionMatrix motifs.bed wellington.tsv.gz homer.tsv.gz
	
	"""
}

process ScoreTFBinding {

	//Computes an integrative score for TF binding based on scores from footprint/peak discovery and from homer

	label "with_cpus"
	tag "ScoreTFBinding on $subset"
	publishDir "$params.tmpDir/$subset", mode: 'copy'
	
	input:
		tuple val(subset), file(motifs), file(wellington), file(homer)
	
	output:
        	tuple val(subset), file("binding.tsv.gz"), emit: BindingScore

	script:
	"""
	
	python3 -m dictys chromatin binding $wellington $homer binding.tsv.gz
	
	"""
	
}

process tssDist {

	//Computes the distance between each (open chromatin) region and each gene's 
	//transcription start site to prioritize nearby pairs that are more likely to 
	//have regulatory effects

	label "with_cpus"
	tag "tssDist on $subset"
	publishDir "$params.tmpDir/$subset", mode: 'copy'
	
	input:
		tuple val(subset), file(motifs), path(WellingtonScores), path(HomerScores), path(ExpressionMatrix)
		path(geneBed)
	
	output:
        	tuple val(subset), file("tssdist.tsv.gz"), emit: TSSdist

	script:
	"""
	
	python3 -m dictys chromatin tssdist $ExpressionMatrix $WellingtonScores $geneBed tssdist.tsv.gz
	
	"""
}

process LinkTFtoGenes {

	//Links TFs to its potential target genes via the relation: TF --- motif --- region --- gene.

	label "with_cpus"
	tag "LinkTFtoGenes on $subset"
	publishDir "$params.tmpDir/$subset", mode: 'copy'
	
	input:
		tuple val(subset), file(tssdist), file(binding)
	
	output:
        	tuple val(subset), file("linking.tsv.gz"), emit: TFtoGenes

	script:
	"""
	
	python3 -m dictys chromatin linking $binding $tssdist linking.tsv.gz
	
	"""
}

process DetermineLinkedGenes {

	//Selects the strongest TF-target gene pairs as a TF binding network that constrains the GRN to be inferred
	
	label "with_cpus"
	tag "DetermineLinkedGenes on $subset"
	publishDir "$params.tmpDir/$subset", mode: 'copy'
	
	input:
		tuple val(subset), file(linking)
	
	output:
        	tuple val(subset), file("binlinking.tsv.gz"), emit: TFtoGenes

	script:
	"""
	
	python3 -m dictys  chromatin binlinking $linking binlinking.tsv.gz 20
	
	"""
}


process ConstructGRN {

	//Uses pyro and pytorch to infer the GRN with stochastic process model under the TF binding network constraint
	
	label "with_gpus"
	tag "ConstructGRN on $subset"
	publishDir "$params.tmpDir/$subset", mode: 'copy'
	
	beforeScript 'module load cuda/11.7.0-gcc-13.2.0'
	
	input:
		tuple val(subset), file(binlinking), file(ExpressionMatrix)
	
	output:
        	tuple val(subset), file("net_weight.tsv.gz"), file("net_meanvar.tsv.gz"), file("net_covfactor.tsv.gz"), file("net_loss.tsv.gz"), file("net_stats.tsv.gz"), emit: GRN

	script:
	"""
	
	python3 -m dictys network reconstruct --device cuda:0 --nth $params.GPUs $ExpressionMatrix $binlinking net_weight.tsv.gz net_meanvar.tsv.gz net_covfactor.tsv.gz net_loss.tsv.gz net_stats.tsv.gz
	
	"""
}

process RefineGRN {

	//Normalizes network edges based on the standard deviation of regulator and target genes. This can overcome biases in the estimation of variance of gene expression due to single-cell sparsity

	label "with_cpus"
	tag "RefineGRN on $subset"
	publishDir "$params.tmpDir/$subset", mode: 'copy'
	
	input:
		tuple val(subset), file(net_weight), file(net_meanvar), file(net_covfactor), file(net_loss), file(net_stats)
	
	output:
        	tuple val(subset), file("net_nweight.tsv.gz"), file(net_meanvar), file(net_covfactor), file(net_loss), file(net_stats), emit: GRN
        	path("${subset}_direct.txt"), emit: Ready

	script:
	"""
	
	python3 -m dictys network normalize --nth $params.threads $net_weight $net_meanvar $net_covfactor net_nweight.tsv.gz
	
	echo "$subset" >> ${subset}_direct.txt
		
	"""
}

process IndirectNetwork {

	//Normalizes network edges based on the standard deviation of regulator and target genes. This can overcome biases in the estimation of variance of gene expression due to single-cell sparsity

	label "with_cpus"
	tag "IndirectNetwork on $subset"
	publishDir "$params.tmpDir/$subset", mode: 'copy'
	
	input:
		tuple val(subset), file(net_weight), file(net_meanvar), file(net_covfactor), file(net_loss), file(net_stats)
	
	output:
        	tuple val(subset), file("net_inweight.tsv.gz"), file("net_iweight.tsv.gz"), file(net_covfactor), file(net_loss), file(net_stats), emit: GRN
        	path("${subset}_indirect.txt"), emit: Ready

	script:
	"""
	
	python3 -m dictys network indirect --nth $params.threads --fi_meanvar $net_meanvar $net_weight $net_covfactor net_iweight.tsv.gz
	
	python3 -m dictys network normalize --nth $params.threads net_iweight.tsv.gz $net_meanvar $net_covfactor net_inweight.tsv.gz
	
	echo "$subset" >> ${subset}_indirect.txt
	
	"""
}


process NetworktoFile {

	//Aggregates all inferred networks to a single output file

	label "with_cpus"
	tag "NetworktoFile"
	publishDir "$params.outDir", mode: 'copy'
	
	input:
		path(dataDir)
		path(tmpDir)
		path(ID_direct)
		path(ID_indirect)
		path(Subsets)
	
	output:
        	file("static.h5")

	script:
	"""
	python3 -m dictys network tofile data/ tmp_static data/subsets.txt static.h5
	
	"""
}

def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow -C nextflowDictys.config run DictysGPU 

       Optional arguments:
        --threads				Number of CPU per subset
        --GPUs					Number of GPU per subset
	--tmpDir 				Directory to store temporary files
	--subsetDir 				Grobbed path to all subset directories
	--outdir				Directory for final network h5 file
	--datadir				Directory containing data for processing
	--peaksDir				Directory containing predetected peak bed files for each subset
	--motifPath				Path of motifs in homer format
	--genomePath				Path of Homer compatible reference genome directory
	--memory				Memory allocation per job
	--time					Job time-limit
        --help                         		Get this usage statement.
        """
}
	
workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}

