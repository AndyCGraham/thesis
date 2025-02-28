// Import Required Modules
include { FastQC; MultiQC } from './modules/QC.nf' addParams(library: 'GEX')
include { FastQC as FastQCATAC; MultiQC as MultiQCATAC } from './modules/QC.nf' addParams(library: 'ATAC')
include { CellRangerARC; CellBender } from './modules/CellRanger.nf'


workflow {
	
	// First we get all paired fastq files (R1 + R2)
	// then put all pairs from the same sample (over different lanes) in the same list
	GEX_read_pair = channel.fromFilePairs( params.GEXreadPairs )
	.groupTuple()
	.map {
	[ it[0] , it.flatten().tail() ]
    	}
    	
    	// Now for the ATAC library
    	ATAC_read_pair = channel.fromFilePairs( params.ATACreadPairs )
	.groupTuple()
	.map {
	[ it[0] , it.flatten().tail() ]
    	}
    	
    	// Get path of reference genome
    	Ref = Channel.value( file(params.Ref) )
    	
    	// Get paths of directories holding fastqs
	GEXpath = Channel.fromPath( file(params.GEXreadPairs).parent, type: 'dir' )
	.unique()
	.map {
	[ it.baseName , it ]
    	}
    	
	libraryPaths = Channel.fromPath( file(params.ATACreadPairs).parent, type: 'dir' )
	.unique()
	.map {
	[ it.baseName , it ]
    	}
    	.combine(GEXpath, by: 0)
    	
    	// Run QC on all fastqs
	FastQC( GEX_read_pair ) 
	FastQCATAC( ATAC_read_pair ) 
	FastQCATAC.out.Results.collect()
	
	//Combine QC results for MultiQC
	QCresults = FastQC.out.Results.combine(FastQCATAC.out.Results, by: 0)
	.map {
	it.flatten().tail()
	}
	.collect()
	
	MultiQC( QCresults )
	
	// Run cellranger ARC
	CellRangerARC( libraryPaths , libraryPaths, Ref )
	
	// Use cellbender to correct for ambient ARN contamination
	CellBender( CellRangerARC.out.CellRangerFM )
  
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
