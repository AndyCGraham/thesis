#!/usr/bin/bash 

#Subset sample BAMs to prefiltered cells

threads=$1

#ClustersToRemove=("DGmossy" "Subiculum")

mkdir data
mkdir data/bams

export threads=$threads

# Pull dictys container
singularity pull --name containers/dictys.sif docker://andycgraham/dictys:1.1.0

for bam in $(ls data/Multiome/CellRangerARC/*/outs/atac_possorted_bam.bam); do
	#Get sample name from path
	sample="$(dirname "$bam")"       # ...get parent
	sample="$(dirname "$sample")"       # ...get grandparent
	sample="$(basename "$sample")"       # ...then remove everything before the last /
	echo $sample
	
	#Subset BAM to only include pre-filtered cells
	singularity exec containers/dictys.sif subset-bam_linux --bam $bam --cell-barcodes samples/$sample/barcodes.tsv --out-bam samples/$sample/bam.bam --cores ${threads}
	
	#Split into per cell BAMs
	singularity exec containers/dictys.sif dictys_helper split_bam.sh samples/$sample/bam.bam samples/$sample/bams/ --buffer_size 50000 --section "CB:Z:"
	
	#Iterate through per cell BAMs, adding sample prefixes to the old barcode, to match scRNA-seq count matrix
	for cell in $(ls samples/$sample/bams/); do
		BC="${cell%.*}"
		newBC=${sample}_${BC}
		singularity exec containers/dictys.sif samtools view -@ ${threads} -h samples/$sample/bams/$cell  | sed "s/${BC}/${newBC}/" | samtools view -@ ${threads} -b > samples/$sample/bams/${newBC}.bam
		rm samples/$sample/bams/$cell
	done
	
	mv -v samples/$sample/bams/* data/bams
	rm -rf samples/$sample/bams/
done

#Convert GEX counts .mtx file to .tsv file for Dictys
singularity exec containers/dictys.sif dictys_helper expression_mtx.py GEXcounts data/expression.tsv.gz

#Make file of cell subsets
subsets="$(tail -n +2 clusters.csv | awk -F , '{print $2}' | sort -u)"
echo "$subsets" | awk '{print $1}' > data/subsets.txt
for subset in $subsets; do

	if [[ "${ClustersToRemove[@]}" =~ ${subset} ]]; then # Exclude clusters which are too small 
		
		echo "excluding $subset"
		grep -v "${subset}" subsets.txt > subsetstmp.txt; mv subsetstmp.txt subsets.txt #Remove small clusters from subsets file
		grep ",$subset"'$' clusters.csv | awk -F , '{print "data/bams/" $1 ".bam"}'| xargs rm #Remove their cell bam files
	
	else 
	
		mkdir -p "data/subsets/$subset"
		grep ",$subset"'$' clusters.csv | awk -F , '{print $1}' > "data/subsets/$subset/names_rna.txt"
		# RNA and ATAC barcodes are the same for joint quantifications
		cp "data/subsets/$subset/names_rna.txt" "data/subsets/$subset/names_atac.txt"
		#sed -i 's/-1/.1/g' data/subsets/$subset/names_rna.txt
		
	
	fi
	
done
#rm clusters.csv

#Move precomputed peaks here
mkdir data/peaks
for BED in ../../multiomeSeurat/clusters/ExNeu/objects/peaks/broadClusters/*.narrowPeak; do
	directory="$(basename "$BED")"
	directory="${directory%_peaks.narrowPeak*}"
	
	if [[ "${ClustersToRemove[@]}" =~ ${directory} ]]; then # Exclude clusters which are too small 
	
		echo "excluding ${directory}"
		
	else #Otherwise move their peaks to dictys directory
	
		directory=${directory}
		mkdir data/peaks/$directory
		cp $BED data/peaks/$directory/peaks.bed
		
	fi
done

#Download mouse HOMOCOCO motifs
wget -q -o /dev/null -O - 'https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/MOUSE/mono/HOCOMOCOv11_full_MOUSE_mono_homer_format_0.0001.motif' | awk -F "\t" 'BEGIN { OFS = "\t"} {if (substr($1,1,1) == ">") $2=substr($2,1,1)tolower(substr($2,2)); print}' > data/motifs.motif

wget -q -o /dev/null -O - 'https://hocomoco12.autosome.org/final_bundle/hocomoco12/H12CORE/formatted_motifs/H12CORE_homer_format_0.001.motif' | awk -F "\t" 'BEGIN { OFS = "\t"} {if (substr($1,1,1) == ">") $2=substr($2,1,1)tolower(substr($2,2)); print}'  > data/motifs.motif

sed -i 's/.h12core/_mouse.h12core/g' data/motifs.motif

#Get genome
singularity exec containers/dictys.sif dictys_helper genome_homer.sh mm10 data/genome

#Change permissions
chmod -R 755 data/genome

#Get mm10 gene locations bed
wget -q -o /dev/null -O data/gene.gtf.gz https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
gunzip data/gene.gtf.gz
singularity exec containers/dictys.sif dictys_helper gene_gtf.sh data/gene.gtf data/gene.bed
rm data/gene.gtf

#Make Makefile
mkdir makefiles
cd makefiles
singularity exec ../containers/dictys.sif dictys_helper makefile_template.sh common.mk config.mk env_none.mk static.mk

cd ../

#Make Check Makefiles
singularity exec containers/dictys.sif dictys_helper makefile_check.py 

