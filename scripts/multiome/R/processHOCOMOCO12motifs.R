#Download necessary files
R.utils::downloadFile("https://hocomoco12.autosome.org/final_bundle/hocomoco12/H12INVIVO/H12INVIVO_pwm.tar.gz", paste0(objects,"motifs/H12INVIVO_pwm.tar.gz"), overwrite=T)
R.utils::gunzip(paste0(objects,"motifs/H12INVIVO_pwm.tar.gz"), overwrite=T)
untar(paste0(objects,"motifs/H12INVIVO_pwm.tar"), exdir=paste0(objects,"motifs/"))
R.utils::downloadFile("https://hocomoco12.autosome.org/final_bundle/hocomoco12/tf_masterlist.tsv", paste0(objects,"motifs/tf_masterlist.tsv"), overwrite=T)

#Read in files
mouseTFs = read.delim(paste0(objects,"motifs/tf_masterlist.tsv"))
humanTFs = mouseTFs[grepl("HUMAN$", mouseTFs$curated.uniprot_id),]
mouseTFs = mouseTFs[grepl("MOUSE$", mouseTFs$curated.uniprot_id),]

#subset human TFs to only those with a mouse orthologue
humanTFs = humanTFs[humanTFs$tfclass.id %in% mouseTFs$tfclass.id,]


#Remove suffix from TF names
humanTFs$curated.uniprot_id = str_extract(humanTFs$curated.uniprot_id, "[^_]+")

#Load in mouse pwms and match the TF names to them
tfFiles = list.files(path = paste0(objects, "motifs/pwm/"))
humanTFs = humanTFs[humanTFs$curated.uniprot_id %in% str_extract(tfFiles, "[^.]+"),]
tfFiles = tfFiles[str_extract(tfFiles, "[^.]+") %in% humanTFs$curated.uniprot_id]
humanTFs = humanTFs[match(str_extract(tfFiles, "[^.]+"), humanTFs$curated.uniprot_id),]
mouseTFs = mouseTFs[match(humanTFs$tfclass.id, mouseTFs$tfclass.id),]
mouseTFs$auto.gene_symbol = make.unique(mouseTFs$auto.gene_symbol)
tfFiles = setNames(tfFiles, mouseTFs$auto.gene_symbol)
humanTFs = humanTFs$curated.uniprot_id

pwmList = purrr::map(tfFiles, function(file){
                                name = str_extract(file, "[^.]+")
                                
                                Matrix <- t(readr::read_table(paste0(objects, "motifs/pwm/",file), 
                                                              col_names = FALSE, skip = 1, progress = F, col_types = cols()))
                                dimnames(Matrix) = list(c("A", "C", "G", "T"))
                                mouseGene = unique(str_extract(mouseTFs[humanTFs == name,]$auto.gene_symbol, "[^.]+"))
                                tfFamily = unique(mouseTFs[humanTFs == name,]$tfclass.family)
                                pfm <- TFBSTools::PWMatrix(ID=mouseGene, name=mouseGene, matrixClass=tfFamily, 
                                                           strand="+", bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                                                           profileMatrix=Matrix)
                                
                              })

pwmList <- do.call(TFBSTools::PWMatrixList, pwmList)
saveRDS(pwmList, file = paste0(objects, "motifs/motifPWMlist.rds"))
