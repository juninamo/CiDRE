# colocalization analysis

suppressMessages(library(GenomicRanges))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(BEDMatrix))
suppressMessages(library(dplyr))
suppressMessages(library(MatrixEQTL))
suppressMessages(library(magrittr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(gassocplot2))
suppressMessages(library(coloc))

setDTthreads(1)
getDTthreads()

run_coloc <- function(eqtl_sumstats, gwas_sumstats){
  tryCatch({
    
    shared = intersect(gwas_sumstats$chr_pos,
                       eqtl_sumstats$chr_pos)
    
    if (length(shared) > 10){
      eqtl_sumstats = eqtl_sumstats %>%
        dplyr::arrange(pos) %>%
        dplyr::filter(chr_pos %in% shared) %>% 
        dplyr::distinct(pos, .keep_all = TRUE) %>% 
        dplyr::mutate_all(funs(ifelse(is.infinite(.),0.1,.))) %>% 
        dplyr::mutate_all(funs(ifelse(is.nan(.),0.1,.))) %>% 
        dplyr::mutate_all(funs(ifelse(is.na(.),0.1,.))) %>% 
        dplyr::mutate(variant_id = as.character(pos))
      gwas_sumstats = gwas_sumstats %>%
        dplyr::arrange(pos) %>%
        dplyr::filter(chr_pos %in% shared) %>% 
        dplyr::distinct(pos, .keep_all = TRUE) %>% 
        dplyr::mutate_all(funs(ifelse(is.infinite(.),0.1,.))) %>% 
        dplyr::mutate_all(funs(ifelse(is.nan(.),0.1,.))) %>% 
        dplyr::mutate_all(funs(ifelse(is.na(.),0.1,.))) %>% 
        dplyr::mutate(variant_id = as.character(pos))
      
      eQTL_dataset = list(beta = eqtl_sumstats$beta,
                          varbeta = eqtl_sumstats$se^2,
                          N = eqtl_sumstats$N, 
                          MAF = eqtl_sumstats$MAF, 
                          type = "quant", 
                          snp = eqtl_sumstats$variant_id)
      gwas_dataset = list(beta = gwas_sumstats$beta,
                          varbeta = gwas_sumstats$se^2,
                          type = "cc", 
                          snp = gwas_sumstats$variant_id,
                          MAF = gwas_sumstats$MAF, 
                          N = gwas_sumstats$N)
      coloc_res = coloc::coloc.abf(dataset1 = eQTL_dataset, dataset2 = gwas_dataset,p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
      res_formatted = dplyr::as_tibble(t(as.data.frame(coloc_res$summary))) %>%
        dplyr::mutate(type = unique(eqtl_sumstats$type),
                      top_qtl_var = dplyr::arrange(eqtl_sumstats, pvalue) %>%
                        dplyr::select(., chr_pos) %>%
                        head(., 1) %>%
                        as.character(), 
                      top_qtl_pval = dplyr::arrange(eqtl_sumstats, pvalue) %>%
                        dplyr::select(., pvalue) %>%
                        head(., 1) %>%
                        as.numeric(), 
                      top_gwas_var = dplyr::arrange(gwas_sumstats, pval) %>%
                        dplyr::select(., chr_pos) %>%
                        head(., 1) %>%
                        as.character(),
                      top_gwas_pval = dplyr::arrange(gwas_sumstats, pval) %>%
                        dplyr::select(., pval) %>%
                        head(., 1) %>%
                        as.numeric() )
      return(res_formatted)
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

zScores <- 
  function (pval, direction = NULL, tails = 2, limit = .Machine$double.xmin){
    if (!is.null(limit)) {
      pval[which(pval < limit)] <- limit
    }
    if (tails == 2) {
      z <- qnorm(pval/2, lower.tail = FALSE)
    }
    else if (tails == 1) {
      z <- qnorm(pval, lower.tail = FALSE)
    }
    else {
      stop("Parameter 'tails' must be set to either 1 or 2.")
    }
    if (!is.null(direction)) {
      z <- z * sign(direction)
    }
    z
  }

# gene information 
gene_coords <- read.table("/data02/home/juninamo/reference/ENSEMBL/GRCh38/GTF/Homo_sapiens.GRCh38.101_gene_annotation_table.txt", header = TRUE)
ref_gr <- GRanges(
  seqnames = paste0("chr",str_split(gene_coords$Chromosome, pattern = ":", simplify = TRUE)[,1]),
  ranges = IRanges(str_split(gene_coords$Chromosome, pattern = ":", simplify = TRUE)[,2], names = gene_coords$GeneSymbol),
  strand = Rle(gene_coords$Strand),
  class = Rle(gene_coords$Class))
ref_df <- as.data.frame(ref_gr, row.names = NULL)
ref_df$gene <- names(ref_gr)
colnames(ref_df) <- gsub("seqnames","chrom",colnames(ref_df))

read2gene <- fread("/data03/inamo/nanopore/RNA/210931_PB29_fastq/PB29_common/GRCh38_all-sr-correct_counts_matrix.tsv") %>%
  .$ids %>%
  stringr::str_split(., pattern = "_", simplify = TRUE) %>%
  as.data.frame() %>%
  dplyr::filter(grepl("ENSG", .$V2)) %>%
  dplyr::mutate(gene_id = stringr::str_split(V2, pattern = "\\.", simplify = TRUE)[,1]) %>%
  dplyr::inner_join(., fread("/data02/home/juninamo/reference/ENSEMBL/GRCh38/GTF/Homo_sapiens.GRCh38.101_gene_annotation_table.txt"), by="gene_id") %>%
  dplyr::select(V1, GeneSymbol) %>%
  dplyr::rename(GENEID = GeneSymbol,query_name = V1)

## EGA
IID_tr <- read.table("/data03/inamo/EGA/sample_Info.txt",header = TRUE)

# covariate file
SUBSET=1:5
POP = "EUR"
for (i in 1:length(SUBSET)){
  for (p in 1:length(POP)){
    tryCatch({
      
      cov = as.data.frame(fread(paste("/data01/EGA/inamo/GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-",SUBSET[i],"_",POP[p],"_PC10.txt",sep=""),header = TRUE))
      cov$SampleID <- str_split(cov$SampleID, pattern = "_", simplify = TRUE)[,ncol(str_split(cov$SampleID, pattern = "_", simplify = TRUE))]
      colnames(cov) <- gsub("SampleID","id",colnames(cov)%>%gsub("^0_","",.))
      assign(paste0("cov_",SUBSET[i],"_",POP[p]), cov)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}



# ImmVar
IID_tr_i <- read.table("/data03/inamo/ImmVar/sample_Info.txt",header = TRUE)

# covariate file
SUBSET <- c("CD4T_Activated","MoDC_unstim","MoDC_FLU","MoDC_IFNb")
POP = "EUR"
for (i in 1:length(SUBSET)){
  for (p in 1:length(POP)){
    tryCatch({
      
      cov = as.data.frame(fread(paste("/data01/ImmVar/inamo/ImmVar/Genotypefiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_",SUBSET[i],"_",POP[p],".hg38_PC10.txt",sep=""),header = TRUE))
      cov$SampleID <- str_split(cov$SampleID, pattern = "_", simplify = TRUE)[,ncol(str_split(cov$SampleID, pattern = "_", simplify = TRUE))]
      colnames(cov) <- gsub("SampleID","id",colnames(cov)%>%gsub("^0_","",.))
      assign(paste0("cov_",SUBSET[i],"_",POP[p]), cov)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}



# DICE
corres = read.table("/data01/DICE/DICE_DNA_RNAseq_correspondance.txt", header = TRUE)

# covariate file
cov = as.data.frame(fread(paste("/data01/DICE/inamo/Genotypefiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup.hg38_PC10.txt",sep=""),header = TRUE)) %>%
  dplyr::mutate(SampleID = str_split(SampleID, pattern = "_", simplify = TRUE)[,ncol(str_split(SampleID, pattern = "_", simplify = TRUE))]) %>%
  magrittr::set_colnames(gsub("SampleID","id",colnames(.)%>%gsub("^0_","",.))) %>%
  tibble::column_to_rownames("id") %>%
  t() %>%
  merge(.,read.table("/data01/DICE/DICE_DNA_RNAseq_correspondance.txt", header = TRUE),by.x="row.names",by.y="IID") %>%
  tibble::column_to_rownames("Run") %>%
  dplyr::select(starts_with("PC")) %>%
  t() %>% as.data.frame() %>%
  dplyr::mutate(id = rownames(.)) %>%
  .[,c("id",as.character(read.table("/data03/inamo/DICE/kallisto/sample_list.txt", header = TRUE)[,1]))]
assign(paste0("cov_DICE"), cov)



# GEUV
IID_tr_g <- 
  read.table("/data03/inamo/GEUV/all_sample.txt", sep="\t") %>%
  .[,c("V1","V3")] %>%
  magrittr::set_colnames(c("IID","stim")) %>%
  magrittr::set_rownames(.$IID) %>%
  .[read.table("/data03/inamo/GEUV/all_sample.txt", sep="\t")[,1],] %>%
  dplyr::mutate(stim = ifelse(grepl("Yoruba", stim), "YRI", "EUR")) %>%
  .[,c("stim","IID")]


# covariate file
POP = "EUR"
for (p in 1:length(POP)){
  tryCatch({
    
    cov = as.data.frame(fread(paste("/data03/inamo/GEUV/VCF/GEUVADIS_GRCh38_",POP[p],"_PC10.txt",sep=""),header = TRUE))
    cov$SampleID <- str_split(cov$SampleID, pattern = "_", simplify = TRUE)[,ncol(str_split(cov$SampleID, pattern = "_", simplify = TRUE))]
    colnames(cov) <- gsub("SampleID","id",colnames(cov)%>%gsub("^0_","",.))
    assign(paste0("cov_GEUV_",POP[p]), cov)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


# normalized expression
isopos_c = fread("/data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/isoform_position.txt") %>% 
  as.data.frame()
psi_c <- fread(paste0("/data03/inamo/sQTL/leafcutter_gencode38/quantile-rank-normalized_psi.txt")) %>%
  dplyr::rename(gid = V1) %>%
  as.data.frame()
peer_factors <- fread(paste0("/data03/inamo/sQTL/leafcutter_gencode38/PEERfactors.txt")) %>%
  as.data.frame()


print (paste0 ("Loading GWAS variants...\n") );
DISEASE="COVID19_HGI_A2"
disease="COVID19_HGI_A2"
pop="EUR"
var_mat = "21_34615210"

print (paste0 ("Loading GWAS variants list...\n") );
var_mat = as.data.frame(fread(paste0("/data02/home/juninamo/GWAS/clumped/",disease,".clumped")))[,c(3,1,4,5)] %>%
  dplyr::filter(SNP == var_mat) %>%
  dplyr::mutate(phenotype = disease) %>%
  dplyr::mutate(chr = paste0("chr",CHR),
                start = BP,
                end = BP) %>%
  dplyr::select(-CHR) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE) %>%
  rtracklayer::liftOver(., rtracklayer::import.chain("/data02/home/juninamo/tools/liftOver/chain_files/hg19ToHg38.over.chain")) %>% 
  as.data.frame() %>%
  dplyr::mutate(SNP_GRCh37 = SNP,
                SNP = paste(seqnames,start,sep="_") %>% gsub("^chr","",.),
                BP = start,
                CHR = seqnames %>% gsub("^chr","",.)) %>%
  dplyr::select(SNP,SNP_GRCh37,CHR,BP,P,phenotype) 


print (paste0 ("number of total independent GWAS-variants;\n") );
nrow(var_mat) %>%
  print()   # number of total independent GWAS-variants
print (paste0 ("number of independent GWAS-variants according to phenotype;\n") );
table(var_mat$phenotype) %>%
  print()  # number of independent GWAS-variants according to phenotype

print (paste0 ("Analyzing ",disease,"...\n") );

if (disease == "COVID19_HGI_A2") {
  var_mat_tmp = var_mat[var_mat$phenotype == disease, ] %>%
    merge(., 
          fread(paste0("/data02/home/juninamo/PRSice/GWAS/",disease,"_ALL_eur_leave_23andme_20210107.b37.txt")), 
          by.x="SNP_GRCh37",
          by.y="chr_pos") %>%
    dplyr::rename(A0 = A2) 
} else {
  var_mat_tmp = var_mat[var_mat$phenotype == disease, ] %>%
    merge(., 
          fread(paste0("/data02/home/juninamo/PRSice/base/alkes/",disease,".txt")), 
          by="SNP_GRCh37",
          by.y="chr_pos") %>%
    dplyr::rename(A0 = A2) 
}
var_mat_tmp = var_mat_tmp %>%
  dplyr::filter(!(chr==6 & pos>25e06 & pos <34e06))

# coloc + Manhattan plot of eQTL and GWAS around GWAS-variants

if (disease == "COVID19_HGI_A2") {
  gwas = fread(paste0("/data02/home/juninamo/PRSice/GWAS/",disease,"_ALL_eur_leave_23andme_20210107.b37.txt")) %>% 
    dplyr::mutate(zscore = qnorm(1-pval),
                  chr = paste0("chr",chr) %>% gsub("chr23","chrX",.),
                  start = pos,
                  end = pos) %>%
    makeGRangesFromDataFrame(., keep.extra.columns = TRUE) %>%
    rtracklayer::liftOver(., rtracklayer::import.chain("/data02/home/juninamo/tools/liftOver/chain_files/hg19ToHg38.over.chain")) %>% 
    as.data.frame() %>%
    dplyr::filter(!(seqnames=="chr6" & start>25e06 & start <34e06)) %>%
    dplyr::filter(seqnames == paste0("chr",var_mat_tmp$CHR) & start>var_mat_tmp$BP-1e06 & start <var_mat_tmp$BP+1e06) %>%
    dplyr::mutate(markers=paste0(seqnames,":",start),
                  chr_pos = gsub("^chr","",markers) %>% gsub(":","_",.),
                  pos = start) %>%
    dplyr::filter(MAF >= 0.05 & MAF <= 0.95 & beta != 0) %>%
    dplyr::select(markers,chr_pos,pos,A2,A1,beta,se,pval,MAF,N,zscore)
  
} else {
  gwas = fread(paste0("/data02/home/juninamo/PRSice/base/alkes/",disease,".txt")) %>% 
    dplyr::mutate(zscore = qnorm(1-pval),
                  chr = paste0("chr",chr) %>% gsub("chr23","chrX",.),
                  start = pos,
                  end = pos) %>%
    makeGRangesFromDataFrame(., keep.extra.columns = TRUE) %>%
    rtracklayer::liftOver(., rtracklayer::import.chain("/data02/home/juninamo/tools/liftOver/chain_files/hg19ToHg38.over.chain")) %>% 
    as.data.frame() %>%
    dplyr::filter(!(seqnames=="chr6" & start>25e06 & start <34e06)) %>%
    dplyr::filter(seqnames == paste0("chr",var_mat_tmp$CHR) & start>var_mat_tmp$BP-1e06 & start <var_mat_tmp$BP+1e06) %>%
    dplyr::mutate(markers=paste0(seqnames,":",start),
                  chr_pos = gsub("^chr","",markers) %>% gsub(":","_",.),
                  pos = start) %>%
    dplyr::filter(MAF >= 0.05 & MAF <= 0.95 & beta != 0) %>%
    dplyr::select(markers,chr_pos,pos,A2,A1,beta,se,pval,MAF,N,zscore)
}

range = 1e05 # range around lead GWAS-variants for investigating genes and variants
lead_variant = var_mat_tmp[1,"SNP"]
lead_rsid = var_mat_tmp[1,"rsid"]
lead_pos = var_mat_tmp[1,"BP"]
chr = var_mat_tmp[1,"CHR"]
lead_A1 = ifelse(var_mat_tmp[1,"beta"] > 0, var_mat_tmp[1,"A1"], var_mat_tmp[1,"A0"])
lead_A0 = ifelse(var_mat_tmp[1,"beta"] > 0, var_mat_tmp[1,"A0"], var_mat_tmp[1,"A1"])

print(lead_variant)
print(lead_rsid)
print(lead_pos)
print(chr)
print(lead_A1)
print(lead_A0)

print (paste0 ("Loading 1000G variants...\n") );
ref <- fread(paste0("/data02/home/juninamo/reference/KGP/AF/all_rs_autosome_EUR_GRCh38_chr",chr,".txt")) %>%
  dplyr::filter(chr_pos %in% gwas$chr_pos)

dir.create(paste0("~/tmp/",disease,"/QTL_peer/splicing/coloc"), showWarnings = F, recursive = T)

# o=1

for ( o in 1:nrow(var_mat_tmp) ) {
  tryCatch({   
    
    range = 1e05 # range around lead GWAS-variants for investigating genes and variants
    
    lead_variant = var_mat_tmp[o,"SNP"]
    lead_rsid = var_mat_tmp[o,"rsid"]
    lead_pos = var_mat_tmp[o,"BP"]
    chr = var_mat_tmp[o,"CHR"]
    lead_A1 = ifelse(var_mat_tmp[o,"beta"] > 0, var_mat_tmp[o,"A1"], var_mat_tmp[o,"A0"])
    lead_A0 = ifelse(var_mat_tmp[o,"beta"] > 0, var_mat_tmp[o,"A0"], var_mat_tmp[o,"A1"])
    
    print(lead_variant)
    print(lead_rsid)
    print(lead_pos)
    print(chr)
    print(lead_A1)
    print(lead_A0)
    
    # EGA
    
    for (stim in 1:5){
      tryCatch({
        
        print(stim)
        variant = fread(paste("/data01/EGA/inamo/GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-",stim,"_",pop,"_chr",chr,".hg38.bim",sep="")) %>%
          dplyr::filter(V1 == chr & V4 > lead_pos-range & V4 < lead_pos+range) %>%
          as.data.frame() %>%
          .[,"V2"]
        snps_i = BEDMatrix(paste("/data01/EGA/inamo/GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-",stim,"_",pop,"_chr",chr,".hg38.bed",sep="")) %>%
          as.matrix() %>%
          magrittr::set_colnames(paste(stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,1],
                                       stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,2],
                                       sep="_")) %>%
          .[, which (colnames(.) %in% variant)] %>%
          as.matrix() %>%
          t()%>%
          data.frame(id = variant, .) %>%
          magrittr::set_colnames(gsub("^X0.","",colnames(.))) %>%
          magrittr::set_colnames(gsub("\\.","@",colnames(.)))
        
        # Position files
        snpspos = data.frame(snp = snps_i$id,
                             chr = paste0("chr",str_split(snps_i$id, pattern = "_", simplify = TRUE)[,1]),
                             pos = as.integer(str_split(snps_i$id, pattern = "_", simplify = TRUE)[,2]))
        
        stimulus = c("Mono_NS","Mono_LPS","Mono_Pam3CSK4","Mono_R848","Mono_IAV")[stim]
        
        cov <- dplyr::bind_rows(eval(parse(text=paste0("cov_",stim,"_",pop))),
                                peer_factors %>%
                                  dplyr::rename(id = PEERfactors) %>%
                                  .[,grepl(paste0("^id$|",stimulus),colnames(.))] %>%
                                  magrittr::set_colnames(gsub(paste0(";",stimulus),"",colnames(.))))
        cov <-  data.frame(id = cov$id, impute::impute.knn(as.matrix(cov[,-1]), k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)$data ) %>%
          magrittr::set_colnames(colnames(cov))
        
        shared_samples = 
          intersect(
            colnames(psi_c) %>% .[grepl(stimulus,.)] %>% stringr::str_split(., pattern = ";", simplify = TRUE) %>% .[,1],
            colnames(snps_i)
          ) %>%
          intersect(
            .,
            colnames(cov)
          )
        
        snps = snps_i[,c("id",shared_samples)]
        expr2 <- psi_c[, c("gid", paste0(shared_samples,";",stimulus)) ] %>%
          magrittr::set_colnames(gsub(paste0(";",stimulus),"",colnames(.))) %>%
          .[ which ( .$gid %in% isopos_c[isopos_c$chr == paste0("chr",chr) & 
                                           isopos_c$s1 > min(snpspos$pos-(5e+05)) & 
                                           isopos_c$s2 < max(snpspos$pos+(5e+05)), "geneid"] ), ]
        cov <- cov[,c("id",shared_samples)]
        
        if ( nrow(expr2) > 0 &&
             all(colnames(snps)[-1] == colnames(expr2)[-1] &&
                 colnames(snps)[-1] == colnames(cov)[-1]) ) {
          
          # print(head(expr2))
          # print(dim(expr2))
          # print(head(snps))
          # print(dim(snps))
          # print(head(cov))
          # print(dim(cov))
          
          # paste("Mean value of expression for gene ",expr2$gid," is ", rowMeans(expr2[, -1]))
          # paste("Standard deviation of expression for gene ", expr2$gid," is ", t(apply(expr2[, -1], 1, sd)))
          # paste("Mean value for SNP ",snps$id," is ", rowMeans(snps[, -1]))
          # paste("Standard deviation for SNP ", snps$id," is ", t(apply(snps[, -1], 1, sd)))
          
          # MatrixQTLでmQTL-SNPsと発現量が相関する遺伝子を抽出
          write.table(snps, paste0("~/tmp/snps_",lead_variant,"_",stim,"_",disease,"_leafcutter_sqtl.txt"), sep = "\t",quote=FALSE, row.names = FALSE)
          write.table(expr2, paste0("~/tmp/exp_",lead_variant,"_",stim,"_",disease,"_leafcutter_sqtl.txt"), sep = "\t",quote=FALSE, row.names = FALSE)
          write.table(cov, paste0("~/tmp/cov_",lead_variant,"_",stim,"_",disease,"_leafcutter_sqtl.txt"), sep = "\t",quote=FALSE, row.names = FALSE)
          
          # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
          useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
          
          # Genotype file name
          SNP_file_name = paste("~/tmp/snps_",lead_variant,"_",stim,"_",disease,"_leafcutter_sqtl.txt", sep="");
          
          # Gene expression file name
          expression_file_name = paste("~/tmp/exp_",lead_variant,"_",stim,"_",disease,"_leafcutter_sqtl.txt", sep="");
          
          # Covariates file name
          # Set to character() for no covariates
          covariates_file_name = paste("~/tmp/cov_",lead_variant,"_",stim,"_",disease,"_leafcutter_sqtl.txt", sep="");
          
          
          # Error covariance matrix
          # Set to numeric() for identity.
          errorCovariance = numeric();
          
          pvOutputThreshold = 1;
          
          # errorCovariance = read.table("Sample_Data/errorCovariance.txt");
          
          # Output file name
          output_file_name = tempfile();
          
          ## Load genotype data
          SNPs = SlicedData$new();
          SNPs$fileDelimiter = "\t";      # the TAB character
          SNPs$fileOmitCharacters = "NA"; # denote missing values;
          SNPs$fileSkipRows = 1;          # one row of column labels
          SNPs$fileSkipColumns = 1;       # one column of row labels
          SNPs$fileSliceSize = 2000;      # read file in slices of 2,000 rows
          SNPs$LoadFile(SNP_file_name);
          
          ## Load gene expression data
          gene = SlicedData$new();
          gene$fileDelimiter = "\t";      # the TAB character
          gene$fileOmitCharacters = "NA"; # denote missing values;
          gene$fileSkipRows = 1;          # one row of column labels
          gene$fileSkipColumns = 1;       # one column of row labels
          gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
          gene$LoadFile(expression_file_name);
          
          ## Load covariates
          cvrt = SlicedData$new();
          cvrt$fileDelimiter = "\t";      # the TAB character
          cvrt$fileOmitCharacters = "NA"; # denote missing values;
          cvrt$fileSkipRows = 1;          # one row of column labels
          cvrt$fileSkipColumns = 1;       # one column of row labels
          if(length(covariates_file_name)>0) {
            cvrt$LoadFile(covariates_file_name);
          }
          
          
          ## Run the analysis
          
          me = Matrix_eQTL_engine(
            snps = SNPs,
            gene = gene,
            cvrt = cvrt,
            output_file_name = output_file_name,
            pvOutputThreshold = pvOutputThreshold,
            useModel = useModel,
            errorCovariance = errorCovariance,
            verbose = FALSE,
            pvalue.hist = TRUE,
            min.pv.by.genesnp = FALSE,
            noFDRsaveMemory = FALSE);
          
          unlink(output_file_name);
          
          
          ## Results:
          cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
          # cat('Detected local eQTLs:', '\n');
          # show(me$cis$eqtls)
          # cat('Detected distant eQTLs:', '\n');
          # show(me$trans$eqtls)
          
          N = read.table(paste("/data01/EGA/inamo/GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-",stim,"_EUR_chr",chr,".hg38.fam",sep="")) %>% nrow() %>% as.integer() 
          
          ## Results:
          
          if (stim == 1){
            manh_data = merge(me$all$eqtls, snpspos, by.x = "snps", by.y = "snp") %>%
              dplyr::mutate(type = paste0(stimulus),
                            N = N)
          } else {
            manh_data = merge(me$all$eqtls, snpspos, by.x = "snps", by.y = "snp") %>%
              dplyr::mutate(type = paste0(stimulus),
                            N = N) %>%
              rbind(.,manh_data)
          }
          
          
        }
        
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    
    
    # DICE
    variant = fread(paste("/data01/DICE/inamo/Genotypefiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_chr",chr,".hg38.bim",sep="")) %>%
      dplyr::filter(V1 == chr & V4 > lead_pos-range & V4 < lead_pos+range) %>%
      as.data.frame() %>%
      .[,"V2"]
    snps_i = BEDMatrix(paste("/data01/DICE/inamo/Genotypefiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_chr",chr,".hg38.bed",sep="")) %>%
      as.matrix() %>%
      magrittr::set_colnames(paste(stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,1],
                                   stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,2],
                                   sep="_")) %>%
      .[, which (colnames(.) %in% variant)] %>%
      as.matrix() %>%
      t()%>%
      data.frame(id = variant, .) %>%
      magrittr::set_colnames(gsub("^X0.","",colnames(.))) %>%
      magrittr::set_colnames(gsub("\\.","-",colnames(.)))
    
    # Position files
    
    snpspos = data.frame(snp = snps_i$id,
                         chr = paste0("chr",str_split(snps_i$id, pattern = "_", simplify = TRUE)[,1]),
                         pos = as.integer(str_split(snps_i$id, pattern = "_", simplify = TRUE)[,2]))
    
    for (stim in c("TFH","TH1","TH17","TH2","TH1-17","TREG_MEMORY","TREG_NAIVE","B_NAIVE","CD4_NAIVE","CD4_N_STIM","CD8_NAIVE","CD8_N_STIM","NONCLASSICAL_MONOCYTES","CLASSICAL_MONOCYTES","NK_CD16POS")){
      tryCatch({
        
        print(stim)
        cov <- dplyr::bind_rows(eval(parse(text=paste0("cov_DICE"))),
                                peer_factors %>%
                                  dplyr::rename(id = PEERfactors) %>%
                                  .[,grepl(paste0("^id$|^SRR"),colnames(.))])
        cov <-  data.frame(id = cov$id, impute::impute.knn(as.matrix(cov[,-1]), k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)$data ) %>%
          magrittr::set_colnames(colnames(cov))
        
        shared_samples = 
          intersect(
            corres[ grepl("_1_RNA_",corres$biospecimen_repository_sample_id) & corres$histological_type == stim , "Run"],
            colnames(psi_c) %>% .[grepl("^SRR",.)]
          )  %>%
          intersect(
            .,
            colnames(cov)
          )
        
        shared_IID = 
          corres %>%
          dplyr::filter(Run %in% shared_samples) %>%
          .$IID %>% as.character() %>%
          intersect(
            .,
            colnames(snps_i)
          )
        shared = 
          corres %>%
          dplyr::filter(Run %in% shared_samples) %>%
          dplyr::filter(IID %in% shared_IID) %>%
          magrittr::set_rownames(.$IID) %>%
          .[shared_IID,]
        
        snps = snps_i[,c("id",shared$IID %>% as.character())] %>%
          magrittr::set_colnames(c("id",shared$Run %>% as.character()))
        expr2 <- psi_c[, c("gid", shared$Run %>% as.character()) ] %>%
          .[ which ( .$gid %in% isopos_c[isopos_c$chr == paste0("chr",chr) & 
                                           isopos_c$s1 > min(snpspos$pos-(5e+05)) & 
                                           isopos_c$s2 < max(snpspos$pos+(5e+05)), "geneid"] ), ]
        
        cov <- cov[,c("id",shared$Run %>% as.character())]
        
        
        if ( nrow(expr2) > 0 &&
             all(colnames(snps)[-1] == colnames(expr2)[-1]) &&
             all(colnames(snps)[-1] == colnames(cov)[-1]) ) {
          
          # print(head(expr2))
          # print(dim(expr2))
          # print(head(snps))
          # print(dim(snps))
          # print(dim(cov))
          
          # paste("Mean value of expression for gene ",expr2$gid," is ", rowMeans(expr2[, -1]))
          # paste("Standard deviation of expression for gene ", expr2$gid," is ", t(apply(expr2[, -1], 1, sd)))
          # paste("Mean value for SNP ",snps$id," is ", rowMeans(snps[, -1]))
          # paste("Standard deviation for SNP ", snps$id," is ", t(apply(snps[, -1], 1, sd)))
          
          # MatrixQTLでmQTL-SNPsと発現量が相関する遺伝子を抽出
          write.table(snps, paste0("~/tmp/snps_dice_",lead_variant,"_",stim,"_",disease,"_leafcutter_sqtl.txt"), sep = "\t",quote=FALSE, row.names = FALSE)
          write.table(expr2, paste0("~/tmp/exp_dice_",lead_variant,"_",stim,"_",disease,"_leafcutter_sqtl.txt"), sep = "\t",quote=FALSE, row.names = FALSE)
          write.table(cov, paste0("~/tmp/cov_dice_",lead_variant,"_",stim,"_",disease,"_leafcutter_sqtl.txt"), sep = "\t",quote=FALSE, row.names = FALSE)
          
          # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
          useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
          
          # Genotype file name
          SNP_file_name = paste("~/tmp/snps_dice_",lead_variant,"_",stim,"_",disease,"_leafcutter_sqtl.txt", sep="");
          
          # Gene expression file name
          expression_file_name = paste("~/tmp/exp_dice_",lead_variant,"_",stim,"_",disease,"_leafcutter_sqtl.txt", sep="");
          
          # Covariates file name
          # Set to character() for no covariates
          covariates_file_name = paste("~/tmp/cov_dice_",lead_variant,"_",stim,"_",disease,"_leafcutter_sqtl.txt", sep="");
          
          
          # Error covariance matrix
          # Set to numeric() for identity.
          errorCovariance = numeric();
          
          pvOutputThreshold = 1;
          
          # errorCovariance = read.table("Sample_Data/errorCovariance.txt");
          
          # Output file name
          output_file_name = tempfile();
          
          ## Load genotype data
          SNPs = SlicedData$new();
          SNPs$fileDelimiter = "\t";      # the TAB character
          SNPs$fileOmitCharacters = "NA"; # denote missing values;
          SNPs$fileSkipRows = 1;          # one row of column labels
          SNPs$fileSkipColumns = 1;       # one column of row labels
          SNPs$fileSliceSize = 2000;      # read file in slices of 2,000 rows
          SNPs$LoadFile(SNP_file_name);
          
          ## Load gene expression data
          gene = SlicedData$new();
          gene$fileDelimiter = "\t";      # the TAB character
          gene$fileOmitCharacters = "NA"; # denote missing values;
          gene$fileSkipRows = 1;          # one row of column labels
          gene$fileSkipColumns = 1;       # one column of row labels
          gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
          gene$LoadFile(expression_file_name);
          
          ## Load covariates
          cvrt = SlicedData$new();
          cvrt$fileDelimiter = "\t";      # the TAB character
          cvrt$fileOmitCharacters = "NA"; # denote missing values;
          cvrt$fileSkipRows = 1;          # one row of column labels
          cvrt$fileSkipColumns = 1;       # one column of row labels
          if(length(covariates_file_name)>0) {
            cvrt$LoadFile(covariates_file_name);
          }
          
          
          ## Run the analysis
          
          me = Matrix_eQTL_engine(
            snps = SNPs,
            gene = gene,
            cvrt = cvrt,
            output_file_name = output_file_name,
            pvOutputThreshold = pvOutputThreshold,
            useModel = useModel,
            errorCovariance = errorCovariance,
            verbose = FALSE,
            pvalue.hist = TRUE,
            min.pv.by.genesnp = FALSE,
            noFDRsaveMemory = FALSE);
          
          unlink(output_file_name);
          
          ## Results:
          cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
          # cat('Detected local eQTLs:', '\n');
          # show(me$cis$eqtls)
          # cat('Detected distant eQTLs:', '\n');
          # show(me$trans$eqtls)
          
          N = read.table(paste("/data01/DICE/inamo/Genotypefiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_chr",chr,".hg38.fam",sep="")) %>% nrow()
          
          ## Results:
          
          manh_data = merge(me$all$eqtls, snpspos, by.x = "snps", by.y = "snp") %>%
            dplyr::mutate(type = stim,
                          N = N ) %>%
            rbind(.,manh_data)
          
          
        }
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    
    
    gene_list = manh_data %>%
      dplyr::group_by(gene,snps) %>%
      dplyr::mutate(count=n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(FDR < 0.05) %>%
      dplyr::select(gene) %>%
      unique() %>%
      as.data.frame() %>%
      .[,"gene"] %>%
      as.character() %>% sort()
    
    for (g in 1:length(gene_list)){
      
      gene_symbol = 
        isopos_c %>%
        dplyr::rename(chrom = chr, start = s1, end = s2) %>%
        dplyr::select(chrom,start,end,geneid) %>%
        unique() %>%
        valr::bed_closest(., ref_df, suffix = c("_sQTL", "_ref")) %>%
        dplyr::filter(geneid_sQTL == gene_list[g]) %>%
        .$gene_ref %>%
        paste0(., collapse = "-")
      
      if (!file.exists(paste0("~/tmp/",disease,"/QTL_peer/splicing/coloc/chr",lead_variant,"-",lead_rsid,"-",lead_A1,"_leafcutter_",gsub(":","-",gene_list[g]),"-",gene_symbol,".txt"))){ 
        
        # Manhattan plot
        manh_data_tmp = manh_data %>%
          dplyr::filter(gene == gene_list[g]) %>%
          dplyr::mutate(zscore = zScores(pvalue, direction=beta, tails=2, limit=.Machine$double.xmin)) %>%
          dplyr::mutate(markers=paste0(.$chr,":",.$pos)) %>%
          dplyr::select(markers,zscore,type) %>%
          tidyr::spread(data = ., key = type, value = zscore) %>%
          merge(.,gwas,by="markers") %>%
          magrittr::set_colnames(gsub("zscore",disease,colnames(.))) %>%
          tibble::column_to_rownames("markers") %>%
          dplyr::select(manh_data$type %>% as.character() %>% unique(), paste0(disease))
        pos_lt = rownames(manh_data_tmp)
        manh_data_tmp = 
          manh_data_tmp %>% 
          dplyr::mutate_all(funs(ifelse(is.infinite(.),0,.))) %>% 
          dplyr::mutate_all(funs(ifelse(is.nan(.),0,.))) %>% 
          dplyr::mutate_all(funs(ifelse(is.na(.),0,.))) %>% 
          magrittr::set_rownames(pos_lt)
        
        dim(manh_data_tmp)
        
        
        # coloc
        
        for (stim in manh_data$type %>% as.character() %>% unique()){
          tryCatch({
            
            print(stim)
            
            if (grepl("LCL", stim)){
              
              coloc_tmp = 
                manh_data %>%
                dplyr::filter(gene == gene_list[g] & type == stim) %>%
                dplyr::mutate(se = dmetar::se.from.p(.$beta, .$pvalue, .$N, effect.size.type = "difference", calculate.g = FALSE)[,"StandardError"] %>% abs()) %>%
                dplyr::left_join(.,
                                 ref %>% dplyr::rename(snps = chr_pos),
                                 by="snps") %>%
                merge(.,
                      fread(paste0("/data03/inamo/GEUV/VCF/GEUVADIS_GRCh38_EUR_chr",chr,".bim")),
                      by.x="snps",by.y="V2") %>%
                dplyr::rename(chr_pos = snps,
                              A1 = V5,
                              A0 = V6) %>%
                dplyr::mutate(MAF = ifelse(is.na(EUR), 0.5, EUR),
                              MAF =  ifelse( (A1 == "A" & (ALT == "A" | ALT == "T")) |
                                               (A1 == "T" & (ALT == "A" | ALT == "T")) | 
                                               (A1 == "C" & (ALT == "C" | ALT == "G")) |
                                               (A1 == "G" & (ALT == "C" | ALT == "G")), 
                                             MAF, 1-MAF))
              assign(paste0("coloc_",gsub("-","_",stim)),coloc_tmp)
              
            } else if (grepl("Mono_", stim)) {
              
              num = grep(stim,c("Mono_NS","Mono_LPS","Mono_Pam3CSK4","Mono_R848","Mono_IAV"))
              coloc_tmp = 
                manh_data %>%
                dplyr::filter(gene == gene_list[g] & type == stim) %>%
                dplyr::mutate(se = dmetar::se.from.p(.$beta, .$pvalue, .$N, effect.size.type = "difference", calculate.g = FALSE)[,"StandardError"] %>% abs()) %>%
                dplyr::left_join(.,
                                 ref %>% dplyr::rename(snps = chr_pos),
                                 by="snps") %>%
                merge(.,
                      fread(paste0("/data01/EGA/inamo/GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-",num,"_EUR_chr",chr,".hg38.bim")),
                      by.x="snps",by.y="V2") %>%
                dplyr::rename(chr_pos = snps,
                              A1 = V5,
                              A0 = V6) %>%
                dplyr::mutate(MAF = ifelse(is.na(EUR), 0.5, EUR),
                              MAF =  ifelse( (A1 == "A" & (ALT == "A" | ALT == "T")) |
                                               (A1 == "T" & (ALT == "A" | ALT == "T")) | 
                                               (A1 == "C" & (ALT == "C" | ALT == "G")) |
                                               (A1 == "G" & (ALT == "C" | ALT == "G")), 
                                             MAF, 1-MAF))
              assign(paste0("coloc_",gsub("-","_",stim)),coloc_tmp)
              
            } else if ( grepl("MoDC_", stim) | grepl("CD4T_", stim) ) {
              
              coloc_tmp = 
                manh_data %>%
                dplyr::filter(gene == gene_list[g] & type == stim) %>%
                dplyr::mutate(se = dmetar::se.from.p(.$beta, .$pvalue, .$N, effect.size.type = "difference", calculate.g = FALSE)[,"StandardError"] %>% abs()) %>%
                dplyr::left_join(.,
                                 ref %>% dplyr::rename(snps = chr_pos),
                                 by="snps") %>%
                merge(.,
                      fread(paste0("/data01/ImmVar/inamo/ImmVar/Genotypefiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_",stim,"_EUR_chr",chr,".hg38.bim")),
                      by.x="snps",by.y="V2") %>%
                dplyr::rename(chr_pos = snps,
                              A1 = V5,
                              A0 = V6) %>%
                dplyr::mutate(MAF = ifelse(is.na(EUR), 0.5, EUR),
                              MAF =  ifelse( (A1 == "A" & (ALT == "A" | ALT == "T")) |
                                               (A1 == "T" & (ALT == "A" | ALT == "T")) | 
                                               (A1 == "C" & (ALT == "C" | ALT == "G")) |
                                               (A1 == "G" & (ALT == "C" | ALT == "G")), 
                                             MAF, 1-MAF))
              assign(paste0("coloc_",gsub("-","_",stim)),coloc_tmp)
              
            } else {
              
              coloc_tmp = 
                manh_data %>%
                dplyr::filter(gene == gene_list[g] & type == stim) %>%
                dplyr::mutate(se = dmetar::se.from.p(.$beta, .$pvalue, .$N, effect.size.type = "difference", calculate.g = FALSE)[,"StandardError"] %>% abs()) %>%
                dplyr::left_join(.,
                                 ref %>% dplyr::rename(snps = chr_pos),
                                 by="snps") %>%
                merge(.,
                      fread(paste0("/data01/DICE/inamo/Genotypefiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_chr",chr,".hg38.bim")),
                      by.x="snps",by.y="V2") %>%
                dplyr::rename(chr_pos = snps,
                              A1 = V5,
                              A0 = V6) %>%
                dplyr::mutate(MAF = ifelse(is.na(EUR), 0.5, EUR),
                              MAF =  ifelse( (A1 == "A" & (ALT == "A" | ALT == "T")) |
                                               (A1 == "T" & (ALT == "A" | ALT == "T")) | 
                                               (A1 == "C" & (ALT == "C" | ALT == "G")) |
                                               (A1 == "G" & (ALT == "C" | ALT == "G")), 
                                             MAF, 1-MAF))
              assign(paste0("coloc_",gsub("-","_",stim)),coloc_tmp)
              
            }
            
          }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        }
        
        eQTL_list = setNames(as.list(c(manh_data$type %>% as.character() %>% unique())), 
                             c(manh_data$type %>% as.character() %>% unique()))
        
        for ( stim in manh_data$type %>% as.character() %>% unique() ){
          tryCatch({
            
            eQTL_list[[stim]] <- eval(parse(text=paste0("coloc_",gsub("-","_",stim))))
            
          }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        }
        
        
        purrr::map_df(eQTL_list, ~run_coloc(., gwas)) %>%
          dplyr::mutate(disease = disease,
                        gene = gene_list[g],
                        gene_symbol = gene_symbol) %>%
          dplyr::arrange(-PP.H4.abf) %>%
          as.data.frame() %>%
          write.table(., paste0("~/tmp/",disease,"/QTL_peer/splicing/coloc/chr",lead_variant,"-",lead_rsid,"-",lead_A1,"_leafcutter_",gsub(":","-",gene_list[g]),"-",gene_symbol,".txt"), sep = "\t",quote=FALSE, row.names = FALSE)
        
      }
    }
    
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
print (paste0 (disease, ", all done!\n") );
lf <- list.files(path = paste0("~/tmp/",disease,"/QTL_peer/splicing/coloc"), pattern = paste0("*_leafcutter_*"), full.names = T)
as.data.frame(do.call(rbind, lapply(lf,fread))) %>%
  dplyr::filter(type == "NONCLASSICAL_MONOCYTES" | type == "CLASSICAL_MONOCYTES" | type == "Mono_NS" ) %>%
  dplyr::filter(gene == "21:33252830:33268394:clu_34140_+") 
as.data.frame(do.call(rbind, lapply(lf,fread))) %>%
  dplyr::arrange(-PP.H4.abf) %>%
  write.table(., paste0("~/coloc.txt"), sep = "\t",quote=FALSE, row.names = FALSE)


# locuscompre

gwas = "COVID19_HGI_A2"
tool = "leafcutter"
lead_variant_gwas = "21_33242905"
sQTL_id = "21:33252830:33268394:clu_34140_+" 
range = 5e05
population = "EUR"
combine = TRUE
legend = TRUE
legend_position = "topleft"
lz_ylab_linebreak = FALSE
genome = "hg38"

chr = sQTL_id %>% stringr::str_split(., pattern = ":", simplify = TRUE) %>% .[,1] %>% paste0("chr",.)
start =  sQTL_id %>% stringr::str_split(., pattern = ":", simplify = TRUE) %>% .[,2] %>% as.integer()
end =  sQTL_id %>% stringr::str_split(., pattern = ":", simplify = TRUE) %>% .[,3] %>% as.integer()

suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(BEDMatrix))
suppressMessages(library(dplyr))
suppressMessages(library(MatrixEQTL))
suppressMessages(library(magrittr))
suppressMessages(library(coloc))
suppressMessages(library(locuscomparer))

peer_factors <- fread(paste0("/data03/inamo/sQTL/leafcutter_gencode38/PEERfactors.txt")) %>%
  as.data.frame()


lead_pos = stringr::str_split(lead_variant_gwas,pattern="_",simplify=TRUE) %>% .[,2] %>% as.integer()
chr = stringr::str_split(lead_variant_gwas,pattern="_",simplify=TRUE) %>% .[,1]

print(lead_pos)
print(chr)

rsid2chrpos = fread(paste0("/data03/inamo/nanopore/RNA/210931_PB29_fastq/shiny/data/dbSNP/rs_ChrPos_",chr,"_MAF0.01.txt")) %>%
  magrittr::set_colnames(c("rsid","snps"))

tmp_snps = tempfile()
tmp_exp = tempfile()
tmp_cov = tempfile()

for (cell in c("NONCLASSICAL_MONOCYTES","CLASSICAL_MONOCYTES","Mono_NS")){
  print(cell)
  # EGA
  if (grepl("Mono_", cell)) {
    print (paste0 ("Loading transcript ratio of EGA...\n") );
    
    stimulus = grep(cell,c("Mono_NS","Mono_LPS","Mono_Pam3CSK4","Mono_R848","Mono_IAV"))
    
    variant = fread(paste("/data01/EGA/inamo/GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-",stimulus,"_EUR_chr",chr,".hg38.bim",sep="")) %>%
      dplyr::filter(V1 == chr & V4 > lead_pos-range & V4 < lead_pos+range) %>%
      as.data.frame() %>%
      .[,"V2"]
    snps_i = BEDMatrix(paste("/data01/EGA/inamo/GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-",stimulus,"_EUR_chr",chr,".hg38.bed",sep="")) %>%
      as.matrix() %>%
      magrittr::set_colnames(paste(stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,1],
                                   stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,2],
                                   sep="_")) %>%
      .[, which (colnames(.) %in% variant)] %>%
      as.matrix() %>%
      t()%>%
      data.frame(id = variant, .) %>%
      magrittr::set_colnames(gsub("^X0.","",colnames(.))) %>%
      magrittr::set_colnames(gsub("\\.","@",colnames(.)))
    
    # Position files
    snpspos = data.frame(snp = snps_i$id,
                         chr = paste0("chr",str_split(snps_i$id, pattern = "_", simplify = TRUE)[,1]),
                         pos = as.integer(str_split(snps_i$id, pattern = "_", simplify = TRUE)[,2]))
    
    cov <- dplyr::bind_rows(fread(paste("/data01/EGA/inamo/GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-",stimulus,"_EUR_PC10.txt",sep=""),header = TRUE) %>%
                              as.data.frame() %>%
                              dplyr::mutate(SampleID = str_split(SampleID, pattern = "_", simplify = TRUE)[,ncol(str_split(SampleID, pattern = "_", simplify = TRUE))]) %>%
                              magrittr::set_colnames(gsub("SampleID","id",colnames(.)%>%gsub("^0_","",.))),
                            peer_factors %>%
                              dplyr::rename(id = PEERfactors) %>%
                              .[,grepl(paste0("^id$|",cell),colnames(.))] %>%
                              magrittr::set_colnames(gsub(paste0(";",cell),"",colnames(.))))
    cov <-  data.frame(id = cov$id, impute::impute.knn(as.matrix(cov[,-1]), k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)$data ) %>%
      magrittr::set_colnames(colnames(cov))
    
    
    tr_ratio =
      fread(paste0("/data03/inamo/sQTL/leafcutter_gencode38/expression_quantile-rank-normalized/isoform_ratio/",tool,"/",sQTL_id,".txt")) %>% 
      dplyr::rename(gid = V1) %>%
      as.data.frame() %>%
      .[,grepl(paste0("^gid|",cell),colnames(.))]
    
    IID_tr <- data.frame(stim = stringr::str_split(colnames(tr_ratio), pattern = ";", simplify = TRUE)[,2][-1], 
                         IID = stringr::str_split(colnames(tr_ratio), pattern = ";", simplify = TRUE)[,1][-1])
    
    shared_samples = 
      intersect(
        colnames(tr_ratio) %>% .[grepl(cell,.)] %>% stringr::str_split(., pattern = ";", simplify = TRUE) %>% .[,1],
        colnames(snps_i)
      ) %>%
      intersect(
        .,
        colnames(cov)
      )
    
    snps = snps_i[,c("id",shared_samples)]
    expr2 <- tr_ratio %>%
      as.data.frame() %>% 
      .[, c("gid", paste0(shared_samples,";",cell)) ] %>%
      magrittr::set_colnames(gsub(paste0(";",cell),"",colnames(.))) %>%
      magrittr::set_colnames(gsub("\\.","@",colnames(.))) 
    cov <- cov[,c("id",shared_samples)]
    
  } else {
    print (paste0 ("Loading transcript ratio of DICE...\n") );
    variant = fread(paste("/data01/DICE/inamo/Genotypefiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_chr",chr,".hg38.bim",sep="")) %>%
      dplyr::filter(V1 == chr & V4 > lead_pos-range & V4 < lead_pos+range) %>%
      as.data.frame() %>%
      .[,"V2"]
    snps_i = BEDMatrix(paste("/data01/DICE/inamo/Genotypefiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_chr",chr,".hg38.bed",sep="")) %>%
      as.matrix() %>%
      magrittr::set_colnames(paste(stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,1],
                                   stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,2],
                                   sep="_")) %>%
      .[, which (colnames(.) %in% variant)] %>%
      as.matrix() %>%
      t()%>%
      data.frame(id = variant, .) %>%
      magrittr::set_colnames(gsub("^X0.","",colnames(.)))
    sample_list = read.table("/data03/inamo/DICE/StringTie_PB29/sample_list.txt", header = TRUE)[,1] %>% as.character()
    IID_tr <- 
      read.table("/data01/DICE/DICE_DNA_RNAseq_correspondance.txt", header = TRUE) %>%
      dplyr::filter(Assay_Type == "RNA-Seq") %>%
      .[,c("Run","histological_type")] %>%
      magrittr::set_colnames(c("IID","stim")) %>%
      magrittr::set_rownames(.$IID) %>%
      .[sample_list,] %>%
      .[,c("stim","IID")]
    
    # Position files
    snpspos = data.frame(snp = snps_i$id,
                         chr = paste0("chr",str_split(snps_i$id, pattern = "_", simplify = TRUE)[,1]),
                         pos = as.integer(str_split(snps_i$id, pattern = "_", simplify = TRUE)[,2]))
    cov <- dplyr::bind_rows(as.data.frame(fread(paste("/data01/DICE/inamo/Genotypefiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup.hg38_PC10.txt",sep=""),header = TRUE)) %>%
                              dplyr::mutate(SampleID = str_split(SampleID, pattern = "_", simplify = TRUE)[,ncol(str_split(SampleID, pattern = "_", simplify = TRUE))]) %>%
                              magrittr::set_colnames(gsub("SampleID","id",colnames(.)%>%gsub("^0_","",.))) %>%
                              tibble::column_to_rownames("id") %>%
                              t() %>%
                              merge(.,read.table("/data01/DICE/DICE_DNA_RNAseq_correspondance.txt", header = TRUE),by.x="row.names",by.y="IID") %>%
                              tibble::column_to_rownames("Run") %>%
                              dplyr::select(starts_with("PC")) %>%
                              t() %>% as.data.frame() %>%
                              dplyr::mutate(id = rownames(.)) %>%
                              .[,c("id",sample_list)],
                            peer_factors %>%
                              dplyr::rename(id = PEERfactors) %>%
                              .[,grepl(paste0("^id$|^SRR"),colnames(.))])
    cov <-  data.frame(id = cov$id, impute::impute.knn(as.matrix(cov[,-1]), k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)$data ) %>%
      magrittr::set_colnames(colnames(cov))
    tr_ratio =
      fread(paste0("/data03/inamo/sQTL/leafcutter_gencode38/expression_quantile-rank-normalized/isoform_ratio/",tool,"/",sQTL_id,".txt")) %>% 
      dplyr::rename(gid = V1) %>%
      as.data.frame() %>%
      .[,grepl(paste0("^gid|^SRR"),colnames(.))]
    
    shared_samples = 
      intersect(
        corres[ grepl("_1_RNA_",corres$biospecimen_repository_sample_id) & corres$histological_type == cell , "Run"],
        colnames(tr_ratio) %>% .[grepl("^SRR",.)]
      )  %>%
      intersect(
        .,
        colnames(cov)
      )
    
    shared_IID = 
      corres %>%
      dplyr::filter(Run %in% shared_samples) %>%
      .$IID %>% as.character() %>%
      intersect(
        .,
        colnames(snps_i)
      )
    shared = 
      corres %>%
      dplyr::filter(Run %in% shared_samples) %>%
      dplyr::filter(IID %in% shared_IID) %>%
      magrittr::set_rownames(.$IID) %>%
      .[shared_IID,]
    
    snps = snps_i[,c("id",shared$IID %>% as.character())] %>%
      magrittr::set_colnames(c("id",shared$Run %>% as.character()))
    expr2 <- tr_ratio %>%
      as.data.frame() %>%
      .[, c("gid", shared$Run %>% as.character()) ]
    cov <- cov[,c("id",shared$Run %>% as.character())]
    
  }
  
  if ( nrow(expr2) > 0 &&
       all(colnames(snps)[-1] == colnames(expr2)[-1] &&
           colnames(snps)[-1] == colnames(cov)[-1]) ) {
    
    # MatrixQTLでmQTL-SNPsと発現量が相関する遺伝子を抽出
    write.table(snps, tmp_snps, sep = "\t",quote=FALSE, row.names = FALSE)
    write.table(expr2, tmp_exp, sep = "\t",quote=FALSE, row.names = FALSE)
    write.table(cov, tmp_cov, sep = "\t",quote=FALSE, row.names = FALSE)
    
    # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
    useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
    
    # Genotype file name
    SNP_file_name = tmp_snps;
    
    # Gene expression file name
    expression_file_name = tmp_exp;
    
    # Covariates file name
    # Set to character() for no covariates
    covariates_file_name = tmp_cov;
    
    # Error covariance matrix
    # Set to numeric() for identity.
    errorCovariance = numeric();
    
    pvOutputThreshold = 1;
    
    # Output file name
    output_file_name = tempfile();
    
    ## Load genotype data
    SNPs = SlicedData$new();
    SNPs$fileDelimiter = "\t";      # the TAB character
    SNPs$fileOmitCharacters = "NA"; # denote missing values;
    SNPs$fileSkipRows = 1;          # one row of column labels
    SNPs$fileSkipColumns = 1;       # one column of row labels
    SNPs$fileSliceSize = 2000;      # read file in slices of 2,000 rows
    SNPs$LoadFile(SNP_file_name);
    
    ## Load gene expression data
    gene = SlicedData$new();
    gene$fileDelimiter = "\t";      # the TAB character
    gene$fileOmitCharacters = "NA"; # denote missing values;
    gene$fileSkipRows = 1;          # one row of column labels
    gene$fileSkipColumns = 1;       # one column of row labels
    gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
    gene$LoadFile(expression_file_name);
    
    ## Load covariates
    cvrt = SlicedData$new();
    cvrt$fileDelimiter = "\t";      # the TAB character
    cvrt$fileOmitCharacters = "NA"; # denote missing values;
    cvrt$fileSkipRows = 1;          # one row of column labels
    cvrt$fileSkipColumns = 1;       # one column of row labels
    if(length(covariates_file_name)>0) {
      cvrt$LoadFile(covariates_file_name);
    }
    
    
    ## Run the analysis
    
    me = Matrix_eQTL_engine(
      snps = SNPs,
      gene = gene,
      cvrt = cvrt,
      output_file_name = output_file_name,
      pvOutputThreshold = pvOutputThreshold,
      useModel = useModel,
      errorCovariance = errorCovariance,
      verbose = FALSE,
      pvalue.hist = TRUE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE);
    
    unlink(output_file_name);
    
    ## Results:
    cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
    
  } else {
    stop("confirm genotype and isoform")
  }
  
  merged = dplyr::left_join(me$all$eqtls, rsid2chrpos, by="snps") %>%
    dplyr::select(rsid,pvalue,snps) %>%
    dplyr::rename(pval = pvalue,
                  chr_pos = snps) %>%
    dplyr::left_join(.,
                     fread(paste0("/data03/inamo/nanopore/RNA/210931_PB29_fastq/shiny/data/gwas/",gwas,"_GRCh38_short.txt")) %>%
                       dplyr::select(c(chr_pos,pval)),
                     by="chr_pos",suffix=c("_qtl","_gwas")) %>%
    dplyr::mutate(chr = stringr::str_split(chr_pos,pattern="_",simplify=TRUE) %>% .[,1],
                  pos = stringr::str_split(chr_pos,pattern="_",simplify=TRUE) %>% .[,2],
                  logp_qtl = -log10(pval_qtl),
                  logp_gwas = -log10(pval_gwas),
                  lead = logp_qtl*logp_gwas) %>%
    dplyr::arrange(-lead) %>%
    dplyr::mutate(rank = dplyr::row_number(-lead),
                  label = ifelse(rank==1,rsid,""),
                  color = ifelse(rank==1,"red","black")) %>%
    na.omit()
  
  g = locuscompare(in_fn1 = merged[,c("rsid","pval_qtl")] %>% dplyr::rename(pval = pval_qtl), 
                   marker_col1 = "rsid", pval_col1 = "pval", 
                   title = paste0(sQTL_id,"-sQTL [",cell,"]"),
                   in_fn2 = merged[,c("rsid","pval_gwas")] %>% dplyr::rename(pval = pval_gwas), 
                   marker_col2 = "rsid", pval_col2 = "pval", 
                   title2 = paste0(gwas,"-GWAS"), 
                   population = population,
                   combine = combine, 
                   legend = legend, legend_position = legend_position, 
                   lz_ylab_linebreak = lz_ylab_linebreak, genome = genome)
  assign(paste0("g_",cell),g)
}
png(file=paste0("~/exp.png"), width = 18000, height = 7000, res=720);
ggpubr::ggarrange(g_NONCLASSICAL_MONOCYTES,g_CLASSICAL_MONOCYTES,g_Mono_NS,
                  ncol=3,nrow=1)
dev.off()


# sQTL plot

suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(BEDMatrix))
suppressMessages(library(dplyr))
suppressMessages(library(magrittr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))

sQTL_id = "21:33252830:33268394:clu_34140_+" 
cell = "NONCLASSICAL_MONOCYTES"
# cell = "CLASSICAL_MONOCYTES"
# cell = "Mono_NS"
variant="21_33242905"
tool = "leafcutter"
chr = sQTL_id %>% stringr::str_split(., pattern = ":", simplify = TRUE) %>% .[,1] %>% paste0("chr",.)
start =  sQTL_id %>% stringr::str_split(., pattern = ":", simplify = TRUE) %>% .[,2] %>% as.integer()
end =  sQTL_id %>% stringr::str_split(., pattern = ":", simplify = TRUE) %>% .[,3] %>% as.integer()

print(sQTL_id)
print(chr)
print(start)
print(end)

print("creating expression data...")
for (cell in c("NONCLASSICAL_MONOCYTES","CLASSICAL_MONOCYTES","Mono_NS")){
  print(cell)
  # GEUV
  if (grepl("LCL", cell)){
    
    snps = BEDMatrix(paste("/data03/inamo/GEUV/VCF/QC_MAF0.05_HWE1E06_nodup_GRCh38_EUR_",chr,sep="")) %>%
      as.matrix() %>%
      magrittr::set_colnames(paste(stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,1],
                                   stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,2],
                                   sep="_")) %>%
      .[, which (
        (colnames(.) %>% stringr::str_split(., pattern = "_", simplify = TRUE) %>% .[,2] %>% as.integer() > start-1e06) &
          (colnames(.) %>% stringr::str_split(., pattern = "_", simplify = TRUE) %>% .[,2] %>% as.integer() < end+1e06)
      )] %>%
      as.matrix() %>%
      t() %>%
      data.frame(id = rownames(.), .) %>%
      magrittr::set_colnames(gsub("^X0.","",colnames(.))) %>%
      magrittr::set_colnames(gsub("\\.","-",colnames(.))) %>%
      tidyr::pivot_longer(cols = -c(id)) %>%
      magrittr::set_colnames(c("snp","sample","genotype"))
    expr2 <- 
      fread(paste0("/data03/inamo/sQTL/leafcutter_gencode38/expression_quantile-rank-normalized/isoform_ratio/",tool,"/",sQTL_id,".txt")) %>% 
      dplyr::rename(gid = V1) %>%
      as.data.frame() %>%
      .[,grepl(paste0("^gid|^ERR"),colnames(.))] %>%
      .[,-1] %>%
      t() %>%
      as.data.frame() %>%
      dplyr::mutate(sample = rownames(.)) %>%
      magrittr::set_colnames(c("expression","sample"))
    
    data_long <- merge(snps, expr2, by="sample") %>%
      merge(., 
            fread(paste("/data03/inamo/GEUV/VCF/QC_MAF0.05_HWE1E06_nodup_GRCh38_EUR_",chr,".bim",sep="")) %>%
              dplyr::filter(V1 == gsub("chr","",chr)) %>%
              dplyr::select(V2,V5,V6) %>%
              magrittr::set_colnames(c("snp","A1","A0")),
            by="snp") %>%
      dplyr::mutate(genotype = as.factor(genotype),
                    condition = cell)
    
    # EGA
  } else if (grepl("Mono_", cell)) {
    
    num = grep(cell,c("Mono_NS","Mono_LPS","Mono_Pam3CSK4","Mono_R848","Mono_IAV"))
    snps = BEDMatrix(paste("/data01/EGA/inamo/GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-",num,"_EUR_",chr,".hg38",sep="")) %>%
      as.matrix() %>%
      magrittr::set_colnames(paste(stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,1],
                                   stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,2],
                                   sep="_")) %>%
      .[, which (
        (colnames(.) %>% stringr::str_split(., pattern = "_", simplify = TRUE) %>% .[,2] %>% as.integer() > start-1e06) &
          (colnames(.) %>% stringr::str_split(., pattern = "_", simplify = TRUE) %>% .[,2] %>% as.integer() < end+1e06)
      )] %>%
      as.matrix() %>%
      t() %>%
      data.frame(id = rownames(.), .) %>%
      magrittr::set_colnames(gsub("^X0.","",colnames(.))) %>%
      magrittr::set_colnames(gsub("\\.","@",colnames(.))) %>%
      tidyr::pivot_longer(cols = -c(id)) %>%
      magrittr::set_colnames(c("snp","sample","genotype"))
    expr2 <- 
      fread(paste0("/data03/inamo/sQTL/leafcutter_gencode38/expression_quantile-rank-normalized/isoform_ratio/",tool,"/",sQTL_id,".txt")) %>% 
      dplyr::rename(gid = V1) %>%
      as.data.frame() %>%
      .[,grepl(paste0("^gid|",cell),colnames(.))] %>%
      .[,-1] %>%
      t() %>%
      as.data.frame() %>%
      dplyr::mutate(sample = stringr::str_split(rownames(.),pattern=";",simplify=TRUE) %>% .[,1]) %>%
      magrittr::set_colnames(c("expression","sample"))
    
    data_long <- merge(snps, expr2, by="sample") %>%
      merge(., 
            fread(paste("/data01/EGA/inamo/GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-",num,"_EUR_",chr,".hg38.bim",sep="")) %>%
              dplyr::filter(V1 == gsub("chr","",chr)) %>%
              dplyr::select(V2,V5,V6) %>%
              magrittr::set_colnames(c("snp","A1","A0")),
            by="snp") %>%
      dplyr::mutate(genotype = as.factor(genotype),
                    condition = cell)
    
    # ImmVar
  } else if ( grepl("MoDC_", cell) | grepl("CD4T_", cell) ) {
    
    snps = BEDMatrix(paste("/data01/ImmVar/inamo/ImmVar/Genotypefiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_",cell,"_EUR_",chr,".hg38",sep="")) %>%
      as.matrix() %>%
      magrittr::set_colnames(paste(stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,1],
                                   stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,2],
                                   sep="_")) %>%
      .[, which (
        (colnames(.) %>% stringr::str_split(., pattern = "_", simplify = TRUE) %>% .[,2] %>% as.integer() > start-1e06) &
          (colnames(.) %>% stringr::str_split(., pattern = "_", simplify = TRUE) %>% .[,2] %>% as.integer() < end+1e06)
      )] %>%
      as.matrix() %>%
      t() %>%
      data.frame(id = rownames(.), .) %>%
      magrittr::set_colnames(gsub("^X0.","",colnames(.))) %>%
      magrittr::set_colnames(gsub("\\.","-",colnames(.))) %>%
      tidyr::pivot_longer(cols = -c(id)) %>%
      magrittr::set_colnames(c("snp","sample","genotype"))
    expr2 <- 
      fread(paste0("/data03/inamo/sQTL/leafcutter_gencode38/expression_quantile-rank-normalized/isoform_ratio/",tool,"/",sQTL_id,".txt")) %>% 
      dplyr::rename(gid = V1) %>%
      as.data.frame() %>%
      .[,grepl(paste0("^gid|",cell),colnames(.))] %>%
      .[,-1] %>%
      t() %>%
      as.data.frame() %>%
      dplyr::mutate(sample = stringr::str_split(rownames(.),pattern=";",simplify=TRUE) %>% .[,1]) %>%
      magrittr::set_colnames(c("expression","sample"))
    
    data_long <- merge(snps, expr2, by="sample") %>%
      merge(., 
            fread(paste("/data01/ImmVar/inamo/ImmVar/Genotypefiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_",cell,"_EUR_",chr,".hg38.bim",sep="")) %>%
              dplyr::filter(V1 == gsub("chr","",chr)) %>%
              dplyr::select(V2,V5,V6) %>%
              magrittr::set_colnames(c("snp","A1","A0")),
            by="snp") %>%
      dplyr::mutate(genotype = as.factor(genotype),
                    condition = cell)
    
    # DICE
  } else {
    
    snps = BEDMatrix(paste("/data01/DICE/inamo/Genotypefiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_",chr,".hg38.bed",sep="")) %>%
      as.matrix() %>%
      magrittr::set_colnames(paste(stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,1],
                                   stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,2],
                                   sep="_")) %>%
      .[, which (
        (colnames(.) %>% stringr::str_split(., pattern = "_", simplify = TRUE) %>% .[,2] %>% as.integer() > start-1e06) &
          (colnames(.) %>% stringr::str_split(., pattern = "_", simplify = TRUE) %>% .[,2] %>% as.integer() < end+1e06)
      )] %>%
      as.matrix() %>%
      t() %>%
      data.frame(id = rownames(.), .) %>%
      magrittr::set_colnames(gsub("^X0.","",colnames(.))) %>%
      magrittr::set_colnames(gsub("\\.","-",colnames(.))) %>%
      tidyr::pivot_longer(cols = -c(id)) %>%
      magrittr::set_colnames(c("snp","sample","genotype"))
    
    shared_samples = 
      intersect(
        corres[ grepl("_1_RNA_",corres$biospecimen_repository_sample_id) & corres$histological_type == cell , "Run"],
        colnames(fread("/data03/inamo/DICE/StringTie_PB29/transcript_ratio.txt")[,-2] %>%
                   dplyr::rename(gid = query_name)) %>% .[grepl("^SRR",.)] %>% stringr::str_split(., pattern = ";", simplify = TRUE) %>% .[,1]
      ) 
    shared_IID = 
      corres %>%
      dplyr::filter(Run %in% shared_samples) %>%
      .$IID %>% as.character() %>%
      intersect(
        .,
        snps$sample
      )
    shared = 
      corres %>%
      dplyr::filter(Run %in% shared_samples) %>%
      dplyr::filter(IID %in% shared_IID) %>%
      magrittr::set_rownames(.$IID) %>%
      .[shared_IID,]
    
    snps = merge(snps,shared[,c("Run","IID")],by.x="sample",by.y="IID") %>%
      dplyr::mutate(sample=Run) %>%
      dplyr::select(-Run)
    expr2 <- 
      fread(paste0("/data03/inamo/sQTL/leafcutter_gencode38/expression_quantile-rank-normalized/isoform_ratio/",tool,"/",sQTL_id,".txt")) %>% 
      dplyr::rename(gid = V1) %>%
      as.data.frame() %>%
      .[,grepl(paste0("^gid|^SRR"),colnames(.))] %>%
      .[,-1] %>%
      t() %>%
      as.data.frame() %>%
      dplyr::mutate(sample = rownames(.)) %>%
      magrittr::set_colnames(c("expression","sample"))
    
    data_long <- merge(snps, expr2, by="sample") %>%
      merge(., 
            fread(paste("/data01/DICE/inamo/Genotypefiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_",chr,".hg38.bim",sep="")) %>%
              dplyr::filter(V1 == gsub("chr","",chr)) %>%
              dplyr::select(V2,V5,V6) %>%
              magrittr::set_colnames(c("snp","A1","A0")),
            by="snp") %>%
      dplyr::mutate(genotype = as.factor(genotype),
                    condition = cell)
    
    
  }
  
  print(head(data_long))
  # if (!file.exists(paste0("/data03/inamo/sQTL/leafcutter_gencode38/expression/",sQTL_id,"_",cell,"_",tool,".txt"))){ 
  #   write.table(data_long, paste0("/data03/inamo/sQTL/leafcutter_gencode38/expression/",sQTL_id,"_",cell,"_",tool,".txt"), quote = FALSE, row.names = FALSE)
  # }
  
  print("plotting...")
  if (variant %in% data_long$snp){
    data_long_tmp = data_long %>%
      dplyr::filter(condition == cell & snp == variant)
  } else {
    stop(paste0("selected variant is not within cis region of ",sQTL_id))
  }
  
  A1 = data_long_tmp %>%
    as.data.frame() %>%
    .[1,"A1"]
  A0 = data_long_tmp %>%
    as.data.frame() %>%
    .[1,"A0"]
  
  if (tool == "StringTie" | tool == "kallisto" | tool == "kallisto-cds-clustered") {
    y_label = paste0(sQTL_id,"\n[normalized isoform ratio]") 
  } else if (tool == "leafcutter") {
    y_label = paste0(sQTL_id,"\n[normalized PSI]") 
  }
  
  fig = ggplot(data_long_tmp, 
               aes(genotype, expression)) +
    geom_jitter(colour = "darkorange",alpha = 0.3, position=position_jitter(width=0.25)) +
    geom_boxplot(alpha = 0.5, outlier.size=0, fill = "steelblue") +
    geom_smooth(method = 'lm',formula = formula, col="darkred", aes(group=1), se=FALSE) +
    annotation_custom(grid::grobTree(grid::textGrob(paste("p=",signif(summary(lm(as.numeric(data_long_tmp$expression)~ as.numeric(data_long_tmp$genotype)))$coefficients[2,4],digits=2)),
                                                    x = unit(0.46, "npc"), y = unit(0.95, "npc"), hjust=0,
                                                    gp=grid::gpar(col="black", fontsize=10, fontface="italic")))) +
    facet_wrap(~condition) + 
    theme_classic() +
    scale_x_discrete(labels= c(paste(A0,A0,sep="/"),paste(A1,A0,sep="/"),paste(A1,A1,sep="/"))) +
    ylab(y_label) +
    xlab(paste0(variant)) +
    theme(strip.text.x=element_text(size=9, color="black", face="bold"),
          strip.text.y=element_text(size=9, color="black", face="bold"),
          legend.position = "bottom",
          plot.title = element_text(size=8),
          axis.title.x = element_text(size=8),
          axis.title.y = element_text(size =8),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          legend.text =  element_text(size = 10), 
          legend.key.size = grid::unit(0.8, "lines"),
          legend.title = element_text(size = 8, hjust = 0))
  assign(paste0("fig_",cell),fig)
}
png(file=paste0("~/exp.png"), width = 9000, height = 7000, res=720);
ggpubr::ggarrange(fig_NONCLASSICAL_MONOCYTES,fig_CLASSICAL_MONOCYTES,fig_Mono_NS,
                  ncol=3,nrow=1)
dev.off()
