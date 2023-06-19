# PEER normalization
library(limma)
library(peer)
library(GenABEL)
library(genefilter)
library(impute)
library(data.table)
library(magrittr)

# covariate file
## EGA
SUBSET=1:5
POP = "EUR"
for (i in 1:length(SUBSET)){
  for (p in 1:length(POP)){
    tryCatch({
      
      stimulus = c("Mono_NS","Mono_LPS","Mono_Pam3CSK4","Mono_R848","Mono_IAV")[i]
      cov = as.data.frame(fread(paste("/data01/EGA/inamo/GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-",SUBSET[i],"_",POP[p],"_PC10.txt",sep=""),header = TRUE))
      cov$SampleID <- stringr::str_split(cov$SampleID, pattern = "_", simplify = TRUE)[,ncol(stringr::str_split(cov$SampleID, pattern = "_", simplify = TRUE))]
      colnames(cov) <- gsub("SampleID","id",colnames(cov)%>%gsub("^0_","",.))
      assign(paste0("cov_",stimulus), cov)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

corres = read.table("/data01/DICE/DICE_DNA_RNAseq_correspondance.txt", header = TRUE)
sample_list = stringr::str_split(list.files("/data03/inamo/DICE/StringTie_PB29/ballgown", pattern = "t_data.ctab", full.names=T, recursive=T), pattern = "/", simplify = TRUE)[,7]
cov = as.data.frame(fread(paste("/data01/DICE/inamo/Genotypefiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup.hg38_PC10.txt",sep=""),header = TRUE)) %>%
  dplyr::mutate(SampleID = stringr::str_split(SampleID, pattern = "_", simplify = TRUE)[,ncol(stringr::str_split(SampleID, pattern = "_", simplify = TRUE))]) %>%
  magrittr::set_colnames(gsub("SampleID","id",colnames(.)%>%gsub("^0_","",.))) %>%
  tibble::column_to_rownames("id") %>%
  t() %>%
  merge(.,read.table("/data01/DICE/DICE_DNA_RNAseq_correspondance.txt", header = TRUE),by.x="row.names",by.y="IID") %>%
  tibble::column_to_rownames("Run") %>%
  dplyr::select(starts_with("PC")) %>%
  t() %>% as.data.frame() %>%
  dplyr::mutate(id = rownames(.)) %>%
  .[,c("id",sample_list)]
assign(paste0("cov_DICE"), cov)

# common intron excision ratio
print (paste0 ("Loading intron excision ratio...\n") );
lf <- list.files(path = "/data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter", pattern = paste0("all_perind.counts.gz.qqnorm_chrchr*"), full.names = T)
data <- as.data.frame(do.call(rbind, lapply(lf,fread)))
data.frame(
  geneid = data[,4],
  chr = paste0("chr",data[,1]),
  s1 = data[,2],
  s2 = data[,3]
) %>%
  write.table(., 
              paste0("/data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/isoform_position.txt"),
              sep = "\t",quote=FALSE, row.names = FALSE)
data[,4:ncol(data)] %>%
  magrittr::set_colnames(gsub("-1_sorted_filtered.bam",";Mono_NS",colnames(.))) %>%
  magrittr::set_colnames(gsub("-2_sorted_filtered.bam",";Mono_LPS",colnames(.))) %>%
  magrittr::set_colnames(gsub("-3_sorted_filtered.bam",";Mono_Pam3CSK4",colnames(.))) %>%
  magrittr::set_colnames(gsub("-4_sorted_filtered.bam",";Mono_R848",colnames(.))) %>%
  magrittr::set_colnames(gsub("-5_sorted_filtered.bam",";Mono_IAV",colnames(.))) %>%
  magrittr::set_colnames(gsub("_sorted_filtered.bam","",colnames(.))) %>%
  # magrittr::set_colnames(gsub("\\-CD4T-",";CD4T_",colnames(.))) %>%
  # magrittr::set_colnames(gsub("\\-MoDC-",";MoDC_",colnames(.))) %>%
  magrittr::set_colnames(gsub("^ID","gid",colnames(.))) %>%
  dplyr::select(-contains(".")) %>%
  write.table(., 
              paste0("/data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/psi_sr-RNAseq.txt"),
              sep = "\t",quote=FALSE, row.names = FALSE)

lf <- list.files(path = "/data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter", pattern = paste0("all_perind.counts.gz.phen_chrchr*"), full.names = T)
data <- as.data.frame(do.call(rbind, lapply(lf,fread)))
data[,4:ncol(data)] %>%
  magrittr::set_colnames(gsub("-1_sorted_filtered.bam",";Mono_NS",colnames(.))) %>%
  magrittr::set_colnames(gsub("-2_sorted_filtered.bam",";Mono_LPS",colnames(.))) %>%
  magrittr::set_colnames(gsub("-3_sorted_filtered.bam",";Mono_Pam3CSK4",colnames(.))) %>%
  magrittr::set_colnames(gsub("-4_sorted_filtered.bam",";Mono_R848",colnames(.))) %>%
  magrittr::set_colnames(gsub("-5_sorted_filtered.bam",";Mono_IAV",colnames(.))) %>%
  magrittr::set_colnames(gsub("_sorted_filtered.bam","",colnames(.))) %>%
  # magrittr::set_colnames(gsub("\\-CD4T-",";CD4T_",colnames(.))) %>%
  # magrittr::set_colnames(gsub("\\-MoDC-",";MoDC_",colnames(.))) %>%
  magrittr::set_colnames(gsub("^ID","gid",colnames(.))) %>%
  dplyr::select(-contains(".")) %>%
  write.table(., 
              paste0("/data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/psi01_sr-RNAseq.txt"),
              sep = "\t",quote=FALSE, row.names = FALSE)

# leafcutter
cohort="EGA"
for (cohort in c(
  "EGA","DICE"
  # ,"GEUV","ImmVar"
  )){
  data <- fread("/data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/psi01_sr-RNAseq.txt") %>%
    tibble::column_to_rownames("gid") %>%
    # magrittr::set_colnames(gsub("\\-CD4T-",";CD4T_",colnames(.))) %>%
    # magrittr::set_colnames(gsub("\\-MoDC-",";MoDC_",colnames(.))) %>%
    magrittr::set_colnames(gsub("-1$",";Mono_NS",colnames(.))) %>%
    magrittr::set_colnames(gsub("-2$",";Mono_LPS",colnames(.))) %>%
    magrittr::set_colnames(gsub("-3$",";Mono_Pam3CSK4",colnames(.))) %>%
    magrittr::set_colnames(gsub("-4$",";Mono_R848",colnames(.))) %>%
    magrittr::set_colnames(gsub("-5$",";Mono_IAV",colnames(.))) %>%
    as.data.frame()
  dim(data)
  
  if (cohort == "EGA"){
    for (stim in 
         c("Mono_NS","Mono_LPS","Mono_Pam3CSK4","Mono_R848","Mono_IAV")){
      tryCatch({
        if (grepl("Mono_", stim)){
          data_imp_sel = data[,grepl(stim,colnames(data))]
          cov <- eval(parse(text=paste0("cov_",stim)))
          shared_samples = 
            intersect(
              colnames(data_imp_sel) %>% stringr::str_split(., pattern = ";", simplify = TRUE) %>% .[,1],
              colnames(cov)
            ) 
          data_imp_sel = data_imp_sel[,paste0(shared_samples,";",stim)]
        } else if (grepl("MoDC_", stim) | grepl("CD4T_", stim)){
          data_imp_sel = data[,grepl(stim,colnames(data))]
          cov <- eval(parse(text=paste0("cov_",stim)))
          shared_samples = 
            intersect(
              colnames(data_imp_sel) %>% stringr::str_split(., pattern = ";", simplify = TRUE) %>% .[,1],
              colnames(cov)
            ) 
          data_imp_sel = data_imp_sel[,paste0(shared_samples,";",stim)]
        } else if (grepl("LCL", stim)){
          data_imp_sel = data[,grepl("^ERR",colnames(data))]
          cov <- eval(parse(text=paste0("cov_GEUV")))
          shared_samples = 
            intersect(
              colnames(data_imp_sel) %>% stringr::str_split(., pattern = ";", simplify = TRUE) %>% .[,1],
              colnames(cov)
            ) 
          data_imp_sel = data_imp_sel[,shared_samples]
        } else {
          data_imp_sel = data[,grepl("^SRR",colnames(data))]
          cov <- eval(parse(text=paste0("cov_DICE")))
          shared_samples = 
            intersect(
              corres[ grepl("_1_RNA_",corres$biospecimen_repository_sample_id) & corres$histological_type == stim , "Run"],
              colnames(data_imp_sel)
            )  %>%
            intersect(
              .,
              colnames(cov)
            )
          data_imp_sel = data_imp_sel[,shared_samples]
        }
        
        cov = cov[,shared_samples]
        # quantile normalization
        exp_q <- normalizeQuantiles(as.matrix(data_imp_sel))
        bin  <- apply( exp_q, 1, function(x){as.character(x);length(unique(x))} )
        var  <- apply( exp_q, 1, var )
        p1  = names( bin[bin>2] )
        p2  = names( var[var>0] )
        p_m = intersect(p1, p2)
        exp_q_f <- exp_q[ p_m, ]
        # rank normalization
        rank_n  <- apply( exp_q_f, 1, function(x){rntransform(x)})  # applyするとr,cが逆転
        rank_n <- t(rank_n)
        rownames(rank_n) <- rownames(exp_q_f)
        colnames(rank_n) <- colnames(exp_q_f)
        # peer normalization
        ## build model
        model = PEER()
        ## Set the maximum number of unobserved factors to model.
        PEER_setNk(model,15)
        ## Set expression data.
        PEER_setPhenoMean(model,t(rank_n))
        ## (Optional) Set observed covariates.
        PEER_setCovariates(model, t(cov))
        ## Train the model, observing convergence:
        PEER_setNmax_iterations(model, 10000)
        PEER_update(model)
        factors = PEER_getX(model)
        dim(factors)
        factors = t(factors) %>%
          magrittr::set_colnames(colnames(data_imp_sel)) %>%
          as.data.frame() %>%
          .[1:15,] %>%
          dplyr::mutate(id = paste0(rep("PEER_factor",15),1:15)) %>%
          .[,c("id",colnames(.)[!grepl("^id$",colnames(.))])]
        weights = PEER_getW(model)
        dim(weights)
        rownames(weights) <- rownames(rank_n)
        colnames(weights) <- paste0(rep("PEER_factor",25),1:25)
        precision = PEER_getAlpha(model)
        dim(precision)
        residuals = PEER_getResiduals(model)
        dim(residuals)
        residuals = t(residuals)
        rownames(residuals) <- rownames(rank_n)
        colnames(residuals) <- colnames(data_imp_sel)
        write.table(rank_n,paste0("/data03/inamo/",cohort,"/leafcutter_gencode38/quantile-rank-normalized_",stim,"_psi.txt"),sep="\t",row.names = TRUE,col.names = NA,quote = FALSE)
        write.table(factors,paste0("/data03/inamo/",cohort,"/leafcutter_gencode38/PEERfactors_",stim,"_psi.txt"),row.names = FALSE,sep="\t",quote = FALSE)
        write.table(weights,paste0("/data03/inamo/",cohort,"/leafcutter_gencode38/PEERweights_",stim,"_psi.txt"),sep="\t",row.names = TRUE,col.names = NA,quote = FALSE)
        write.table(precision,paste0("/data03/inamo/",cohort,"/leafcutter_gencode38/PEERprecision_",stim,"_psi.txt"),sep="\t",row.names = TRUE,col.names = NA,quote = FALSE)
        write.table(residuals,paste0("/data03/inamo/",cohort,"/leafcutter_gencode38/PEERresiduals_",stim,"_psi.txt"),sep="\t",row.names = TRUE,col.names = NA,quote = FALSE)
        print(paste0(stim,", done..."))
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
  } else {
    for (stim in 
         c(
           # "TFH","TH1","TH17","TH1-17","TH2","TREG_MEMORY","TREG_NAIVE","B_NAIVE","CD4_NAIVE","CD4_N_STIM","CD8_NAIVE","CD8_N_STIM","NK_CD16POS",
           "NONCLASSICAL_MONOCYTES","CLASSICAL_MONOCYTES"
           )){
      tryCatch({
        if (grepl("Mono_", stim)){
          data_imp_sel = data[,grepl(stim,colnames(data))]
          cov <- eval(parse(text=paste0("cov_",stim)))
          shared_samples = 
            intersect(
              colnames(data_imp_sel) %>% stringr::str_split(., pattern = ";", simplify = TRUE) %>% .[,1],
              colnames(cov)
            ) 
          data_imp_sel = data_imp_sel[,paste0(shared_samples,";",stim)]
        } else if (grepl("MoDC_", stim) | grepl("CD4T_", stim)){
          data_imp_sel = data[,grepl(stim,colnames(data))]
          cov <- eval(parse(text=paste0("cov_",stim)))
          shared_samples = 
            intersect(
              colnames(data_imp_sel) %>% stringr::str_split(., pattern = ";", simplify = TRUE) %>% .[,1],
              colnames(cov)
            ) 
          data_imp_sel = data_imp_sel[,paste0(shared_samples,";",stim)]
        } else if (grepl("LCL", stim)){
          data_imp_sel = data[,grepl("^ERR",colnames(data))]
          cov <- eval(parse(text=paste0("cov_GEUV")))
          shared_samples = 
            intersect(
              colnames(data_imp_sel) %>% stringr::str_split(., pattern = ";", simplify = TRUE) %>% .[,1],
              colnames(cov)
            ) 
          data_imp_sel = data_imp_sel[,shared_samples]
        } else {
          data_imp_sel = data[,grepl("^SRR",colnames(data))]
          cov <- eval(parse(text=paste0("cov_DICE")))
          shared_samples = 
            intersect(
              corres[ grepl("_1_RNA_",corres$biospecimen_repository_sample_id) & corres$histological_type == stim , "Run"],
              colnames(data_imp_sel)
            )  %>%
            intersect(
              .,
              colnames(cov)
            )
          data_imp_sel = data_imp_sel[,shared_samples]
        }
        cov = cov[,shared_samples]
        # quantile normalization
        exp_q <- normalizeQuantiles(as.matrix(data_imp_sel))
        bin  <- apply( exp_q, 1, function(x){as.character(x);length(unique(x))} )
        var  <- apply( exp_q, 1, var )
        p1  = names( bin[bin>2] )
        p2  = names( var[var>0] )
        p_m = intersect(p1, p2)
        exp_q_f <- exp_q[ p_m, ]
        # rank normalization
        rank_n  <- apply( exp_q_f, 1, function(x){rntransform(x)})  # applyするとr,cが逆転
        rank_n <- t(rank_n)
        rownames(rank_n) <- rownames(exp_q_f)
        colnames(rank_n) <- colnames(exp_q_f)
        # peer normalization
        ## build model
        model = PEER()
        ## Set the maximum number of unobserved factors to model.
        PEER_setNk(model,15)
        ## Set expression data.
        PEER_setPhenoMean(model,t(rank_n))
        ## (Optional) Set observed covariates.
        PEER_setCovariates(model, t(cov))
        ## Train the model, observing convergence:
        PEER_setNmax_iterations(model, 10000)
        PEER_update(model)
        factors = PEER_getX(model)
        dim(factors)
        factors = t(factors) %>%
          magrittr::set_colnames(colnames(data_imp_sel)) %>%
          as.data.frame() %>%
          .[1:15,] %>%
          dplyr::mutate(id = paste0(rep("PEER_factor",15),1:15)) %>%
          .[,c("id",colnames(.)[!grepl("^id$",colnames(.))])]
        weights = PEER_getW(model)
        dim(weights)
        rownames(weights) <- rownames(rank_n)
        colnames(weights) <- paste0(rep("PEER_factor",25),1:25)
        precision = PEER_getAlpha(model)
        dim(precision)
        residuals = PEER_getResiduals(model)
        dim(residuals)
        residuals = t(residuals)
        rownames(residuals) <- rownames(rank_n)
        colnames(residuals) <- colnames(data_imp_sel)
        colnames(rank_n) <- colnames(data_imp_sel)
        write.table(rank_n,paste0("/data03/inamo/",cohort,"/leafcutter_gencode38/quantile-rank-normalized_",stim,"_psi.txt"),sep="\t",row.names = TRUE,col.names = NA,quote = FALSE)
        write.table(factors,paste0("/data03/inamo/",cohort,"/leafcutter_gencode38/PEERfactors_",stim,"_psi.txt"),sep="\t",row.names = FALSE,quote = FALSE)
        write.table(weights,paste0("/data03/inamo/",cohort,"/leafcutter_gencode38/PEERweights_",stim,"_psi.txt"),sep="\t",row.names = TRUE,col.names = NA,quote = FALSE)
        write.table(precision,paste0("/data03/inamo/",cohort,"/leafcutter_gencode38/PEERprecision_",stim,"_psi.txt"),sep="\t",row.names = TRUE,col.names = NA,quote = FALSE)
        write.table(residuals,paste0("/data03/inamo/",cohort,"/leafcutter_gencode38/PEERresiduals_",stim,"_psi.txt"),sep="\t",row.names = TRUE,col.names = NA,quote = FALSE)
        print(paste0(stim,", done..."))
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
  }
}

for (stim in 
     c(
       # "CD4T_Activated","MoDC_unstim","MoDC_FLU","MoDC_IFNb",
       "Mono_NS",
       # "Mono_LPS","Mono_Pam3CSK4","Mono_R848","Mono_IAV",
       # "LCL"
       # "TFH","TH1","TH17","TH1-17","TH2","TREG_MEMORY","TREG_NAIVE","B_NAIVE","CD4_NAIVE","CD4_N_STIM","CD8_NAIVE","CD8_N_STIM","NK_CD16POS",
       "NONCLASSICAL_MONOCYTES","CLASSICAL_MONOCYTES"
       )){
  
  if (grepl("LCL", stim)){
    cohort="GEUV"
  } else if (grepl("Mono_", stim)){
    cohort="EGA"
  } else if ( grepl("MoDC_", stim) | grepl("CD4T_", stim) ) {
    cohort="ImmVar"
  } else {
    cohort="DICE"
  }
  
  if (stim=="Mono_NS"){
    print(stim)
    tmp = fread(paste0("/data03/inamo/",cohort,"/leafcutter_gencode38/quantile-rank-normalized_",stim,"_psi.txt"))
    print(dim(tmp))
  } else {
    print(stim)
    tmp = dplyr::full_join(tmp,
                           fread(paste0("/data03/inamo/",cohort,"/leafcutter_gencode38/quantile-rank-normalized_",stim,"_psi.txt")),
                           by="V1")
    print(dim(tmp))
  }
}
write.table(tmp, paste0("/data03/inamo/sQTL/leafcutter_gencode38/quantile-rank-normalized_psi.txt"),sep="\t",row.names = FALSE, quote = FALSE)

for (stim in 
     c(
       # "CD4T_Activated","MoDC_unstim","MoDC_FLU","MoDC_IFNb",
       "Mono_NS",
       # "Mono_LPS","Mono_Pam3CSK4","Mono_R848","Mono_IAV",
       # "LCL"
       # "TFH","TH1","TH17","TH1-17","TH2","TREG_MEMORY","TREG_NAIVE","B_NAIVE","CD4_NAIVE","CD4_N_STIM","CD8_NAIVE","CD8_N_STIM","NK_CD16POS",
       "NONCLASSICAL_MONOCYTES","CLASSICAL_MONOCYTES"
     )){
  
  if (grepl("LCL", stim)){
    cohort="GEUV"
  } else if (grepl("Mono_", stim)){
    cohort="EGA"
  } else if ( grepl("MoDC_", stim) | grepl("CD4T_", stim) ) {
    cohort="ImmVar"
  } else {
    cohort="DICE"
  }
  if (stim=="Mono_NS"){
    print(stim)
    tmp = fread(paste0("/data03/inamo/",cohort,"/leafcutter_gencode38/PEERfactors_",stim,"_psi.txt")) 
    print(dim(tmp))
  } else {
    print(stim)
    tmp = dplyr::full_join(tmp,
                           fread(paste0("/data03/inamo/",cohort,"/leafcutter_gencode38/PEERfactors_",stim,"_psi.txt")),
                           by="id")
    fread(paste0("/data03/inamo/",cohort,"/leafcutter_gencode38/PEERfactors_",stim,"_psi.txt")) 
    print(dim(tmp))
  }
}
write.table(tmp %>%
              dplyr::rename(PEERfactors = id),
            paste0("/data03/inamo/sQTL/leafcutter_gencode38/PEERfactors.txt"),sep="\t",row.names = FALSE,quote = FALSE)




