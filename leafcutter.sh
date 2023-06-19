#!/bin/sh
#$ -S /bin/sh

# sQTL

## gencode 38

## extract exon and make annotation code
grep -w "exon" /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/SQANTI3_output/gencode38_classification.filtered_lite.gtf | awk '{print $1 " " $4 " " $5 " " $7 " " $12}' | sed -e 's/"//g' | sed -e 's/;//g' > /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/SQANTI3_output/gencode38_classification.filtered_lite.exon.txt
~/tools/leafcutter/leafviz/gtf2leafcutter.pl -o /data02/home/juninamo/tools/leafcutter/annotation_codes/gencode38 /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/SQANTI3_output/gencode38_classification.filtered_lite.gtf


## # ImmVar
## for i in CD4T_Activated MoDC_unstim MoDC_FLU MoDC_IFNb; do
##     mkdir --parents /data03/inamo/ImmVar/leafcutter_gencode38/${i}
## 
##     cat /data01/ImmVar/inamo/ImmVar/${i}_RNAseq_IID.list | while read line; do
##         
##         tmp=$(echo -e "${i}")
##         cell=$(echo -e ${tmp} | sed -e "s/_/-/g");
##         LINE=$(cat /data01/ImmVar/inamo/ImmVar/SRA_IID_table.txt | grep -w ${line} | grep ${cell} | awk '{print $2}');
## 
##         if [ ! -e "/data03/inamo/ImmVar/BAM/${LINE}_STAR_GRCh38/${LINE}_sorted_GRCh38.bam.junc" ] && [ -e "/data03/inamo/ImmVar/BAM/${LINE}_STAR_GRCh38/${LINE}_sorted_GRCh38.bam.bai" ]; then
##         echo -e "regtools junctions extract: ${LINE}_sorted_GRCh38.bam.junc ..."
## 
##         regtools junctions extract \
##         -s 0 -a 8 -m 50 -M 500000 \
##         /data03/inamo/ImmVar/BAM/${LINE}_STAR_GRCh38/${LINE}_sorted_GRCh38.bam \
##         -o /data03/inamo/ImmVar/BAM/${LINE}_STAR_GRCh38/${LINE}_sorted_GRCh38.bam.junc
## 
##         cat /data03/inamo/ImmVar/BAM/${LINE}_STAR_GRCh38/${LINE}_sorted_GRCh38.bam.junc | grep "^chr" | grep -v "^chrM" > /data03/inamo/ImmVar/BAM/${LINE}_STAR_GRCh38/${LINE}_sorted_filtered.bam.junc
##         
##         fi
## 
##     done
## 
## done


# EGA
for i in {1..5}; do
    mkdir --parents /data03/inamo/EGA/leafcutter_gencode38/Mono_${i}

    cat /data01/EGA/inamo/fq_sample.list | grep -e "-${i}$" | while read line; do
        
        if [ ! -e "/data03/inamo/EGA/BAM/${line}_STAR_GRCh38/${line}_sorted_GRCh38.bam.junc" ] && [ -e "/data03/inamo/EGA/BAM/${line}_STAR_GRCh38/${line}_sorted_GRCh38.bam.bai" ]; then
        echo -e "regtools junctions extract: ${line}_sorted_GRCh38.bam.junc ..."

        regtools junctions extract \
        -s 0 -a 8 -m 50 -M 500000 \
        /data03/inamo/EGA/BAM/${line}_STAR_GRCh38/${line}_sorted_GRCh38.bam \
        -o /data03/inamo/EGA/BAM/${line}_STAR_GRCh38/${line}_sorted_GRCh38.bam.junc

        cat /data03/inamo/EGA/BAM/${line}_STAR_GRCh38/${line}_sorted_GRCh38.bam.junc | grep "^chr" | grep -v "^chrM" > /data03/inamo/EGA/BAM/${line}_STAR_GRCh38/${line}_sorted_filtered.bam.junc
        
        fi

    done

done


# DICE
cat /data01/DICE/inamo/cell_type.list | while read cell; do
    mkdir --parents /data03/inamo/DICE/leafcutter_gencode38/${cell}

    cat /data01/DICE/inamo/${cell}_sra.list | while read sra; do
        
        if [ ! -e "/data03/inamo/DICE/BAM/${cell}/${sra}_STAR_GRCh38/${sra}_sorted_GRCh38.bam.junc" ] && [ -e "/data03/inamo/DICE/BAM/${cell}/${sra}_STAR_GRCh38/${sra}_sorted_GRCh38.bam.bai" ]; then
        echo -e "regtools junctions extract: ${cell}, ${sra}_sorted_GRCh38.bam.junc ..."

        regtools junctions extract \
        -s 0 -a 8 -m 50 -M 500000 \
        /data03/inamo/DICE/BAM/${cell}/${sra}_STAR_GRCh38/${sra}_sorted_GRCh38.bam \
        -o /data03/inamo/DICE/BAM/${cell}/${sra}_STAR_GRCh38/${sra}_sorted_GRCh38.bam.junc

        cat /data03/inamo/DICE/BAM/${cell}/${sra}_STAR_GRCh38/${sra}_sorted_GRCh38.bam.junc | grep "^chr" | grep -v "^chrM" > /data03/inamo/DICE/BAM/${cell}/${sra}_STAR_GRCh38/${sra}_sorted_filtered.bam.junc
        
        fi

    done

done


## # GEUV
## cat /data03/inamo/GEUV/all_sample.txt | cut -f 1 | while read sample; do
##         
##     if [ ! -e "/data03/inamo/GEUV/BAM/${sample}_STAR_GRCh38/${sample}_sorted_GRCh38.bam.junc" ] && [ -e "/data03/inamo/GEUV/BAM/${sample}_STAR_GRCh38/${sample}_sorted_GRCh38.bam.bai" ]; then
##     echo -e "regtools junctions extract: ${sample}_sorted_GRCh38.bam.junc ..."
## 
##     regtools junctions extract \
##     -s 0 -a 8 -m 50 -M 500000 \
##     /data03/inamo/GEUV/BAM/${sample}_STAR_GRCh38/${sample}_sorted_GRCh38.bam \
##     -o /data03/inamo/GEUV/BAM/${sample}_STAR_GRCh38/${sample}_sorted_GRCh38.bam.junc
## 
##     cat /data03/inamo/GEUV/BAM/${sample}_STAR_GRCh38/${sample}_sorted_GRCh38.bam.junc | grep "^chr" | grep -v "^chrM" > /data03/inamo/GEUV/BAM/${sample}_STAR_GRCh38/${sample}_sorted_filtered.bam.junc
##         
##     fi
## 
## done



# EGA and ImmVar and DICE and GEUV

touch /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/all_juncfiles.txt
ls /data03/inamo/EGA/BAM/*_STAR_GRCh38/*_sorted_filtered.bam.junc >> /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/all_juncfiles.txt
## ls /data03/inamo/ImmVar/BAM/*_STAR_GRCh38/*_sorted_filtered.bam.junc >> /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/all_juncfiles.txt
ls /data03/inamo/DICE/BAM/*/*_STAR_GRCh38/*_sorted_filtered.bam.junc >> /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/all_juncfiles.txt
## ls /data03/inamo/GEUV/BAM/*_STAR_GRCh38/*_sorted_filtered.bam.junc >> /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/all_juncfiles.txt

python tools/leafcutter/clustering/leafcutter_cluster_regtools.py \
-j /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/all_juncfiles.txt \
-r /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/ \
-o all

export PYTHONPATH="$HOME/.pyenv/versions/2.7.18/lib/python2.7/site-packages"
.pyenv/versions/2.7.18/bin/python tools/leafcutter/scripts/prepare_phenotype_table_modified.py \
/data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/all_perind.counts.gz -p 10

# make group files
## scripts/nanopore/make_groups_leafcutter.R

## # ImmVar
## for cell in CD4T_Activated MoDC_unstim MoDC_FLU MoDC_IFNb; do
##     if [ ! -e "/data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_ImmVar_${cell}_cluster_significance.txt" ]; then
##     echo -e "ImmVar, ${cell}"
##     ~/tools/leafcutter/scripts/leafcutter_ds_modified.R \
##     --num_threads 30 \
##     --output_prefix /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_ImmVar_${cell} \
##     --exon_file /data02/home/juninamo/tools/leafcutter/annotation_codes/gencode38/gencode38_all_exons.txt.gz \
##     /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/all_perind_numers.counts.gz \
##     /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/ImmVar_${cell}_groups.txt
##     fi

##     if [ ! -e "/data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_ImmVar_${cell}_leafviz.RData" ]; then
##     echo -e "ImmVar, ${cell}"
##     ~/tools/leafcutter/leafviz/prepare_results_modified.R \
##     /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/all_perind_numers.counts.gz \
##     /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_ImmVar_${cell}_cluster_significance.txt \
##     /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_ImmVar_${cell}_effect_sizes.txt \
##     /data02/home/juninamo/tools/leafcutter/annotation_codes/gencode38/gencode38 \
##     --meta_data_file /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/ImmVar_${cell}_groups.txt \
##     --code ${cell}_vs_others \
##     --num_threads 30 \
##     --ntop_pca 500 \
##     --output /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_ImmVar_${cell}_leafviz.RData
##     fi
## done

# EGA
for cell in NS LPS Pam3CSK4 R848 IAV; do
    if [ ! -e "/data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_EGA_${cell}_cluster_significance.txt" ]; then
    echo -e "EGA, ${cell}"
    ~/tools/leafcutter/scripts/leafcutter_ds_modified.R \
    --num_threads 30 \
    --output_prefix /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_EGA_${cell} \
    --exon_file /data02/home/juninamo/tools/leafcutter/annotation_codes/gencode38/gencode38_all_exons.txt.gz \
    /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/all_perind_numers.counts.gz \
    /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/EGA_${cell}_groups.txt
    fi

    if [ ! -e "/data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_EGA_${cell}_leafviz.RData" ]; then
    echo -e "EGA, ${cell}"
    ~/tools/leafcutter/leafviz/prepare_results_modified.R \
    /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/all_perind_numers.counts.gz \
    /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_EGA_${cell}_cluster_significance.txt \
    /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_EGA_${cell}_effect_sizes.txt \
    /data02/home/juninamo/tools/leafcutter/annotation_codes/gencode38/gencode38 \
    --meta_data_file /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/EGA_${cell}_groups.txt \
    --code ${cell}_vs_others \
    --num_threads 10 \
    --ntop_pca 500 \
    --output /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_EGA_${cell}_leafviz.RData
    fi
done


# DICE
cat /data01/DICE/inamo/cell_type.list | while read cell; do
    if [ ! -e "/data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_DICE_${cell}_cluster_significance.txt" ]; then
    echo -e "DICE, ${cell}"
    ~/tools/leafcutter/scripts/leafcutter_ds_modified.R \
    --num_threads 30 \
    --output_prefix /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_DICE_${cell} \
    --exon_file /data02/home/juninamo/tools/leafcutter/annotation_codes/gencode38/gencode38_all_exons.txt.gz \
    /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/all_perind_numers.counts.gz \
    /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/DICE_${cell}_groups.txt
    fi

    if [ ! -e "/data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_DICE_${cell}_leafviz.RData" ]; then
    echo -e "DICE, ${cell}"
    ~/tools/leafcutter/leafviz/prepare_results_modified.R \
    /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/all_perind_numers.counts.gz \
    /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_DICE_${cell}_cluster_significance.txt \
    /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_DICE_${cell}_effect_sizes.txt \
    /data02/home/juninamo/tools/leafcutter/annotation_codes/gencode38/gencode38 \
    --meta_data_file /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/DICE_${cell}_groups.txt \
    --code ${cell}_vs_others \
    --num_threads 30 \
    --ntop_pca 500 \
    --output /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_DICE_${cell}_leafviz.RData
    fi
done

## # GEUV
## if [ ! -e "/data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_GEUV_cluster_significance.txt" ]; then
## echo -e "GEUV"
## ~/tools/leafcutter/scripts/leafcutter_ds_modified.R \
## --num_threads 30 \
## --output_prefix /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_GEUV \
## --exon_file /data02/home/juninamo/tools/leafcutter/annotation_codes/gencode38/gencode38_all_exons.txt.gz \
## /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/all_perind_numers.counts.gz \
## /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/GEUV_CEU1YRI0_groups.txt
## fi

## if [ ! -e "/data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_GEUV_CEU1YRI0_leafviz.RData" ]; then
## echo -e "ImmVar, ${cell}"
## ~/tools/leafcutter/leafviz/prepare_results_modified.R \
## /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/all_perind_numers.counts.gz \
## /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_GEUV_cluster_significance.txt \
## /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_GEUV_effect_sizes.txt \
## /data02/home/juninamo/tools/leafcutter/annotation_codes/gencode38/gencode38 \
## --meta_data_file /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/GEUV_CEU1YRI0_groups.txt \
## --code CEU_vs_YRI \
## --num_threads 30 \
## --ntop_pca 500 \
## --output /data02/home/juninamo/reference/GENCODE/GRCh38/release_38/leafcutter/leafcutter_ds_GEUV_leafviz.RData
## fi



