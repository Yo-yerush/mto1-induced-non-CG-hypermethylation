library(dplyr)
library(ggplot2)
library(cowplot)
library(grid)
library(GenomicRanges)

    ###### TEs results
    TE_1 <- rbind(
        read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CG/Transposable_Elements_CG_genom_annotations.csv"),
        read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHG/Transposable_Elements_CHG_genom_annotations.csv"),
        read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHH/Transposable_Elements_CHH_genom_annotations.csv")
    )
    #dplyr::rename(gene_id = Transposon_Name) %>%
    #if (TE_SF_loop != "all") {
    #    TE_1 <- TE_1 %>% filter(grepl(TE_SF_loop, Transposon_Super_Family))
    #}
    TE_1 <- makeGRangesFromDataFrame(TE_1, keep.extra.columns = T) %>%
        sort()


    ################################
    ###### DEGs
    DEGs <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv")
    names(DEGs) <- gsub("locus_tag", "gene_id", names(DEGs))
    symbol_indx <- DEGs %>%
        dplyr::select(gene_id, gene, log2FoldChange, pValue, short_description) %>%
        dplyr::rename(Symbol = gene)
    symbol_indx$log2FoldChange <- round(symbol_indx$log2FoldChange, 3)

    upregulated <- DEGs %>%
        filter(log2FoldChange > 0 & padj < 0.05) %>% # & gene_model_type == "protein_coding") %>%
        dplyr::select(gene_id)

    downregulated <- DEGs %>%
        filter(log2FoldChange < 0 & padj < 0.05) %>% # gene_model_type == "protein_coding") %>%
        dplyr::select(gene_id)

    non_sig <- DEGs %>%
        filter(padj > 0.5 & gene_model_type == "protein_coding") %>%
        dplyr::select(gene_id)


    ################################
    ###### TAIR10 annotations
    gff3 <- rtracklayer::import.gff3("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10/TAIR10 gff3/TAIR10_GFF3_genes.gff")
    gff3_df <- gff3 %>%
        as.data.frame() %>%
        dplyr::rename(gene_id = ID) %>%
        filter(type == "gene") %>%
        dplyr::select(seqnames, start, end, width, strand, gene_id)


    ################################
    ################################
    ###### overlapp with TEs up/down-stream

    # overlap with TEs function
    overlap_fun <- function(x.gr, x.TE = TE_1) {
        m <- findOverlaps(x.gr, x.TE)
        x.out <- x.gr[queryHits(m)] %>%
            as.data.frame() %>%
            mutate(tmp = paste(.$seqnames, .$start, .$end, .$strand, .$gene_id, sep = "_")) %>%
            distinct(tmp, .keep_all = T) %>%
            dplyr::select(gene_id)

        # each gene and its TEs overlapped family
        TE2gene_df <- data.frame(
            gene_id = x.gr[queryHits(m)]$gene_id,
            TE_family = x.TE[subjectHits(m)]$Transposon_Family,
            TE_super_family = x.TE[subjectHits(m)]$Transposon_Super_Family
        ) %>%
            mutate(tmp = paste(.$gene_id, .$TE_family, .$TE_super_family, sep = "_del_")) %>%
            distinct(tmp, .keep_all = T) %>%
            dplyr::select(-tmp)

        return(list(
            overlapIDs = x.out,
            TE2gene = TE2gene_df
        ))
    }

    # remove suplicates fanilies or super-families in columns
    remove_dup <- function(y) {
        y <- as.character(unique(unlist(strsplit(y, "; "))))
        paste(y, collapse = "; ")
    }

    #######
    ## main function
    DEG_near_TEs <- function(deg_type, TE.upstream = F, TE.downstream = F, TE.both_sides = F) {
        deg.gr <- merge(deg_type, gff3_df) %>%
            relocate(gene_id, .after = last_col()) %>%
            makeGRangesFromDataFrame(., keep.extra.columns = T)
        upstream.l <- shift(deg.gr, 2000)
        downstream.l <- shift(deg.gr, -2000)

        up_overlap <- overlap_fun(upstream.l)$overlapIDs
        down_overlap <- overlap_fun(downstream.l)$overlapIDs

        upstream_TE_fam <- overlap_fun(upstream.l)$TE2gene
        downstream_TE_fam <- overlap_fun(downstream.l)$TE2gene
        names(upstream_TE_fam) <- gsub("TE", "upstream", names(upstream_TE_fam))
        names(downstream_TE_fam) <- gsub("TE", "downstream", names(downstream_TE_fam))

        if (TE.both_sides) {
            return_df <- merge(up_overlap, down_overlap) %>%
                merge(., upstream_TE_fam) %>%
                merge(., downstream_TE_fam)
        } else if (TE.upstream) {
            return_df <- up_overlap %>%
                merge(., upstream_TE_fam)
        } else if (TE.downstream) {
            return_df <- down_overlap %>%
                merge(., downstream_TE_fam)
        }

        # group families and super-families by 'gene_id'
        return_df <- return_df %>%
            #mutate(tmp = paste(return_df[, 1], return_df[, 2], return_df[, 3], sep = "_")) %>%
            group_by(gene_id) %>%
            summarise(
                across(
                    contains("stream_family") | contains("stream_super_family"),
                    ~ remove_dup(paste(., collapse = "; "))
                )
            ) %>%
            merge(symbol_indx, ., by = "gene_id") %>%
            as.data.frame() %>%
            #dplyr::select(-tmp) %>%
            arrange(pValue)

        return(return_df)
    }

    xl_list <- list(
        Up.Stream = rbind(
            DEG_near_TEs(upregulated, TE.upstream = T),
            DEG_near_TEs(downregulated, TE.upstream = T),
            DEG_near_TEs(non_sig, TE.upstream = T)
        ),
        Down.Stream = rbind(
            DEG_near_TEs(upregulated, TE.downstream = T),
            DEG_near_TEs(downregulated, TE.downstream = T),
            DEG_near_TEs(non_sig, TE.downstream = T)
        ) # ,
        # Both.Sides = rbind(
        #    DEG_near_TEs(upregulated, TE.both_sides = T),
        #    DEG_near_TEs(downregulated, TE.both_sides = T)
        # )
    )

    ########################

    ########################

    final_df = data.frame(TE_superfamily = c("Gypsy", "Copia", "LINE", "Helitron", "DNA", "SINE"), Upregulated = NA, Downregulated = NA, Non.DE = NA)

    for (TE_SF_loop in final_df$TE_superfamily) {
        ########################
        # upstream <- readxl::read_excel("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/TEs_near_DEGs/TEs_near_DEGs_201124.xlsx", sheet = "Up.Stream")

        upstream <- xl_list[["Up.Stream"]]
        if (TE_SF_loop != "all") {
            upstream <- upstream %>% filter(grepl(TE_SF_loop, upstream_super_family))
        }
        upstream <- upstream %>% dplyr::select(gene_id)

        # downstream <- readxl::read_excel("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/TEs_near_DEGs/TEs_near_DEGs_201124.xlsx", sheet = "Down.Stream")
        downstream <- xl_list[["Down.Stream"]]
        if (TE_SF_loop != "all") {
            downstream <- downstream %>% filter(grepl(TE_SF_loop, downstream_super_family))
        }
        downstream <- downstream %>% dplyr::select(gene_id)
        ########################


        RNAseq <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
            dplyr::rename(gene_id = locus_tag, DEG_log2FC = log2FoldChange) %>%
            dplyr::select(gene_id, DEG_log2FC, padj, pValue)

        overlapped_TE_context <- function(context, TE_SF) {
            TE <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", context, "/Transposable_Elements_", context, "_genom_annotations.csv"))
            if (TE_SF != "all") {
                TE <- TE %>% filter(grepl(TE_SF, Transposon_Super_Family))
            }
            TE <- makeGRangesFromDataFrame(TE, keep.extra.columns = TRUE)

            overlapped <- data.frame()

            gene_feature <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", context, "/Genes_", context, "_genom_annotations.csv")) %>%
                makeGRangesFromDataFrame(keep.extra.columns = TRUE)

            m <- findOverlaps(gene_feature, TE)
            m_df <- gene_feature[queryHits(m)]
            m_df$overlapping_TE <- paste(TE$Transposon_Super_Family[subjectHits(m)], TE$Transposon_Family[subjectHits(m)], sep = ", ")
            # m_df$Transposon_Super_Family <- TE$Transposon_Super_Family[subjectHits(m)]
            # m_df$Transposon_Family <- TE$Transposon_Family[subjectHits(m)]

            overlapped <- rbind(overlapped, as.data.frame(m_df)) %>%
                dplyr::select(gene_id)

            return(overlapped)
        }

        ################
        remove_dup_DMR <- function(y) {
            y <- as.character(unique(unlist(strsplit(y, ","))))
            paste(y, collapse = ",")
        }
        ################

        ################
        overlapped_TE_sig = rbind(
            overlapped_TE_context("CG", TE_SF_loop),
            overlapped_TE_context("CHG", TE_SF_loop),
            overlapped_TE_context("CHH", TE_SF_loop),
            upstream,
            downstream
        ) %>%
            distinct(gene_id) %>%
            merge.data.frame(RNAseq, ., by = "gene_id") %>%
            #filter(padj < 0.05) %>%
            dplyr::select(gene_id, DEG_log2FC, padj)
        
        final_df[final_df$TE_superfamily == TE_SF_loop, "Upregulated"] = filter(overlapped_TE_sig, DEG_log2FC > 0) %>% filter(padj < 0.05) %>% nrow()

        final_df[final_df$TE_superfamily == TE_SF_loop, "Downregulated"] = filter(overlapped_TE_sig, DEG_log2FC < 0) %>% filter(padj < 0.05) %>% nrow()
        
        final_df[final_df$TE_superfamily == TE_SF_loop, "Non.DE"] = filter(overlapped_TE_sig, padj > 0.05) %>% nrow()
    }

write.csv(final_df, "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/Overlapped_TEs_with_DMRs_n_DEGs_summarized.csv", row.names = F)
