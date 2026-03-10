library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(ggbreak)

#####################################################################
# if want to try fisher with 'less' argument, edit in line 98       #
# fisher <- fisher.test(contingency_table, alternative = "greater") #
#####################################################################


TE <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10/TAIR10 transposable elements/TAIR10_Transposable_Elements.txt", sep = "\t")
TE$Transposon_ID <- TE$Transposon_Name
for (chr.i in 1:5) {
    TE$Transposon_Name <- gsub(paste0("AT", chr.i, "TE.*"), paste0("Chr", chr.i), TE$Transposon_Name)
}

filter_fun <- function(db, fam.name, fam.col, is.grep = F) {
    if (is.grep) {
        x <- db[grep(fam.name, fam.col), ]
    } else {
        x <- db[fam.col == fam.name, ]
    }
    return(x)
}


## this function claculate the enrichment pvalue and score for unique IDs in familly copmare to super-family
## 'use.grep' - how to filter.. if by exact name or partial matching
TE_DMRs_overlap <- function(TE_family, TE_super_family, TE_df = TE, is.tairs = F, use.grep = F) {
    ###############################
    ###### TEs
    names(TE_df)[1:4] <- c("seqnames", "strand", "start", "end")
    # TE_df$strand = ifelse(TE_df$strand == "true","+","-")
    TE_df$strand <- "*"
    # TE_gr = makeGRangesFromDataFrame(TE_df, keep.extra.columns = T)

    if (is.tairs) {
        TE_tair_family <- filter_fun(TE_df, TE_family, TE_df$Transposon_ID, is.grep = T)
        TE_tair_super_family <- filter_fun(TE_df, TE_super_family, TE_df$Transposon_Super_Family, is.grep = T)
    } else if (TE_super_family == "all") { # if its test for super-family group compare to all groups
        TE_tair_family <- filter_fun(TE_df, TE_family, TE_df$Transposon_Super_Family, is.grep = T)
        TE_tair_super_family <- TE_df
    } else {
        TE_tair_family <- filter_fun(TE_df, TE_family, TE_df$Transposon_Family, is.grep = use.grep)
        TE_tair_super_family <- filter_fun(TE_df, TE_super_family, TE_df$Transposon_Super_Family, is.grep = T)
    }

    ###############################
    ###### DMRs
    CG <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CG/Transposable_Elements_CG_genom_annotations.csv")
    CHG <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHG/Transposable_Elements_CHG_genom_annotations.csv")
    CHH <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHH/Transposable_Elements_CHH_genom_annotations.csv")

    DMRs <- rbind(CG, CHG, CHH) %>%
        dplyr::rename(Transposon_Name = gene_id) %>%
        mutate(tmp_col = paste(seqnames, start, end, Transposon_Name, sep = "XXX")) %>%
        distinct(tmp_col, .keep_all = T) %>%
        select(-tmp_col)

    if (TE_super_family != "all") {
        if (is.tairs) {
            family_overlap <- filter_fun(DMRs, TE_family, DMRs$Transposon_Name, is.grep = T) # DMRs overlap with family
            super_family_overlap <- filter_fun(DMRs, TE_super_family, DMRs$Transposon_Super_Family, is.grep = T) # DMRs overlap with super-family
        } else {
            # if its test for family group compare to super-family
            family_overlap <- filter_fun(DMRs, TE_family, DMRs$Transposon_Family, is.grep = use.grep) # DMRs overlap with family
            super_family_overlap <- filter_fun(DMRs, TE_super_family, DMRs$Transposon_Super_Family, is.grep = T) # DMRs overlap with super-family
        }
    } else {
        # if its test for super-family group compare to all groups
        family_overlap <- filter_fun(DMRs, TE_family, DMRs$Transposon_Super_Family, is.grep = T) # DMRs overlap with super-family
        super_family_overlap <- DMRs # all DMRs
    }


    ###############################
    ###### enrichment analysis
    # https://doi.org/10.1016/j.ygeno.2017.05.008
    # https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html

    # op2 - 'family' compare to 'super-family'
    a <- length(unique(family_overlap$Transposon_Name)) # overlapped in family
    b <- nrow(TE_tair_family) - a # non-overlapped in family
    c.0 <- length(unique(super_family_overlap$Transposon_Name)) # overlapped in super-family
    d.0 <- nrow(TE_tair_super_family) - c.0 # non-overlapped in super-family
    c <- c.0 - a # without family
    d <- d.0 - b # without family

    contingency_table <- matrix(c(a, c, b, d),
        nrow = 2,
        dimnames = list(
            c(TE_family, TE_super_family),
            c("Overlapped_IDs", "Non-overlapped_IDs")
        )
    )

    fisher <- fisher.test(contingency_table, alternative = "greater")
    # TE_score = as.numeric(fisher$estimate)
    TE_score <- (a / b) / (c / d)


    ###############################
    ###### "gain" or "loss" (direction)  data frame
    direction_family_df <- data.frame(
        hyper = nrow(family_overlap[family_overlap$direction == 1, ]),
        hypo = nrow(family_overlap[family_overlap$direction == -1, ])
    )
    direction_super_family_df <- data.frame(
        hyper = nrow(super_family_overlap[super_family_overlap$direction == 1, ]),
        hypo = nrow(super_family_overlap[super_family_overlap$direction == -1, ])
    )

    ###############################
    return(list(
        overlapped = nrow(family_overlap),
        overlapped_IDs = a,
        annotated = nrow(TE_tair_family),
        pValue = fisher$p.value,
        score = TE_score,
        contingency_table = contingency_table,
        direction_family = direction_family_df,
        direction_super_family = direction_super_family_df
    ))
}

#####################################################
#####################################################

#####################################################
#####################################################
###### all families within each super-family
{
    family_tests <- function(superFamily_name, te.df = TE) {
        unique.f <- unique(te.df[grep(superFamily_name, te.df$Transposon_Super_Family), "Transposon_Family"])
        x <- data.frame(family = NA, pValue = NA, score = NA)
        for (i in 1:length(unique.f)) {
            x[i, 1] <- unique.f[i]
            try({
                x[i, 2] <- TE_DMRs_overlap(unique.f[i], superFamily_name, te.df)$pValue
            })
            try({
                x[i, 3] <- TE_DMRs_overlap(unique.f[i], superFamily_name, te.df)$score
            })
        }
        x <- x %>%
            arrange(desc(score)) %>%
            filter(pValue < 0.05)
        return(x)
    }

    ## plot
    for (family.loop in c("Gypsy", "Copia", "LINE", "SINE|Rath", "DNA", "Helitron")) {
        plot_df <- family_tests(family.loop)
        if (family.loop == "DNA") {
            family.loop <- "TIR"
            height.l <- 5
        } else if (family.loop == "LINE") {
            height.l <- 1.5
        } else {
            height.l <- 2.75
        }

        if (nrow(plot_df) != 0) {
            plot_df$family <- factor(plot_df$family, levels = plot_df$family)
            f.p <- ggplot(plot_df, aes(x = score, y = family)) +
                geom_bar(stat = "identity", width = 0.5, fill = "gray60", colour = "black") +
                theme_classic() +
                theme( # panel.spacing = unit(2, "lines"),
                    text = element_text(family = "serif"),
                    title = element_text(face = "bold"),
                    axis.text.y = element_text(face = "bold", hjust = 0.9, size = 9),
                    axis.text.x = element_text(face = "bold", size = 9),
                    axis.title.x = element_text(size = 12, face = "bold"),
                    axis.line = element_blank(),
                    axis.ticks = element_line(size = 0.75),
                    axis.ticks.length.y = unit(-0.3, "cm"),
                    axis.ticks.length.x = unit(-0.1, "cm")
                    # strip.text = element_blank(),
                    # axis.text.x = element_blank(),
                    # axis.ticks.x = element_blank()
                    # panel.border = element_rect(colour = "black", fill=NA, size=1)
                ) +
                # scale_y_break(c(8.25,34), scales = "fixed", ticklabels = c(0,2,4,6,8,35,36)) +
                labs(y = "", x = "Score", title = gsub("\\|.*", "", family.loop)) +
                geom_text(aes(label = ifelse(pValue < 0.05, ifelse(pValue < 0.01, ifelse(pValue < 0.001, "***", "**"), "*"), "")),
                    hjust = -0.3, size = 4.5
                ) +
                geom_rect(aes(xmin = 0, xmax = max(score) * 1.175, ymin = 0, ymax = nrow(plot_df) + 0.5),
                    fill = "transparent", color = "black", size = 0.75
                ) +
                scale_x_continuous(breaks = seq(0, round(max(plot_df$score)), by = 1))


            svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/", family.loop, "_families_enrichment.svg"),
                width = 3.5, height = height.l
            )
            print(f.p)
            dev.off()
        } else {
            message("*** there is non sig. in '", family.loop, "' super-family")
        }
    }
}

#####################################################
#####################################################
###### each super-families compare to all
{
    families <- list(
        Gypsy = TE_DMRs_overlap("Gypsy", "all"),
        Copia = TE_DMRs_overlap("Copia", "all"),
        LINE = TE_DMRs_overlap("LINE", "all"),
        SINE = TE_DMRs_overlap("SINE|Rath", "all"),
        TIR = TE_DMRs_overlap("DNA", "all"),
        Helitron = TE_DMRs_overlap("Helitron", "all"),
        Unassigned = TE_DMRs_overlap("Unassigned", "all")
    )

    # Extract the relevant information for each family using lapply
    total_sp <- sapply(families, function(fam) fam$annotated)
    total_DMRs_sp <- sapply(families, function(fam) fam$overlapped)
    total_DMRs <- sapply(families, function(fam) (fam$direction_super_family$hyper + fam$direction_super_family$hypo))
    hyper_sp <- sapply(families, function(fam) fam$direction_family$hyper)
    hypo_sp <- sapply(families, function(fam) fam$direction_family$hypo)
    pValue_sp <- sapply(families, function(fam) fam$pValue)
    score_sp <- sapply(families, function(fam) fam$score)

    n.rep <- length(families)
    plot_2_df <- data.frame(
        family = rep(names(families), 2),
        direction = c(rep("Hyper", n.rep), rep("Hypo", n.rep)),
        presentage = c(hyper_sp / total_DMRs, hypo_sp / total_DMRs) * 100,
        pValue = c(pValue_sp, rep(1, n.rep)),
        score = score_sp
    )

    level_order <- unique(plot_2_df[order((hyper_sp + hypo_sp), decreasing = T), 1])
    sf.p <- ggplot(plot_2_df, aes(factor(family, levels = level_order), y = presentage, fill = direction)) +
        geom_bar(stat = "identity", width = 0.5, colour = "black") +
        scale_fill_manual(values = c("Hyper" = "#d96c6c", "Hypo" = "#6c96d9")) +
        theme_classic() +
        theme( # panel.spacing = unit(2, "lines"),
            legend.key.size = unit(0.4, "cm"),
            text = element_text(family = "serif"),
            axis.text.x = element_text(angle = 45, vjust = 0.75, hjust = 0.5, size = 9), # , face="bold"),
            axis.text.y = element_text(size = 10), # , face="bold"),
            axis.title.y = element_text(size = 12), # , face="bold"),
            axis.line = element_blank(),
            axis.ticks = element_line(color = "black", size = 0.75),
            axis.ticks.length.x = unit(-0.225, "cm"),
            axis.ticks.length.y = unit(-0.05, "cm")
        ) +
        scale_y_continuous(breaks = seq(0, round(max(plot_2_df$presentage)) * 1.1, by = 5)) +
        labs(y = "DMRs Count (%)", x = "", fill = "Direction") +
        geom_text(
            aes(
                y = c(((hyper_sp + hypo_sp) / total_DMRs) * 100, rep(0, n.rep)),
                label = ifelse(pValue < 0.05, ifelse(pValue < 0.01, ifelse(pValue < 0.001, "***", "**"), "*"), "")
            ),
            vjust = 0.25, size = 4
        ) +
        geom_rect(aes(xmin = 0.5, xmax = n.rep + 0.5, ymin = 0, ymax = round(max(presentage)) * 1.1),
            fill = "transparent", color = "black", size = 1
        )

    svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/Super-family_enrichment.svg"),
        width = 3, height = 3
    )
    print(sf.p)
    dev.off()


    ################################################
    ###  prepare table for supplementary
    ### run untill 'plot_2_df' in the above script
    super_families_results <- rbind(
        do.call(rbind, lapply(families, function(x) {
            data.frame(Overlapped_IDs = x$contingency_table[1, 1], Non_overlapped_IDs = x$contingency_table[1, 2])
        }))
    ) %>%
        cbind(plot_2_df[1:7, ]) %>%
        relocate(family, .before = Overlapped_IDs) %>%
        relocate(pValue, .after = Non_overlapped_IDs) %>%
        select(-direction) %>%
        rename(hyper_to_total = presentage) %>%
        mutate(hypo_to_total = plot_2_df[8:14, 3],
        total_DMRs = total_DMRs_sp,
        hyper = (hyper_sp / total_DMRs_sp) * 100,
        hypo = (hypo_sp / total_DMRs_sp) * 100
        ) %>%
        relocate(total_DMRs, .after = pValue) %>%
        relocate(hyper, .after = total_DMRs) %>%
        relocate(hypo, .after = hyper) %>%
        relocate(hypo_to_total, .after = hyper_to_total) %>%
        relocate(score, .before = pValue)

    # total TEs row
    all_TE_df <- rbind(
        read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CG/Transposable_Elements_CG_genom_annotations.csv"),
        read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHG/Transposable_Elements_CHG_genom_annotations.csv"),
        read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHH/Transposable_Elements_CHH_genom_annotations.csv")
    ) %>%
        mutate(tmp = paste(seqnames, start, end, gene_id, sep = "_")) %>%
        distinct(tmp, .keep_all = T)

    all_overlap <- all_TE_df %>%
        distinct(gene_id) %>%
        nrow()
    all_non_overlap <- (unique(TE$Transposon_ID) %>% length()) - all_overlap

    all_hyper <- (families[[1]]$direction_super_family["hyper"] / nrow(all_TE_df)) * 100
    all_hypo <- (families[[1]]$direction_super_family["hypo"] / nrow(all_TE_df)) * 100

    write.csv(
        rbind(
            super_families_results,
            data.frame(family = "Total TEs", Overlapped_IDs = all_overlap, Non_overlapped_IDs = all_non_overlap, pValue = "", score = "", total_DMRs = unique(total_DMRs), hyper = all_hyper, hypo = all_hypo, hyper_to_total = "", hypo_to_total = "")
        ),
        "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/Super-family_enrichment_results.csv",
        row.names = F
    )
}










#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
###### LTR retrotransposon groups
## ONSEN (DOI: 10.1038/nature09861)
ONSEN_df <- data.frame(
    Name = paste0("ONSEN", 1:8),
    gene_id = c("AT1G11265", "AT3G61330", "AT5G13205", "AT1G58140", "AT1G48710", "AT3G59720", "AT1G21945", "AT3G32415"),
    Transposon_ID = c("AT1TE12295", "AT3TE92525", "AT5TE15240", "AT1TE71045", "AT1TE59755", "AT3TE89830", "AT1TE24850", "AT3TE54550")
)
ONSEN <- TE_DMRs_overlap(paste(ONSEN_df$Transposon_ID, collapse = "|"), "LTR/Copia", is.tairs = T)

## EVADE (DOI: 10.1038/nature08328)
evd <- TE_DMRs_overlap("AT5TE20395", "LTR/Copia", is.tairs = T) # maybe plot the DMR-pValue

## SISYPHUS
sis <- TE_DMRs_overlap("AT3TE76225", "LTR/Copia", is.tairs = T)

## ATHILA2
ATHILA2 <- TE_DMRs_overlap("ATHILA2", "LTR/Gypsy")

## 'TE_DMRs_overlap' function results for groups
LTR_DMRs_res <- data.frame(
    Name = c("ONSEN", "EVADE", "SISYPHUS", "ATHILA2"),
    IDs = c(ONSEN$annotated, evd$annotated, sis$annotated, ATHILA2$annotated),
    Overlapped.IDs = c(ONSEN$overlapped_IDs, evd$overlapped_IDs, sis$overlapped_IDs, ATHILA2$overlapped_IDs),
    sig.TEGs = NA,
    Overlapped.DMRs = c(ONSEN$overlapped, evd$overlapped, sis$overlapped, ATHILA2$overlapped)
)


## df for RNA
sig_LTR_0 <- data.frame(
    Name = c("EVADE", "SISYPHUS"),
    gene_id = c("AT5G17125", "AT3G50625"),
    Transposon_ID = c("AT5TE20395", "AT3TE76225")
) %>%
    rbind(ONSEN_df, .)

###########
## merge with RNAseq
RNA <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    rename(gene_id = locus_tag) %>%
    select(gene_id, log2FoldChange, pValue)

TEG <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/Methylome.At_description_file.csv.gz") %>%
    select(gene_id, Derives_from) %>%
    filter(!is.na(Derives_from)) %>%
    rename(Transposon_ID = Derives_from) %>%
    merge.data.frame(., RNA, by = "gene_id") %>%
    merge.data.frame(., TE[, c("Transposon_Family", "Transposon_Super_Family", "Transposon_ID")]) %>%
    arrange(pValue)



## LTR-TEGs data frame
sig_LTR <- data.frame(
    Name = "ATHILA2",
    gene_id = TEG[TEG$Transposon_Family == "ATHILA2", "gene_id"],
    Transposon_ID = TEG[TEG$Transposon_Family == "ATHILA2", "Transposon_ID"]
) %>%
    rbind(sig_LTR_0, .) %>%
    merge.data.frame(., TEG[, -grep("gene_id", names(TEG))], by = "Transposon_ID") %>%
    arrange(pValue)

TEGs_group <- function(x) {
    ann <- sig_LTR[grep(x, sig_LTR$Name), ]
    sig <- ann[ann$pValue < 0.05, ]
    xx <- paste0(nrow(sig), "/", nrow(ann))
    return(xx)
}
LTR_DMRs_res$sig.TEGs <- c(
    TEGs_group("ONSEN"),
    TEGs_group("EVADE"),
    TEGs_group("SISYPHUS"),
    TEGs_group("ATHILA2")
)
