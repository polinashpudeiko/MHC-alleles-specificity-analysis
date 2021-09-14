library(seqinr)
library(immunotation)
library(RCurl)
library(ggpointdensity)
library(bios2mds)
library(parallel)
library(ggplot2)
library(reshape2)



# Processing of alignments ------------------------------------------------

msf.res  <- read.alignment(file = "./A_prot.msf",
                           format = "msf")

# drop sequences started with N
msf.res$seq <- msf.res$seq[!grepl("N",msf.res$nam)]
msf.res$nam <- msf.res$nam[!grepl("N",msf.res$nam)]

# drop duplicates
msf.res$nam <- msf.res$nam[!duplicated(msf.res$seq)]
msf.res$seq <- msf.res$seq[!duplicated(msf.res$seq)]

# only longest sequence with different mutations left
names.subset <- unlist(purrr::map(purrr::map(msf.res$nam, 
                                             .f = ~strsplit(.x, ':')), 
                                  .f = ~paste(.x[[1]][1],.x[[1]][2], sep=':')), 
                       use.names = FALSE)
msf.res$nam <- msf.res$nam[!duplicated(names.subset)]
msf.res$seq <- msf.res$seq[!duplicated(names.subset)]

# approach with cutting sequence
msf.res$seq <- purrr::map(msf.res$seq, .f = ~substring(.x, first=28, last=115))
msf.res$nb <- 4064

###
###
###
# Identity/Similarity estimation
###
###
###

# it's possible to choose matrix = "similarity" for similarity estimation
msf.alignment <-  as.matrix(dist.alignment(msf.res, matrix = "identity" ))

# write.csv(msf.alignment, 'identity_hlaa_longest_cutted_a1.csv')

###
###
###
# BLOSUM/PAM usage
###
###
###

msf_al_bios <- import.msf("./A_prot.msf", 
                          aa.to.upper = TRUE, gap.to.dash = TRUE)

msf.res  <- read.alignment(file = "./A_prot.msf",
                           format = "msf")

msf_al_bios_sub <- msf_al_bios[!grepl("N",msf.res$nam)]

msf.res$seq <- msf.res$seq[!grepl("N",msf.res$nam)]
msf.res$nam <- msf.res$nam[!grepl("N",msf.res$nam)]

msf_al_bios_sub <- msf_al_bios_sub[!duplicated(msf.res$seq)]

alleles_comb <- t(combn(unique(names(msf_al_bios_sub)), 2))

numCores <- detectCores()

#PAM
dissim_for_pair_pam <- function(n) {
  dissim <- dis(msf_al_bios_sub[[alleles_comb[n, 1]]], 
                msf_al_bios_sub[[alleles_comb[n, 2]]], sub.mat.id = "PAM40", 
                gap = NULL)
  return(dissim)
}

results_pam <- mclapply(1:9101511, dissim_for_pair_pam, mc.cores = numCores)

alleles_comb_res_pam <- data.frame(alleles_comb)
alleles_comb_res_pam$results_pam <- results_pam

# dealing with problems from unmelting table
alleles_comb_res_df <- dcast(alleles_comb_res_pam, X1 ~ X2)
rownames(alleles_comb_res_df) <- alleles_comb_res_df$X1
row <- alleles_comb_res_df[c(1),-c(1)]
alleles_comb_res_df$X1 <- t(row)

row_to_add <- alleles_comb_res_df[["A*80:07"]]
alleles_comb_res_df <- data.frame(t(alleles_comb_res_df))

alleles_comb_res_df[["A*80:07"]] <- c(row_to_add, 0)
alleles_comb_res_df <- data.frame(t(alleles_comb_res_df))

alleles_comb_res_df[lower.tri(alleles_comb_res_df)] <- t(alleles_comb_res_df)[lower.tri(alleles_comb_res_df)]
alleles_comb_res_df[is.null(alleles_comb_res_df)] <- 0
colnames(alleles_comb_res_df) <- rownames(alleles_comb_res_df)

# write.csv(as.matrix(alleles_comb_res_df), 'allele_pam40.csv', quote = F)

#JTT
alleles_comb_res <- data.frame(alleles_comb)
results <- mclapply(1:9101511, dissim_for_pair, mc.cores = numCores)
alleles_comb_res$results_jtt <- results

# dealing with problems from unmelting table
alleles_comb_res_df <- dcast(alleles_comb_res, X1 ~ X2)
rownames(alleles_comb_res_df) <- alleles_comb_res_df$X1
row <- alleles_comb_res_df[c(1),-c(1)]
alleles_comb_res_df$X1 <- t(row)

row_to_add <- alleles_comb_res_df[["A*80:07"]]
alleles_comb_res_df <- data.frame(t(alleles_comb_res_df))

alleles_comb_res_df[["A*80:07"]] <- c(row_to_add, 0)
alleles_comb_res_df <- data.frame(t(alleles_comb_res_df))

alleles_comb_res_df[lower.tri(alleles_comb_res_df)] <- t(alleles_comb_res_df)[lower.tri(alleles_comb_res_df)]
alleles_comb_res_df[is.null(alleles_comb_res_df)] <- 0
colnames(alleles_comb_res_df) <- rownames(alleles_comb_res_df)

# write.csv(as.matrix(alleles_comb_res_df), 'allele_jtt.csv', quote = F)

#BLOSSUM
alleles_comb_res <- data.frame(alleles_comb)
results <- mclapply(1:9101511, dissim_for_pair, mc.cores = numCores)
alleles_comb_res$results_bloss <- unlist(results)

# dealing with problems from unmelting table
alleles_comb_res_df <- dcast(alleles_comb_res, X1 ~ X2)
rownames(alleles_comb_res_df) <- alleles_comb_res_df$X1
row <- alleles_comb_res_df[c(1),-c(1)]
alleles_comb_res_df$X1 <- t(row)

row_to_add <- alleles_comb_res_df[["A*80:07"]]
alleles_comb_res_df <- data.frame(t(alleles_comb_res_df))

alleles_comb_res_df[["A*80:07"]] <- c(row_to_add, 0)
alleles_comb_res_df <- data.frame(t(alleles_comb_res_df))

alleles_comb_res_df[lower.tri(alleles_comb_res_df)] <- t(alleles_comb_res_df)[lower.tri(alleles_comb_res_df)]
alleles_comb_res_df[is.null(alleles_comb_res_df)] <- 0
colnames(alleles_comb_res_df) <- rownames(alleles_comb_res_df)

# write.csv(as.matrix(alleles_comb_res_df), 'allele_blossum80.csv', quote = F)


# Generate script for netMHCpan -------------------------------------------

allele_list <- read.delim("~/Desktop/Internship/HLA-alignment/allele_list.txt", 
                          header=FALSE)

c <- 1
file_bash_script<-file('script_mhcI.sh', 'w')
for (i in seq(1, length(allele_list$V1))){
  if (nchar(paste(allele_list$V1[c:i], collapse=",")) >= 900){
    output_xlsI <- paste(c('mhcI', c, 'out.xls'), collapse = '_')
    writeLines(paste(
      "./netMHCpan -f ../immunotation/inst/extdata/virus_example.fasta -v -BA -xls -a",
      paste(allele_list$V1[c:i-1], collapse=","),"-xlsfile", output_xlsI,sep=' '), 
      file_bash_script)
    c <- i
  }
}
close(file_bash_script)


# Obtaining matrices with BA/EL scores ------------------------------------

# BA-scores

mhcI_conc.xls <- readr::read_delim('mhcI_conc.csv',delim = ",",col_names = TRUE)
mhcI_conc.xls <- mhcI_conc.xls[,-c(1)]

allele_specif_cor_matr_BA <- data.frame(matrix(ncol = length(unique(proper_format)), 
                                               nrow = length(unique(proper_format))))
colnames(allele_specif_cor_matr_BA) <- unique(proper_format)
rownames(allele_specif_cor_matr_BA) <- unique(proper_format)

alleles_comb <- t(combn(unique(proper_format), 2))

for (n in 1:5302396){
  allele1_core <- paste('core', alleles_comb[n, 1], sep='.')
  allele2_core <- paste('core', alleles_comb[n, 2], sep='.')
  m1 <-ifelse(mhcI_conc.xls[[allele1_core]] == mhcI_conc.xls[[allele2_core]], 
              mhcI_conc.xls[[allele1_core]], NA)
  m2 <-ifelse(mhcI_conc.xls[[allele2_core]] == mhcI_conc.xls[[allele1_core]], 
              mhcI_conc.xls[[allele2_core]], NA)
  allele1_score <- paste('BA-score', alleles_comb[n, 1], sep='.')
  allele2_score <- paste('BA-score', alleles_comb[n, 2], sep='.')
  allele_specif_cor_matr_BA[alleles_comb[n, 1], alleles_comb[n, 2]] <-
    cor(mhcI_conc.xls[[allele1_score]][!is.na(m1)], 
        mhcI_conc.xls[[allele2_score]][!is.na(m2)])
  allele_specif_cor_matr_BA[alleles_comb[n, 2], alleles_comb[n, 1]] <-
    cor(mhcI_conc.xls[[allele1_score]][!is.na(m1)], 
        mhcI_conc.xls[[allele2_score]][!is.na(m2)])
}

# write.csv(allele_specif_cor_matr_BA, 'allele_specif_cor_matr_BA.csv', quote = F)

# EL-scores

allele_specif_cor_matr <- data.frame(matrix(ncol = length(unique(proper_format)), 
                                            nrow = length(unique(proper_format))))
colnames(allele_specif_cor_matr) <- unique(proper_format)
rownames(allele_specif_cor_matr) <- unique(proper_format)

alleles_comb <- t(combn(unique(proper_format), 2))

for (n in 1:5302396){
  allele1_core <- paste('core', alleles_comb[n, 1], sep='.')
  allele2_core <- paste('core', alleles_comb[n, 2], sep='.')
  m1 <-ifelse(mhcI_conc.xls[[allele1_core]] == mhcI_conc.xls[[allele2_core]], 
              mhcI_conc.xls[[allele1_core]], NA)
  m2 <-ifelse(mhcI_conc.xls[[allele2_core]] == mhcI_conc.xls[[allele1_core]], 
              mhcI_conc.xls[[allele2_core]], NA)
  allele1_score <- paste('BA-score', alleles_comb[n, 1], sep='.')
  allele2_score <- paste('BA-score', alleles_comb[n, 2], sep='.')
  allele_specif_cor_matr[alleles_comb[n, 1], alleles_comb[n, 2]] <-
    cor(mhcI_conc.xls[[allele1_score]][!is.na(m1)], 
        mhcI_conc.xls[[allele2_score]][!is.na(m2)])
  allele_specif_cor_matr[alleles_comb[n, 2], alleles_comb[n, 1]] <-
    cor(mhcI_conc.xls[[allele1_score]][!is.na(m1)], 
        mhcI_conc.xls[[allele2_score]][!is.na(m2)])
}

# write.csv(allele_specif_cor_matr, 'allele_specif_cor_matr.csv', quote = F)


# Similarity vs.  netMHCpan scores ----------------------------------------

sim_al <- readr::read_delim('similarity_hlaa_wo_dupl.csv',delim = ",",
                            col_names = TRUE)
sim_al <- data.frame(sim_al)
rownames(sim_al) <- sim_al$X1
sim_al <- sim_al[,-c(1)]

first_line_al <- strsplit(readr::read_lines('similarity_hlaa_wo_dupl.csv', 
                                            n_max = 1), '\",\"')[[1]]
colnames(sim_al) <- first_line_al[2:length(first_line_al)]

proper_format <- unlist(purrr::map(first_line_al[2:length(first_line_al)], 
                                   .f = ~get_mhcpan_input(.x, 'MHC-I')), 
                        use.names = FALSE)

transf <- unlist(purrr::map(purrr::map(gsub('A\\*', 'HLA-A', colnames(sim_al)), 
                                       .f = ~strsplit(.x, ':')), 
                            .f = ~paste(.x[[1]][1],.x[[1]][2], sep=':')), 
                 use.names = FALSE)

m <- match(transf, proper_format)

sim_al_sub <- sim_al[!is.na(transf[m]),!is.na(transf[m])]

# write.csv(sim_al_sub, 'similarity_hlaa_only_netmhcpan_proper_names_wo_dupl.csv', quote = F)


sim_al_sub <- readr::read_delim('similarity_hlaa_only_netmhcpan_proper_names_wo_dupl.csv',
                                delim = ",",col_names = TRUE)
sim_al_sub <- data.frame(sim_al_sub)
rownames(sim_al_sub) <- sim_al_sub$X1
sim_al_sub <- sim_al_sub[,-c(1)]

sim_al_sub$format <- transf[!is.na(transf[m])]

# Approach with sum duplicates - the worst 

# v <- aggregate(sim_al_sub[-(3452)], list(sim_al_sub[,3452]), 
#                function(x) sum(unique(x)))
# rownames(v) <- v$Group.1
# v <- v[,-c(1)]
# v_t <- data.frame(t(v))
# v_t$format <- transf[!is.na(transf[m])]
# v_t_agg <- aggregate(v_t[-(3250)], list(v_t[,3250]), function(x) sum(unique(x)))
# rownames(v_t_agg) <- v_t_agg$Group.1
# v_t_agg <- v_t_agg[,-c(1)]
# rownames(v_t_agg) <- gsub('-', '\\.', rownames(v_t_agg))
# rownames(v_t_agg) <- gsub(':', '\\.', rownames(v_t_agg))
# 
# m <- match(colnames(v_t_agg), colnames(allele_specif_cor_matr))
# v_t_agg <- v_t_agg[!is.na(colnames(v_t_agg)[m]),
#                    !is.na(colnames(v_t_agg)[m])]


# Approach with taking the mean for duplicates
v_mean <- aggregate(sim_al_sub[-(3452)], list(sim_al_sub[,3452]), 
                    function(x) mean(unique(x)))
rownames(v_mean) <- v_mean$Group.1
v_mean <- v_mean[,-c(1)]
m <- match(unlist(purrr::map(purrr::map(gsub('A\\*', 'HLA-A', colnames(sim_al)), 
                                        .f = ~strsplit(.x, ':')), 
                             .f = ~paste(.x[[1]][1],.x[[1]][2], sep=':')), 
                  use.names = FALSE), proper_format)
v_mean_t <- data.frame(t(v_mean))
v_mean_t$format <- transf[!is.na(transf[m])]

v_mean_t_t_agg <- aggregate(v_mean_t[-(3250)], list(v_mean_t[,3250]), 
                            function(x) mean(unique(x)))
rownames(v_mean_t_t_agg) <- v_mean_t_t_agg$Group.1
v_mean_t_t_agg <- v_mean_t_t_agg[,-c(1)]
rownames(v_mean_t_t_agg) <- gsub('-', '\\.', rownames(v_mean_t_t_agg))
rownames(v_mean_t_t_agg) <- gsub(':', '\\.', rownames(v_mean_t_t_agg))

###
###
###
# Loading EL-scores
###
###
###

allele_specif_cor_matr <- readr::read_delim('allele_specif_cor_matr.csv',
                                            delim = ",",col_names = TRUE)
allele_specif_cor_matr <- data.frame(allele_specif_cor_matr)
rownames(allele_specif_cor_matr) <- allele_specif_cor_matr$X1
allele_specif_cor_matr <- allele_specif_cor_matr[,-c(1)]
allele_specif_cor_matr[is.na(allele_specif_cor_matr)] <- 1

m <- match(colnames(allele_specif_cor_matr), colnames(v_mean_t_t_agg))
allele_specif_cor_matr <- allele_specif_cor_matr[!is.na(colnames(allele_specif_cor_matr)[m]),
                                                 !is.na(colnames(allele_specif_cor_matr)[m])]
rownames(allele_specif_cor_matr) <- colnames(allele_specif_cor_matr)
rownames(allele_specif_cor_matr) <- gsub('-', '\\.', 
                                         rownames(allele_specif_cor_matr))
rownames(allele_specif_cor_matr) <- gsub(':', '\\.', 
                                         rownames(allele_specif_cor_matr))
allele_specif_cor_matr <- allele_specif_cor_matr[colnames(v_mean_t_t_agg),
                                                 colnames(v_mean_t_t_agg)]


res_cor_mean <- cor(c(as.matrix(allele_specif_cor_matr)), 
                    c(as.matrix(v_mean_t_t_agg)))

allele_specif_cor_matr_m <- melt(replace(allele_specif_cor_matr, 
                                         lower.tri(allele_specif_cor_matr, 
                                                   TRUE), NA), na.rm = TRUE)
v_mean_t_t_agg_m <- melt(replace(v_mean_t_t_agg, 
                                 lower.tri(v_mean_t_t_agg, 
                                           TRUE), NA), na.rm = TRUE)

merged_df_mean <- merge(allele_specif_cor_matr_m, v_mean_t_t_agg_m,by='row.names')

# distribution of alignment's scores/netMHCpan scores
ggplot(merged_df_mean, aes(value.y)) +
  geom_density()+theme_classic()

ggplot(merged_df_mean, aes(value.x, value.y)) +
  geom_hex(bins = 100)+theme_classic()

# ggsave('merged_df_mean_sim_al.png', dpi = 300, width = 6, height = 4)

###
###
###
# Loading BA-scores
###
###
###

allele_specif_cor_matr_BA <- readr::read_delim('allele_specif_cor_matr_BA.csv',
                                               delim = ",",col_names = TRUE)
allele_specif_cor_matr_BA <- data.frame(allele_specif_cor_matr_BA)
rownames(allele_specif_cor_matr_BA) <- allele_specif_cor_matr_BA$X1
allele_specif_cor_matr_BA <- allele_specif_cor_matr_BA[,-c(1)]
allele_specif_cor_matr_BA[is.na(allele_specif_cor_matr_BA)] <- 1

m <- match(colnames(allele_specif_cor_matr_BA), colnames(v_mean_t_t_agg))
allele_specif_cor_matr_BA <- allele_specif_cor_matr_BA[!is.na(colnames(allele_specif_cor_matr_BA)[m]),
                                                       !is.na(colnames(allele_specif_cor_matr_BA)[m])]
rownames(allele_specif_cor_matr_BA) <- colnames(allele_specif_cor_matr_BA)
rownames(allele_specif_cor_matr_BA) <- gsub('-', '\\.', rownames(allele_specif_cor_matr_BA))
rownames(allele_specif_cor_matr_BA) <- gsub(':', '\\.', rownames(allele_specif_cor_matr_BA))
allele_specif_cor_matr_BA <- allele_specif_cor_matr_BA[colnames(v_mean_t_t_agg),
                                                       colnames(v_mean_t_t_agg)]

res_cor_mean <- cor(c(as.matrix(allele_specif_cor_matr_BA)), 
                    c(as.matrix(v_mean_t_t_agg)))

allele_specif_cor_matr_BA_m <- melt(replace(allele_specif_cor_matr_BA, 
                                            lower.tri(allele_specif_cor_matr_BA, 
                                                      TRUE), NA), na.rm = TRUE)
v_mean_t_t_agg_m <- melt(replace(v_mean_t_t_agg, 
                                 lower.tri(v_mean_t_t_agg, 
                                           TRUE), NA), na.rm = TRUE)

merged_df_mean <- merge(allele_specif_cor_matr_BA_m, v_mean_t_t_agg_m,
                        by='row.names')

# distribution of alignment's scores/netMHCpan scores
ggplot(merged_df_mean, aes(value.y)) +
  geom_density()+theme_classic()

ggplot(merged_df_mean, aes(value.x, value.y)) +
  geom_hex(bins = 100)+theme_classic()

# ggsave('merged_df_mean_sim_al.png', dpi = 300, width = 6, height = 4)



# Identity vs.netMHCpan scores -------------------------------------------


id_al <- readr::read_delim('identity_hlaa_longest_cutted.csv',
                           delim = ",", 
                           col_names = TRUE)
id_al <- data.frame(id_al)
rownames(id_al) <- id_al$X1
id_al <- id_al[,-c(1)]

rownames(id_al) <- colnames(id_al)

first_line_al <- strsplit(readr::read_lines('identity_hlaa_longest_cutted.csv', 
                                            n_max = 1), '\",\"')[[1]]

proper_format <- unlist(purrr::map(first_line_al[2:length(first_line_al)], 
                                   .f = ~get_mhcpan_input(.x, 'MHC-I')), 
                        use.names = FALSE)
transf_id <- unlist(purrr::map(purrr::map(gsub('A\\.', 'HLA-A', colnames(id_al)), 
                                          .f = ~strsplit(.x, '\\.')), 
                               .f = ~paste(.x[[1]][1],.x[[1]][2], sep=':')), 
                    use.names = FALSE)

m <- match(transf_id, proper_format)

id_al_sub <- id_al[!is.na(transf_id[m]),!is.na(transf_id[m])]
# write.csv(id_al_sub, 'identityt_hlaa_only_netmhcpan_proper_names_longest_cutted.csv', quote = F)

id_al_sub <- readr::read_delim('identityt_hlaa_only_netmhcpan_proper_names_longest_cutted.csv',
                               delim = ",",col_names = TRUE)
id_al_sub <- data.frame(id_al_sub)
rownames(id_al_sub) <- id_al_sub$X1
id_al_sub <- id_al_sub[,-c(1)]
rownames(id_al_sub) <- colnames(id_al_sub)

id_al_sub$format <- transf_id[!is.na(transf_id[m])]

rownames(id_al_sub) <- gsub('-', '\\.', rownames(id_al_sub))
rownames(id_al_sub) <- gsub(':', '\\.', rownames(id_al_sub))


###
###
###
# Loading EL-scores
###
###
###

allele_specif_cor_matr <- readr::read_delim('allele_specif_cor_matr.csv',
                                            delim = ",",col_names = TRUE)
allele_specif_cor_matr <- data.frame(allele_specif_cor_matr)
rownames(allele_specif_cor_matr) <- allele_specif_cor_matr$X1
allele_specif_cor_matr <- allele_specif_cor_matr[,-c(1)]
allele_specif_cor_matr[is.na(allele_specif_cor_matr)] <- 1

m <- match(colnames(allele_specif_cor_matr), colnames(id_al_sub))
allele_specif_cor_matr <- allele_specif_cor_matr[!is.na(colnames(allele_specif_cor_matr)[m]),
                                                 !is.na(colnames(allele_specif_cor_matr)[m])]
rownames(allele_specif_cor_matr) <- colnames(allele_specif_cor_matr)
rownames(allele_specif_cor_matr) <- gsub('-', '\\.', 
                                         rownames(allele_specif_cor_matr))
rownames(allele_specif_cor_matr) <- gsub(':', '\\.', 
                                         rownames(allele_specif_cor_matr))
allele_specif_cor_matr <- allele_specif_cor_matr[colnames(id_al_sub),
                                                 colnames(id_al_sub)]


res_cor_mean <- cor(c(as.matrix(allele_specif_cor_matr)), 
                    c(as.matrix(id_al_sub)))

allele_specif_cor_matr_m <- melt(replace(allele_specif_cor_matr, 
                                         lower.tri(allele_specif_cor_matr, 
                                                   TRUE), NA), na.rm = TRUE)
v_mean_t_t_agg_m <- melt(replace(id_al_sub, 
                                 lower.tri(id_al_sub, 
                                           TRUE), NA), na.rm = TRUE)

merged_df_mean <- merge(allele_specif_cor_matr_m, v_mean_t_t_agg_m,
                        by='row.names')

# distribution of alignment's scores/netMHCpan scores
ggplot(merged_df_mean, aes(value.y)) +
  geom_density()+theme_classic()

ggplot(merged_df_mean, aes(value.x, value.y)) +
  geom_hex(bins = 100)+theme_classic()

# ggsave('merged_df_mean_sim_al.png', dpi = 300, width = 6, height = 4)

###
###
###
# Loading BA-scores
###
###
###

allele_specif_cor_matr_BA <- readr::read_delim('allele_specif_cor_matr_BA.csv',
                                               delim = ",",col_names = TRUE)
allele_specif_cor_matr_BA <- data.frame(allele_specif_cor_matr_BA)
rownames(allele_specif_cor_matr_BA) <- allele_specif_cor_matr_BA$X1
allele_specif_cor_matr_BA <- allele_specif_cor_matr_BA[,-c(1)]
allele_specif_cor_matr_BA[is.na(allele_specif_cor_matr_BA)] <- 1

m <- match(colnames(allele_specif_cor_matr_BA), colnames(id_al_sub))
allele_specif_cor_matr_BA <- allele_specif_cor_matr_BA[!is.na(colnames(allele_specif_cor_matr_BA)[m]),
                                                       !is.na(colnames(allele_specif_cor_matr_BA)[m])]
rownames(allele_specif_cor_matr_BA) <- colnames(allele_specif_cor_matr_BA)
rownames(allele_specif_cor_matr_BA) <- gsub('-', '\\.', 
                                            rownames(allele_specif_cor_matr_BA))
rownames(allele_specif_cor_matr_BA) <- gsub(':', '\\.', 
                                            rownames(allele_specif_cor_matr_BA))
allele_specif_cor_matr_BA <- allele_specif_cor_matr_BA[colnames(id_al_sub),
                                                       colnames(id_al_sub)]

res_cor_mean <- cor(c(as.matrix(allele_specif_cor_matr_BA)), c(as.matrix(id_al_sub)))

allele_specif_cor_matr_BA_m <- melt(replace(allele_specif_cor_matr_BA, 
                                            lower.tri(allele_specif_cor_matr_BA, 
                                                      TRUE), NA), na.rm = TRUE)
v_mean_t_t_agg_m <- melt(replace(id_al_sub, 
                                 lower.tri(id_al_sub, 
                                           TRUE), NA), na.rm = TRUE)

merged_df_mean <- merge(allele_specif_cor_matr_BA_m, v_mean_t_t_agg_m,
                        by='row.names')


# distribution of alignment's scores/netMHCpan scores
ggplot(merged_df_mean, aes(value.y)) +
  geom_density()+theme_classic()

ggplot(merged_df_mean, aes(value.x, value.y)) +
  geom_hex(bins = 100)+theme_classic()

# ggsave('merged_df_mean_sim_al.png', dpi = 300, width = 6, height = 4)

# find particular example of alleles
merged_df_mean[merged_df_mean$value.x >0.8,][merged_df_mean
                                             [merged_df_mean$value.x >0.8,
                                               ]$value.y >0.35,]

rownames(allele_specif_cor_matr)[grepl('0.821716', 
                  sapply(allele_specif_cor_matr[,'HLA.A25.05'], as.character))]

v_mean_t_t_agg["HLA.A23.46",'HLA.A25.05']
allele_specif_cor_matr['HLA.A25.05',"HLA.A23.46"]
allele_specif_cor_matr_BA['HLA.A25.05',"HLA.A23.46"]


# JTT vs.netMHCpan scores -------------------------------------------------


alleles_comb_res_df <- readr::read_delim('allele_jtt.csv', delim = ",", 
                                         skip=1,col_names = F)
alleles_comb_res_df <- data.frame(alleles_comb_res_df)
rownames(alleles_comb_res_df) <- alleles_comb_res_df$X1
alleles_comb_res_df <- alleles_comb_res_df[,-c(1)]
colnames(alleles_comb_res_df) <- rownames(alleles_comb_res_df)
alleles_comb_res_df[is.na(alleles_comb_res_df)] <- 0

proper_format <- unlist(purrr::map(gsub('\\.', ':', gsub('A\\.', 'A\\*', 
                                    colnames(alleles_comb_res_df))), 
                                   .f = ~get_mhcpan_input(.x, 'MHC-I')), 
                        use.names = FALSE)
transf <- unlist(purrr::map(purrr::map(gsub('A\\.', 'HLA-A', 
                                            colnames(alleles_comb_res_df)), 
                                       .f = ~strsplit(.x, '\\.')), 
                            .f = ~paste(.x[[1]][1],.x[[1]][2], sep=':')), 
                 use.names = FALSE)
m <- match(transf, proper_format)

alleles_comb_res_df <- alleles_comb_res_df[!is.na(transf[m]),!is.na(transf[m])]
# write.csv(as.matrix(alleles_comb_res_df), 'allele_jtt_prop_names.csv', quote = F)


alleles_comb_res_df <- readr::read_delim('allele_jtt_prop_names.csv', 
                                         delim = ",", skip=1,col_names = F)
alleles_comb_res_df <- data.frame(alleles_comb_res_df)
rownames(alleles_comb_res_df) <- alleles_comb_res_df$X1
alleles_comb_res_df <- alleles_comb_res_df[,-c(1)]
colnames(alleles_comb_res_df) <- rownames(alleles_comb_res_df)
alleles_comb_res_df <- apply(alleles_comb_res_df,2,as.numeric)
alleles_comb_res_df <- data.frame(alleles_comb_res_df)
rownames(alleles_comb_res_df) <- colnames(alleles_comb_res_df)

transf <- unlist(purrr::map(purrr::map(gsub('A\\.', 'HLA-A', 
                                            colnames(alleles_comb_res_df)), 
                                       .f = ~strsplit(.x, '\\.')), 
                            .f = ~paste(.x[[1]][1],.x[[1]][2], sep=':')), 
                 use.names = FALSE)
proper_format <- unlist(purrr::map(gsub('\\.', ':', gsub('A\\.', 'A\\*', 
                                      colnames(alleles_comb_res_df))), 
                                   .f = ~get_mhcpan_input(.x, 'MHC-I')), 
                        use.names = FALSE)
m <- match(transf, proper_format)

alleles_comb_res_df <- alleles_comb_res_df[!is.na(transf[m]),!is.na(transf[m])]
alleles_comb_res_df$format <- transf[!is.na(transf[m])]
dim(alleles_comb_res_df)

v_mean <- aggregate(alleles_comb_res_df[-(3452)], list(alleles_comb_res_df[,3452]), 
                    function(x) mean(unique(x)))
rownames(v_mean) <- v_mean$Group.1
v_mean <- v_mean[,-c(1)]
m <- match(transf, proper_format)
v_mean_t <- data.frame(t(v_mean))
v_mean_t$format <- transf[!is.na(transf[m])]

v_mean_t_t_agg <- aggregate(v_mean_t[-(3250)], list(v_mean_t[,3250]), 
                            function(x) mean(unique(x)))
rownames(v_mean_t_t_agg) <- v_mean_t_t_agg$Group.1
v_mean_t_t_agg <- v_mean_t_t_agg[,-c(1)]
rownames(v_mean_t_t_agg) <- gsub('-', '\\.', rownames(v_mean_t_t_agg))
rownames(v_mean_t_t_agg) <- gsub(':', '\\.', rownames(v_mean_t_t_agg))

m <- match(colnames(v_mean_t_t_agg), colnames(allele_specif_cor_matr_BA))
v_mean_t_t_agg <- v_mean_t_t_agg[!is.na(colnames(v_mean_t_t_agg)[m]),
                                 !is.na(colnames(v_mean_t_t_agg)[m])]

###
###
###
# Loading EL-scores
###
###
###

allele_specif_cor_matr <- readr::read_delim('allele_specif_cor_matr.csv',
                                            delim = ",",col_names = TRUE)
allele_specif_cor_matr <- data.frame(allele_specif_cor_matr)
rownames(allele_specif_cor_matr) <- allele_specif_cor_matr$X1
allele_specif_cor_matr <- allele_specif_cor_matr[,-c(1)]
allele_specif_cor_matr[is.na(allele_specif_cor_matr)] <- 1

m <- match(colnames(allele_specif_cor_matr), colnames(id_al_sub))
allele_specif_cor_matr <- allele_specif_cor_matr[!is.na(colnames(allele_specif_cor_matr)[m]),
                                                 !is.na(colnames(allele_specif_cor_matr)[m])]
rownames(allele_specif_cor_matr) <- colnames(allele_specif_cor_matr)
rownames(allele_specif_cor_matr) <- gsub('-', '\\.', 
                                         rownames(allele_specif_cor_matr))
rownames(allele_specif_cor_matr) <- gsub(':', '\\.', 
                                         rownames(allele_specif_cor_matr))
allele_specif_cor_matr <- allele_specif_cor_matr[colnames(id_al_sub),
                                                 colnames(id_al_sub)]


res_cor_mean <- cor(c(as.matrix(allele_specif_cor_matr)), 
                    c(as.matrix(id_al_sub)))

allele_specif_cor_matr_m <- melt(replace(allele_specif_cor_matr, 
                                         lower.tri(allele_specif_cor_matr, 
                                                   TRUE), NA), na.rm = TRUE)
v_mean_t_t_agg_m <- melt(replace(id_al_sub, 
                                 lower.tri(id_al_sub, 
                                           TRUE), NA), na.rm = TRUE)

merged_df_mean <- merge(allele_specif_cor_matr_m, v_mean_t_t_agg_m,
                        by='row.names')

# distribution of alignment's scores/netMHCpan scores
ggplot(merged_df_mean, aes(value.y)) +
  geom_density()+theme_classic()

ggplot(merged_df_mean, aes(value.x, value.y)) +
  geom_hex(bins = 100)+theme_classic()

# ggsave('merged_df_mean_sim_al.png', dpi = 300, width = 6, height = 4)

###
###
###
# Loading BA-scores
###
###
###

allele_specif_cor_matr_BA <- readr::read_delim('allele_specif_cor_matr_BA.csv',
                                               delim = ",",col_names = TRUE)
allele_specif_cor_matr_BA <- data.frame(allele_specif_cor_matr_BA)
rownames(allele_specif_cor_matr_BA) <- allele_specif_cor_matr_BA$X1
allele_specif_cor_matr_BA <- allele_specif_cor_matr_BA[,-c(1)]
allele_specif_cor_matr_BA[is.na(allele_specif_cor_matr_BA)] <- 1

m <- match(colnames(allele_specif_cor_matr_BA), colnames(id_al_sub))
allele_specif_cor_matr_BA <- allele_specif_cor_matr_BA[!is.na(colnames(allele_specif_cor_matr_BA)[m]),
                                                       !is.na(colnames(allele_specif_cor_matr_BA)[m])]
rownames(allele_specif_cor_matr_BA) <- colnames(allele_specif_cor_matr_BA)
rownames(allele_specif_cor_matr_BA) <- gsub('-', '\\.', 
                                            rownames(allele_specif_cor_matr_BA))
rownames(allele_specif_cor_matr_BA) <- gsub(':', '\\.', 
                                            rownames(allele_specif_cor_matr_BA))
allele_specif_cor_matr_BA <- allele_specif_cor_matr_BA[colnames(id_al_sub),
                                                       colnames(id_al_sub)]

res_cor_mean <- cor(c(as.matrix(allele_specif_cor_matr_BA)), 
                    c(as.matrix(id_al_sub)))

allele_specif_cor_matr_BA_m <- melt(replace(allele_specif_cor_matr_BA, 
                                            lower.tri(allele_specif_cor_matr_BA, 
                                                      TRUE), NA), na.rm = TRUE)
v_mean_t_t_agg_m <- melt(replace(id_al_sub, 
                                 lower.tri(id_al_sub, 
                                           TRUE), NA), na.rm = TRUE)

merged_df_mean <- merge(allele_specif_cor_matr_BA_m, v_mean_t_t_agg_m,
                        by='row.names')


# distribution of alignment's scores/netMHCpan scores
ggplot(merged_df_mean, aes(value.y)) +
  geom_density()+theme_classic()

ggplot(merged_df_mean, aes(value.x, value.y)) +
  geom_hex(bins = 100)+theme_classic()

# ggsave('merged_df_mean_sim_al.png', dpi = 300, width = 6, height = 4)


# PAM40 vs.netMHCpan scores -----------------------------------------------


allele_pam40 <- readr::read_delim('allele_pam40.csv', delim = ",", 
                                  skip=1,col_names = F)
allele_pam40 <- data.frame(allele_pam40)
rownames(allele_pam40) <- allele_pam40$X1
allele_pam40 <- allele_pam40[,-c(1)]
colnames(allele_pam40) <- rownames(allele_pam40)

proper_format <- unlist(purrr::map(gsub('\\.', ':', gsub('A\\.', 'A\\*', 
                                                         colnames(allele_pam40))), 
                                   .f = ~get_mhcpan_input(.x, 'MHC-I')), 
                        use.names = FALSE)
transf <- unlist(purrr::map(purrr::map(gsub('A\\.', 'HLA-A', 
                                            colnames(allele_pam40)), 
                                       .f = ~strsplit(.x, '\\.')), 
                            .f = ~paste(.x[[1]][1],.x[[1]][2], sep=':')), 
                 use.names = FALSE)
m <- match(transf, proper_format)

allele_pam40 <- allele_pam40[!is.na(transf[m]),!is.na(transf[m])]
# write.csv(allele_pam40, 'allele_pam40_prop_names.csv', quote = F)

alleles_comb_res_df <- readr::read_delim('allele_pam40_prop_names.csv', 
                                         delim = ",", skip=1,col_names = F)
alleles_comb_res_df <- data.frame(alleles_comb_res_df)
rownames(alleles_comb_res_df) <- alleles_comb_res_df$X1
alleles_comb_res_df <- alleles_comb_res_df[,-c(1)]
colnames(alleles_comb_res_df) <- rownames(alleles_comb_res_df)

alleles_comb_res_df <- data.frame(alleles_comb_res_df)
rownames(alleles_comb_res_df) <- colnames(alleles_comb_res_df)

transf <- unlist(purrr::map(purrr::map(gsub('A\\.', 'HLA-A', 
                                            colnames(alleles_comb_res_df)), 
                                       .f = ~strsplit(.x, '\\.')), 
                            .f = ~paste(.x[[1]][1],.x[[1]][2], sep=':')), 
                 use.names = FALSE)
proper_format <- unlist(purrr::map(gsub('\\.', ':', gsub('A\\.', 'A\\*', 
                                    colnames(alleles_comb_res_df))), 
                                   .f = ~get_mhcpan_input(.x, 'MHC-I')), 
                        use.names = FALSE)
m <- match(transf, proper_format)

alleles_comb_res_df <- alleles_comb_res_df[!is.na(transf[m]),!is.na(transf[m])]
alleles_comb_res_df$format <- transf[!is.na(transf[m])]
dim(alleles_comb_res_df)

alleles_comb_res_df[, -c(3452)] <- sapply(alleles_comb_res_df[, -c(3452)], 
                                          as.numeric)

v_mean <- aggregate(alleles_comb_res_df[-(3452)], list(alleles_comb_res_df[,3452]), 
                    function(x) mean(unique(x)))
rownames(v_mean) <- v_mean$Group.1
v_mean <- v_mean[,-c(1)]
m <- match(transf, proper_format)
v_mean_t <- data.frame(t(v_mean))
v_mean_t$format <- transf[!is.na(transf[m])]

v_mean_t_t_agg <- aggregate(v_mean_t[-(3250)], list(v_mean_t[,3250]), 
                            function(x) mean(unique(x)))
rownames(v_mean_t_t_agg) <- v_mean_t_t_agg$Group.1
v_mean_t_t_agg <- v_mean_t_t_agg[,-c(1)]
rownames(v_mean_t_t_agg) <- gsub('-', '\\.', rownames(v_mean_t_t_agg))
rownames(v_mean_t_t_agg) <- gsub(':', '\\.', rownames(v_mean_t_t_agg))

m <- match(colnames(v_mean_t_t_agg), colnames(allele_specif_cor_matr_BA))
# m <- match(colnames(v_mean_t_t_agg), colnames(allele_specif_cor_matr))
v_mean_t_t_agg <- v_mean_t_t_agg[!is.na(colnames(v_mean_t_t_agg)[m]),
                                 !is.na(colnames(v_mean_t_t_agg)[m])]
v_mean_t_t_agg[is.na(v_mean_t_t_agg)] <- 0

###
###
###
# Loading EL-scores
###
###
###

allele_specif_cor_matr <- readr::read_delim('allele_specif_cor_matr.csv',
                                            delim = ",",col_names = TRUE)
allele_specif_cor_matr <- data.frame(allele_specif_cor_matr)
rownames(allele_specif_cor_matr) <- allele_specif_cor_matr$X1
allele_specif_cor_matr <- allele_specif_cor_matr[,-c(1)]
allele_specif_cor_matr[is.na(allele_specif_cor_matr)] <- 1

m <- match(colnames(allele_specif_cor_matr), colnames(id_al_sub))
allele_specif_cor_matr <- allele_specif_cor_matr[!is.na(colnames(allele_specif_cor_matr)[m]),
                                                 !is.na(colnames(allele_specif_cor_matr)[m])]
rownames(allele_specif_cor_matr) <- colnames(allele_specif_cor_matr)
rownames(allele_specif_cor_matr) <- gsub('-', '\\.', 
                                         rownames(allele_specif_cor_matr))
rownames(allele_specif_cor_matr) <- gsub(':', '\\.', 
                                         rownames(allele_specif_cor_matr))
allele_specif_cor_matr <- allele_specif_cor_matr[colnames(id_al_sub),
                                                 colnames(id_al_sub)]


res_cor_mean <- cor(c(as.matrix(allele_specif_cor_matr)), 
                    c(as.matrix(id_al_sub)))

allele_specif_cor_matr_m <- melt(replace(allele_specif_cor_matr, 
                                         lower.tri(allele_specif_cor_matr, 
                                                   TRUE), NA), na.rm = TRUE)
v_mean_t_t_agg_m <- melt(replace(id_al_sub, 
                                 lower.tri(id_al_sub, 
                                           TRUE), NA), na.rm = TRUE)

merged_df_mean <- merge(allele_specif_cor_matr_m, v_mean_t_t_agg_m,
                        by='row.names')

# distribution of alignment's scores/netMHCpan scores
ggplot(merged_df_mean, aes(value.y)) +
  geom_density()+theme_classic()

ggplot(merged_df_mean, aes(value.x, value.y)) +
  geom_hex(bins = 100)+theme_classic()

# ggsave('merged_df_mean_sim_al.png', dpi = 300, width = 6, height = 4)

###
###
###
# Loading BA-scores
###
###
###

allele_specif_cor_matr_BA <- readr::read_delim('allele_specif_cor_matr_BA.csv',
                                               delim = ",",col_names = TRUE)
allele_specif_cor_matr_BA <- data.frame(allele_specif_cor_matr_BA)
rownames(allele_specif_cor_matr_BA) <- allele_specif_cor_matr_BA$X1
allele_specif_cor_matr_BA <- allele_specif_cor_matr_BA[,-c(1)]
allele_specif_cor_matr_BA[is.na(allele_specif_cor_matr_BA)] <- 1

m <- match(colnames(allele_specif_cor_matr_BA), colnames(id_al_sub))
allele_specif_cor_matr_BA <- allele_specif_cor_matr_BA[!is.na(colnames(allele_specif_cor_matr_BA)[m]),
                                                       !is.na(colnames(allele_specif_cor_matr_BA)[m])]
rownames(allele_specif_cor_matr_BA) <- colnames(allele_specif_cor_matr_BA)
rownames(allele_specif_cor_matr_BA) <- gsub('-', '\\.', 
                                            rownames(allele_specif_cor_matr_BA))
rownames(allele_specif_cor_matr_BA) <- gsub(':', '\\.', 
                                            rownames(allele_specif_cor_matr_BA))
allele_specif_cor_matr_BA <- allele_specif_cor_matr_BA[colnames(id_al_sub),
                                                       colnames(id_al_sub)]

res_cor_mean <- cor(c(as.matrix(allele_specif_cor_matr_BA)), c(as.matrix(id_al_sub)))

allele_specif_cor_matr_BA_m <- melt(replace(allele_specif_cor_matr_BA, 
                                            lower.tri(allele_specif_cor_matr_BA, 
                                                      TRUE), NA), na.rm = TRUE)
v_mean_t_t_agg_m <- melt(replace(id_al_sub, 
                                 lower.tri(id_al_sub, 
                                           TRUE), NA), na.rm = TRUE)

merged_df_mean <- merge(allele_specif_cor_matr_BA_m, v_mean_t_t_agg_m,
                        by='row.names')


# distribution of alignment's scores/netMHCpan scores
ggplot(merged_df_mean, aes(value.y)) +
  geom_density()+theme_classic()

ggplot(merged_df_mean, aes(value.x, value.y)) +
  geom_hex(bins = 100)+theme_classic()

# ggsave('merged_df_mean_sim_al.png', dpi = 300, width = 6, height = 4)



# BLOSUM80 vs.netMHCpan scores --------------------------------------------


allele_blosum80 <- readr::read_delim('allele_blossum80.csv', delim = ",", 
                                     skip=1,col_names = F)
allele_blosum80 <- data.frame(allele_blosum80)
rownames(allele_blosum80) <- allele_blosum80$X1
allele_blosum80 <- allele_blosum80[,-c(1)]
colnames(allele_blosum80) <- rownames(allele_blosum80)

proper_format <- unlist(purrr::map(gsub('\\.', ':', gsub('A\\.', 'A\\*', 
                                                         colnames(allele_blosum80))), 
                                   .f = ~get_mhcpan_input(.x, 'MHC-I')), 
                        use.names = FALSE)
transf <- unlist(purrr::map(purrr::map(gsub('A\\.', 'HLA-A', 
                                            colnames(allele_blosum80)), 
                                       .f = ~strsplit(.x, '\\.')), 
                            .f = ~paste(.x[[1]][1],.x[[1]][2], sep=':')), 
                 use.names = FALSE)
m <- match(transf, proper_format)

allele_blosum80 <- allele_blosum80[!is.na(transf[m]),!is.na(transf[m])]
# write.csv(allele_blosum80, 'allele_blosum80_prop_names.csv', quote = F)

alleles_comb_res_df <- readr::read_delim('allele_blosum80_prop_names.csv', 
                                         delim = ",", skip=1,col_names = F)
alleles_comb_res_df <- data.frame(alleles_comb_res_df)
rownames(alleles_comb_res_df) <- alleles_comb_res_df$X1
alleles_comb_res_df <- alleles_comb_res_df[,-c(1)]
colnames(alleles_comb_res_df) <- rownames(alleles_comb_res_df)

transf <- unlist(purrr::map(purrr::map(gsub('A\\.', 'HLA-A', 
                                            colnames(alleles_comb_res_df)), 
                                       .f = ~strsplit(.x, '\\.')), 
                            .f = ~paste(.x[[1]][1],.x[[1]][2], sep=':')), 
                 use.names = FALSE)

m <- match(transf, proper_format)

alleles_comb_res_df <- alleles_comb_res_df[!is.na(transf[m]),!is.na(transf[m])]
alleles_comb_res_df$format <- transf[!is.na(transf[m])]


v_mean <- aggregate(alleles_comb_res_df[-(3444)], list(alleles_comb_res_df[,3444]), 
                    function(x) mean(unique(x)))
rownames(v_mean) <- v_mean$Group.1
v_mean <- v_mean[,-c(1)]
m <- match(transf, proper_format)
v_mean_t <- data.frame(t(v_mean))
v_mean_t$format <- transf[!is.na(transf[m])]

v_mean_t_t_agg <- aggregate(v_mean_t[-(3242)], list(v_mean_t[,3242]), 
                            function(x) mean(unique(x)))
rownames(v_mean_t_t_agg) <- v_mean_t_t_agg$Group.1
v_mean_t_t_agg <- v_mean_t_t_agg[,-c(1)]
dim(v_mean_t_t_agg)
rownames(v_mean_t_t_agg) <- gsub('-', '\\.', rownames(v_mean_t_t_agg))
rownames(v_mean_t_t_agg) <- gsub(':', '\\.', rownames(v_mean_t_t_agg))

m <- match(colnames(v_mean_t_t_agg), colnames(allele_specif_cor_matr_BA))
v_mean_t_t_agg <- v_mean_t_t_agg[!is.na(colnames(v_mean_t_t_agg)[m]),
                                 !is.na(colnames(v_mean_t_t_agg)[m])]
v_mean_t_t_agg[is.na(v_mean_t_t_agg)] <- 0

###
###
###
# Loading EL-scores
###
###
###

allele_specif_cor_matr <- readr::read_delim('allele_specif_cor_matr.csv',
                                            delim = ",",col_names = TRUE)
allele_specif_cor_matr <- data.frame(allele_specif_cor_matr)
rownames(allele_specif_cor_matr) <- allele_specif_cor_matr$X1
allele_specif_cor_matr <- allele_specif_cor_matr[,-c(1)]
allele_specif_cor_matr[is.na(allele_specif_cor_matr)] <- 1

m <- match(colnames(allele_specif_cor_matr), colnames(id_al_sub))
allele_specif_cor_matr <- allele_specif_cor_matr[!is.na(colnames(allele_specif_cor_matr)[m]),
                                                 !is.na(colnames(allele_specif_cor_matr)[m])]
rownames(allele_specif_cor_matr) <- colnames(allele_specif_cor_matr)
rownames(allele_specif_cor_matr) <- gsub('-', '\\.', 
                                         rownames(allele_specif_cor_matr))
rownames(allele_specif_cor_matr) <- gsub(':', '\\.', 
                                         rownames(allele_specif_cor_matr))
allele_specif_cor_matr <- allele_specif_cor_matr[colnames(id_al_sub),
                                                 colnames(id_al_sub)]


res_cor_mean <- cor(c(as.matrix(allele_specif_cor_matr)), 
                    c(as.matrix(id_al_sub)))

allele_specif_cor_matr_m <- melt(replace(allele_specif_cor_matr, 
                                         lower.tri(allele_specif_cor_matr, 
                                                   TRUE), NA), na.rm = TRUE)
v_mean_t_t_agg_m <- melt(replace(id_al_sub, 
                                 lower.tri(id_al_sub, 
                                           TRUE), NA), na.rm = TRUE)

merged_df_mean <- merge(allele_specif_cor_matr_m, v_mean_t_t_agg_m,by='row.names')

# distribution of alignment's scores/netMHCpan scores
ggplot(merged_df_mean, aes(value.y)) +
  geom_density()+theme_classic()

ggplot(merged_df_mean, aes(value.x, value.y)) +
  geom_hex(bins = 100)+theme_classic()

# ggsave('merged_df_mean_sim_al.png', dpi = 300, width = 6, height = 4)

###
###
###
# Loading BA-scores
###
###
###

allele_specif_cor_matr_BA <- readr::read_delim('allele_specif_cor_matr_BA.csv',delim = ",",col_names = TRUE)
allele_specif_cor_matr_BA <- data.frame(allele_specif_cor_matr_BA)
rownames(allele_specif_cor_matr_BA) <- allele_specif_cor_matr_BA$X1
allele_specif_cor_matr_BA <- allele_specif_cor_matr_BA[,-c(1)]
allele_specif_cor_matr_BA[is.na(allele_specif_cor_matr_BA)] <- 1

m <- match(colnames(allele_specif_cor_matr_BA), colnames(id_al_sub))
allele_specif_cor_matr_BA <- allele_specif_cor_matr_BA[!is.na(colnames(allele_specif_cor_matr_BA)[m]),
                                                       !is.na(colnames(allele_specif_cor_matr_BA)[m])]
rownames(allele_specif_cor_matr_BA) <- colnames(allele_specif_cor_matr_BA)
rownames(allele_specif_cor_matr_BA) <- gsub('-', '\\.', 
                                            rownames(allele_specif_cor_matr_BA))
rownames(allele_specif_cor_matr_BA) <- gsub(':', '\\.', 
                                            rownames(allele_specif_cor_matr_BA))
allele_specif_cor_matr_BA <- allele_specif_cor_matr_BA[colnames(id_al_sub),
                                                       colnames(id_al_sub)]

res_cor_mean <- cor(c(as.matrix(allele_specif_cor_matr_BA)), c(as.matrix(id_al_sub)))

allele_specif_cor_matr_BA_m <- melt(replace(allele_specif_cor_matr_BA, 
                                            lower.tri(allele_specif_cor_matr_BA, 
                                                      TRUE), NA), na.rm = TRUE)
v_mean_t_t_agg_m <- melt(replace(id_al_sub, 
                                 lower.tri(id_al_sub, 
                                           TRUE), NA), na.rm = TRUE)

merged_df_mean <- merge(allele_specif_cor_matr_BA_m, v_mean_t_t_agg_m,by='row.names')


# distribution of alignment's scores/netMHCpan scores
ggplot(merged_df_mean, aes(value.y)) +
  geom_density()+theme_classic()

ggplot(merged_df_mean, aes(value.x, value.y)) +
  geom_hex(bins = 100)+theme_classic()

# ggsave('merged_df_mean_sim_al.png', dpi = 300, width = 6, height = 4)

# BA vs. EL scores --------------------------------------------------------

cor(c(as.matrix(allele_specif_cor_matr)), 
    c(as.matrix(allele_specif_cor_matr_BA)))

allele_specif_cor_matr_BA_m <- melt(replace(allele_specif_cor_matr_BA, 
                                            lower.tri(allele_specif_cor_matr_BA, 
                                                      TRUE), NA), na.rm = TRUE)
allele_specif_cor_matr_m <- melt(replace(allele_specif_cor_matr, 
                                         lower.tri(allele_specif_cor_matr, 
                                                   TRUE), NA), na.rm = TRUE)
merged_df <- merge(allele_specif_cor_matr_BA_m, 
                   allele_specif_cor_matr_m,by='row.names')
ggplot(merged_df, aes(value.x, value.y)) +
  geom_hex(bins = 100)+theme_classic()

ggplot(merged_df, aes(value.x)) +
  geom_density()+theme_classic()

# ggsave('el_ba_cor.png', dpi = 300, width = 6, height = 4)