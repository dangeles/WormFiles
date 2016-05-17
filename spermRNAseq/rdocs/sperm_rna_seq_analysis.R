#source("http://bioconductor.org/biocLite.R")
#biocLite("rhdf5")
#biocLite("biomaRt")

setwd('~/WormFiles/spermRNAseq/rdocs')
#gene info for sleuth
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "celegans_gene_ensembl")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)


#install.packages("devtools")

#say no to binary files
#devtools::install_github("pachterlab/sleuth")

library("sleuth")

#point to your directory+
base_dir <- "~/WormFiles/spermRNAseq/rdocs/sleuth"

#get ids
sample_id <- dir(file.path(base_dir, "results"))
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "results", id, "kallisto"))
kal_dirs


s2c <- read.table(file.path(base_dir, "sperm_rnaseq_info.txt"), header = TRUE, stringsAsFactors= FALSE)
s2c <- dplyr::select(s2c, sample = experiment, genotype, age)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)


#prepend and make object, state maximum model here
so <- sleuth_prep(s2c, ~ genotype+age+sperm, target_mapping= t2g)

#fit the models
so <- sleuth_fit(so,~ genotype + age, fit_name = 'minimal')
so <- sleuth_fit(so,~ genotype + age+sperm, fit_name = 'full')

#Wald test implementations
#full
so <- sleuth_wt(so, which_beta = 'ageold', which_model = 'full')
so <- sleuth_wt(so, which_beta = 'genotypemt', which_model = 'full')
so <- sleuth_wt(so, which_beta = 'spermhollow', which_model = 'full')

#minimal
so <- sleuth_wt(so, which_beta = 'ageold', which_model = 'minimal')
so <- sleuth_wt(so, which_beta = 'genotypemt', which_model = 'minimal')


#likelihood test
so <- sleuth_lrt(so, 'minimal', 'full')

#if you want to look at shiny
sleuth_live(so)


#following line is to publish...
#saveRDS(so, file = '~/shiny/sleuth/aging_fog2_AngelesAndLeighton_2016.rds')

#write results to tables
results_table <- sleuth_results(so, 'ageold','full', test_type= 'wt')
write.csv(results_table, "~/WormFiles/spermRNAseq/input/agebeta_wt.csv")

results_table <- sleuth_results(so, 'genotypemt','full', test_type= 'wt')
write.csv(results_table, "~/WormFiles/spermRNAseq/input/genotypebeta_wt.csv")

results_table <- sleuth_results(so, 'spermhollow','full', test_type= 'wt')
write.csv(results_table, "~/WormFiles/spermRNAseq/input/spermhollowbeta_wt.csv")

results_table <- sleuth_results(so, 'minimal:full', test_type='lrt')
write.csv(results_table, "~/WormFiles/spermRNAseq/input/lrt.csv")


