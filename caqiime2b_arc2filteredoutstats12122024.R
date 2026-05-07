# 12/12/2024
# Colleen Ahern

# merging my file of info on what sample IDs were filtered out with Malte's QC report to look at the reads of what samples
# got filtered out after rarefaction/filtering

metadata_a_filteredout <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/ampliseq_fungiplast/metadata_a_12122024_filteredout.txt", header=T, stringsAsFactors=F, sep="\t")
malte_qc <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/ampliseq_fungiplast/output_arc2/overall_summary.tsv", header=T, stringsAsFactors=F, sep="\t")


mafo <- metadata_a_filteredout %>%
  filter(sample.id %in% filtered.out)

mafo <- mafo[,1:11]
mafo <- mafo %>%
  rename('sample' = 'sample.id')

df <- merge(mafo, malte_qc, by = 'sample')
