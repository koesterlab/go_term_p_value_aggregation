# Analyses adapted along the lines of:
# * https://github.com/pachterlab/aggregationDE/blob/master/SRP100701/R/tx_pipeline.R
# * https://github.com/pachterlab/aggregationDE/blob/master/SRP100701/R/GO_analysis.R
#
# As published in:
#   Yi, Lynn, Harold Pimentel, Nicolas L. Bray, and Lior Pachter. 2018.
#   “Gene-Level Differential Analysis at Transcript-Level Resolution.”
#   Genome Biology 19 (1): 53. https://doi.org/10.1186/s13059-018-1419-z.
#

# gene-centric annotation of the Mus musculus genome
library(org.Mm.eg.db) # bioconductor-org.Mm.eg.db
# GO Term explanations
library(GO.db) # bioconductor-go.db
library(stats)
# provide sane and readable data manipulation syntax and lots of useful features
library(tidyverse) # r-tidyverse
# provide the Lancaster aggregation function
library(aggregation) # r-aggregation


# load data
print('loading data')
tx_results <- readRDS("../sleuth/diffexp/modelXYZ.rds") %>%
  filter( !is.na(ens_gene))

# create a map from GO IDs to the GO terms explaining them
go_map <- mapIds(GO.db,
                keys=keys(GO.db, keytype="GOID"),
                column="TERM",
                keytype="GOID") %>%
          enframe(name = "go_id", value = "go_term")

# create a mapping from ENSEMBL Transcript IDs to GO terms
transcripts_to_go_terms <- AnnotationDbi::select(org.Mm.eg.db,
                                  keys = keys(org.Mm.eg.db, keytype = "ENSEMBLTRANS"),
                                  columns = c("GO", "ONTOLOGY"),
                                  keytype = "ENSEMBLTRANS"
                                  ) %>%
    rename( target_id = ENSEMBLTRANS, go_id = GO ) %>%
    # only include GO Terms from the BP (biological process) ontology set
    # alternatives are: CC (cellular component) and MF (molecular function)
    # also see: http://geneontology.org/docs/ontology-documentation/
    # further filtering could e.g. leverage the EVIDENCE annotation codes as 
    # described at: http://geneontology.org/docs/guide-go-evidence-codes/
    filter( ONTOLOGY == "BP")

# do direct GO term p-value aggregation over transcripts
tx_results_go <- full_join(tx_results %>%
                             filter(test_stat >= 0) %>%
                             select(target_id, pval, mean_obs) %>%
                             mutate(target_id = str_replace(target_id, "\\..+$", "")),
                           transcripts_to_go_terms,
                           by = "target_id") %>%
    drop_na() %>%
    group_by(go_id) %>%
    # Lancaster p-value aggregation, weighted by the mean_obs of each transcript
    summarise( lan = lancaster(pval, mean_obs)) %>%
    # as the GO term aggregation does not incorporate the hierarchical and overlapping
    # structure of GO terms, p-values are not independent and only a multiple
    # testing control that does not assume this can be used
    mutate( lan_bonferroni = p.adjust(lan, method = 'bonferroni') ) %>%
    # add GO term explanations onto the results
    left_join(go_map, by="go_id")


#fisher exact test to test enrichment of GO terms in GO analysis
calculate_enrichment <- function(results, term)
{
  terms <- results %>%
    filter( grepl(term, go_term) )
  print(terms$go_term)
  n_total_all      <- nrow(results)
  n_selected_all   <- results %>% summarise( sum(lan_bonferroni < .05, na.rm=TRUE) )
  n_total_terms    <- nrow(terms)
  n_selected_terms <- terms %>% summarise( sum(lan_bonferroni <.05, na.rm=TRUE) )
  
  n_selected_all_not_terms = n_selected_all - n_selected_terms
  n_terms_not_selected     = n_total_terms  - n_selected_terms
  remaining                = n_total_all    - n_selected_terms - n_selected_all_not_terms - n_terms_not_selected
  x <- matrix( c( as.numeric(n_selected_terms),     as.numeric(n_selected_all_not_terms),
                  as.numeric(n_terms_not_selected), as.numeric(remaining)),
                  2)
  print(x)
  p <- fisher.test(x, alternative = "greater")
  p
}

print('Calculate enrichment for all GO terms containing "immun".')
enrichment <- calculate_enrichment(tx_results_go, 'immun')
