# sleuth p-value aggregation over GO terms

This is just a repo to save code-work implementing a p-value aggregation of sleuth differential expression p-values of individual transcripts over GO terms. So this is work in progress, but should work on sleuth results in RDS format as saved with `sleuth_save()`.

CAVEAT:
GO terms are not independent, but instead overlap in various ways through their dependencies.
Thus, assumptions of the p-value aggregation are not valid.
A way to work around this might be the information theoretic approach presented in:

Alterovitz, Gil, Michael Xiang, Mamta Mohan, and Marco F. Ramoni. 2007. “GO PaD: The Gene Ontology Partition Database.” Nucleic Acids Research 35 (suppl_1): D322–27. https://doi.org/10.1093/nar/gkl799.

However, `GO PaD` does not seem to exist any more, so the partitioning of the GO database would have to be re-implemented.

