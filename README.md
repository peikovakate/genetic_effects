# genetic_effects

Workflow:

1. Build connected components from fine-mapped credible sets
2. Fetch for variants from connected components in eQTL sumstat files for each qtl group
3. Combine eqtl susmstats across all qtl groups (merge tsvs into one)
4. Select eqtls that are most frequent (or at least in 95% of qtlGroups)
5. Extract effects (beta, se, p-value) into sep tibble
6. Find sharing and similarities