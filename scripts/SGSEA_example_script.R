# =========SGSEA Example Script============#
# This script demonstrates the step-by-step workflow of SGSEA analysis.

# Load SGSEA package
library(SGSEA)

# Step 1: Load Example Data
data("KIRC", package = "SGSEA")
head(KIRC)

# Step 2: Extract survival data and gene expression data
survTime <- as.numeric(KIRC[, 2])    # 2nd column = survival time
survStatus <- as.numeric(KIRC[, 3])  # 3rd column = survival status
gene_expr <- KIRC[, 4:ncol(KIRC)]    # Columns 4+ are gene expression

# Step 3: Normalize the Gene Expression Data (EXCLUDING Survival Columns)
# User should filter the genes based on their own choice before this step
normalized_gene_expr <- getNorm(gene_expr, plot = FALSE)

# Step 4: Compute Log Hazard Ratio (LHR)
# Ensure `survTime` and `survStatus` are Valid
cat("Checking survival variables:\n")
print(summary(survTime))  # Ensure times are valid
print(table(survStatus, useNA = "always"))  # Ensure status contains only 0 or 1

# Convert `survStatus` from (1 = alive, 2 = dead) to (0 = alive, 1 = dead)
survStatus[survStatus == 1] <- 0
survStatus[survStatus == 2] <- 1

# Ensure survStatus contains only 0 and 1
valid_status <- survStatus %in% c(0, 1)
survTime <- survTime[valid_status]
survStatus <- survStatus[valid_status]
normalized_gene_expr <- normalized_gene_expr[valid_status, ]  # Keep only valid rows


lhr_results <- getLHR(normalizedData = normalized_gene_expr,
                      survTime = survTime,
                      survStatus = survStatus)


# Step 5: Retrieve Pathways
# Option 1: Get Reactome Pathways (for general pathway analysis)
rpathways <- getReactome(species = "human")

# Option 2: Get GO Terms (for Gene Ontology analysis)
# gopathways <- getGO(species = "human")

# Step 6: Run SGSEA Analysis
sgsea_results <- getSGSEA(pathways = rpathways,
                          stats = lhr_results,
                          minGenes = 5,
                          maxGenes = 500)

# Step 7: Generate an Enrichment Plot for a Specific Pathway
# Selecting the first pathway from the SGSEA results table
# or typing the name directly: "Homo sapiens: Cell Cycle, Mitotic"
selected_pathway <- sgsea_results$rpathway[1]

# Calling getEnrichPlot() with the correct argument structure
getEnrichPlot(pathways = rpathways,
              pathwayName = selected_pathway,
              stats = lhr_results)

# Step 8: Get the Top 20 Enriched Pathways (10 top and 10 bottom)
top_pathways <- getTop10(sgsea_results, rpathways, lhr_results, plotParam = 0.15)

# Save Results to CSV
write.csv(sgsea_results, "SGSEA_Results.csv", row.names = FALSE)


# ===========  GSEA Example Script=============#
# Step 1 : load the data
data("KIRC_DEA", package = "SGSEA")
head(KIRC_DEA)

# Step 2: compute log fold change
lfc_results <-getLFC(KIRC_DEA,cancer,normal,70,70)

# Step 3:  Retrieve Pathways
# Option 1: Get Reactome Pathways (for general pathway analysis)
rpathways <- getReactome(species = "human")

# Option 2: Get GO Terms (for Gene Ontology analysis)
# gopathways <- getGO(species = "human")

# Step 4: Run GSEA Analysis
dea_results <- getSGSEA(pathways = rpathways,
                     stats = lfc_results,
                     minGenes = 5,
                     maxGenes = 500)

# Step 5: Get the Top 20 Enriched Pathways (10 top and 10 bottom)
getTop10(dea_results, rpathways,lfc_results, 0.15)

# Step 6: Generate an Enrichment Plot for a Specific Pathway
getEnrichPlot(rpathways,"Homo sapiens: Cell Cycle, Mitotic",lfc_results)
