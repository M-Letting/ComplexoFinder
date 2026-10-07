################################################################################
## Handler of lookup tables                                                   ##
## If required lookup tables are not found, they will be downloaded from      ##
## the specified online sources                                               ##
################################################################################

# Required lookup tables and their sources
# all_human_sequeces.rds / human_proteome_with_isoforms.fasta
#     Source: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
#     Source: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz
# ebi_cp_lookup.rds
#     Source: https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv
# HUMAN_9606_idmapping.dat
#     Source: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz
# human_proteome.rds
#     Source: https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=tsv&query=(organism_id:9606)%20AND%20(reviewed:true)&fields=accession,gene_names

# Libraries
library(data.table)

# Set working directory to 00-LookupTables
setwd(here::here())
setwd("ComplexoFinder/00-LookupTables")

# Check which lookup tables are available
lookup_tables <- data.table(
  Name = c(
    "all_human_sequences",
    "ebi_cp_lookup",
    "uniprot_idmapping",
    "human_proteome"
  ),
  File = c(
    "all_human_sequences.rds",
    "ebi_cp_lookup.rds",
    "HUMAN_9606_idmapping.dat",
    "human_proteome.rds"
  ),
  Available = c(
    file.exists("all_human_sequences.rds"),
    file.exists("ebi_cp_lookup.rds"),
    file.exists("HUMAN_9606_idmapping.dat"),
    file.exists("human_proteome.rds")
  )
)

# Status message
if (all(lookup_tables$Available)) {
  # Load lookup tables into environment
  cat(paste0(rep("#", 80), collapse = ""), "\n")
  cat("All lookup tables are available. No download needed.\n")
  cat("Loading accession to gene name mapping...\n")

  # Load UniProt mapping file
  uniprot_map <- fread(
    "HUMAN_9606_idmapping.dat",
    sep = "\t",
    header = FALSE
  )
  names(uniprot_map) <- c("UniProtKB_AC", "ID_type", "ID_value")
  uniprot_map <- uniprot_map[ID_type == "Gene_Name"]
  gene_lookup <- setNames(uniprot_map$ID_value, uniprot_map$UniProtKB_AC)

  cat("  Loaded", length(gene_lookup), "UniProt to gene name mappings\n")

  # Load Protein Sequence Library for Peptide Position Finding
  sequence_library_file <- "all_human_sequences.rds"
  sequence_library <- readRDS(sequence_library_file)

  # Create lookup for fast access
  sequence_lookup <- setNames(
    sequence_library$sequence,
    sequence_library$accession
  )

  cat("  Loaded", nrow(sequence_library), "sequences from RDS\n")

  # Load EBI Complex Portal lookup
  EBI_LT <- readRDS("ebi_cp_lookup.rds")

  cat("  Loaded EBI Complex Portal lookup with", nrow(EBI_LT), "rows\n")

  # Load human proteome lookup
  human_proteome <- readRDS("human_proteome.rds")

  cat("  Loaded human proteome lookup with", nrow(human_proteome), "rows\n")

  # Set working directory back to project root
  setwd(here::here())

  # Stop execution since all tables are loaded
  stop("All lookup tables loaded successfully!")
} else {
  cat("Some lookup tables are missing. Now downloading and processing...\n")
  cat("This might take some time...\n")
}

####################################
## Download missing lookup tables ##
####################################

# Check for all human sequences lookup and download if not available
if (!lookup_tables[Name == "all_human_sequences", Available]) {
  cat("Downloading all human sequences from UniProt...\n")

  # Get zip files from UniProt and read fasta
  con <- gzcon(url(
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
  ))
  lines <- readLines(con)
  close(con)
  con <- gzcon(url(
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz"
  ))
  lines <- c(lines, readLines(con))
  close(con)

  # Parse FASTA
  header_idx <- which(startsWith(lines, ">"))
  if (length(header_idx) == 0L) {
    stop("No FASTA headers found in UniProt download.")
  }

  accessions <- sub("^>[^|]*\\|([^|]+)\\|.*", "\\1", lines[header_idx])
  start_idx <- header_idx + 1L
  end_idx <- c(header_idx[-1L] - 1L, length(lines))

  n_entries <- length(header_idx)
  sequences <- character(n_entries)
  pb <- txtProgressBar(min = 0, max = n_entries, style = 3)

  for (i in seq_len(n_entries)) {
    if (start_idx[i] <= end_idx[i]) {
      sequences[i] <- paste0(lines[start_idx[i]:end_idx[i]], collapse = "")
    } else {
      sequences[i] <- ""
    }
    if (i %% 250L == 0L || i == n_entries) {
      setTxtProgressBar(pb, i)
    }
  }

  close(pb)

  # Save sequence library as data.table with columns "accession" and "sequence"
  sequence_library <- data.table(accession = accessions, sequence = sequences)
  cat("  Parsed", nrow(sequence_library), "protein sequences from FASTA\n")

  saveRDS(sequence_library, "all_human_sequences.rds")
  cat("  Saved sequence library to all_human_sequences.rds\n")
}

# Check for EBI Complex Portal lookupt and download if not available
if (!lookup_tables[Name == "ebi_cp_lookup", Available]) {
  cat("Downloading EBI Complex Portal lookup table...\n")

  # Read EBI CP file
  ebi_cp_data <- data.table::fread(
    "https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv",
    sep = "\t",
    header = TRUE
  )

  # Remove non-relevant columns
  ebi_cp_data <- ebi_cp_data[, .(
    complex_id = `#Complex ac`,
    complex_name = `Recommended name`,
    uniprot_id = `Identifiers (and stoichiometry) of molecules in complex`
  )]

  # Convert to long format by splitting uniprot_id by |
  ebi_cp_long <- ebi_cp_data[,
    .(uniprot_id = unlist(strsplit(uniprot_id, "\\|"))),
    by = .(complex_id, complex_name)
  ]

  # Add stoichiometry column by extracting number in parentheses from uniprot_id
  ebi_cp_long$Stoichiometry <- as.integer(gsub(
    ".*\\((\\d+)\\)",
    "\\1",
    ebi_cp_long$uniprot_id
  ))

  # Remove stoichiometry from uniprot_id
  ebi_cp_long$uniprot_id <- gsub("\\s*\\(\\d+\\)", "", ebi_cp_long$uniprot_id)

  # Remove all rows with CHEBI in uniprot_id
  ebi_cp_long <- ebi_cp_long[!grepl("CHEBI", uniprot_id)]

  # If [ is in uniprot_id, remove [ and ]
  ebi_cp_long$uniprot_id <- gsub("\\[|\\]", "", ebi_cp_long$uniprot_id)

  # If , in uniprot_id, split into multiple rows by ,
  ebi_cp_long <- ebi_cp_long[,
    .(uniprot_id = unlist(strsplit(uniprot_id, ","))),
    by = .(complex_id, complex_name, Stoichiometry)
  ]

  # Trim whitespace from uniprot_id
  ebi_cp_long$uniprot_id <- trimws(ebi_cp_long$uniprot_id)

  # Reorder columns
  ebi_cp_long <- ebi_cp_long[, .(
    complex_id,
    complex_name,
    uniprot_id,
    Stoichiometry
  )]

  # Map UniProt IDs to gene names using the UniProt ID mapping file
  if (!file.exists("HUMAN_9606_idmapping.dat")) {
    cat("  Downloading UniProt ID mapping for gene name lookup...\n")
    idmapping_dl <- data.table::fread(
      "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz",
      sep = "\t", header = FALSE
    )
    data.table::fwrite(idmapping_dl, "HUMAN_9606_idmapping.dat", sep = "\t", col.names = FALSE)
    rm(idmapping_dl)
    cat("  Saved UniProt ID mapping file\n")
  }

  idmapping_raw <- data.table::fread(
    "HUMAN_9606_idmapping.dat", sep = "\t", header = FALSE
  )
  names(idmapping_raw) <- c("UniProtKB_AC", "ID_type", "ID_value")
  gene_map <- idmapping_raw[ID_type == "Gene_Name"][
    , .(gene_name = ID_value[1L]), by = .(base_uniprot = UniProtKB_AC)
  ]
  rm(idmapping_raw)

  # Strip isoform suffix for lookup (O60706-1 -> O60706) but preserve original accession
  ebi_cp_long[, base_uniprot := sub("-\\d+$", "", uniprot_id)]
  ebi_cp_long <- gene_map[ebi_cp_long, on = "base_uniprot"]
  ebi_cp_long[, base_uniprot := NULL]

  # save ebi_cp_long as rds
  saveRDS(ebi_cp_long, "ebi_cp_lookup.rds")
}

# Check for UniProt ID mapping file and download if not available
if (!lookup_tables[Name == "uniprot_idmapping", Available]) {
  cat("Downloading UniProt ID mapping file...\n")
  idmapping <- data.table::fread(
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz",
    sep = "\t"
  )

  # Save as is without processing
  fwrite(idmapping, "HUMAN_9606_idmapping.dat", sep = "\t")
}

# Check for human proteome lookup and download if not available
if (!lookup_tables[Name == "human_proteome", Available]) {
  cat("Downloading human proteome from UniProt...\n")

  url <- "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=tsv&query=(organism_id:9606)%20AND%20(reviewed:true)&fields=accession,gene_names"

  # Download and read TSV (faster than FASTA parsing)
  temp <- tempfile()
  download.file(url, temp)
  human_proteome_df <- read.delim(gzfile(temp), stringsAsFactors = FALSE)
  unlink(temp)

  # Rename columns
  names(human_proteome_df) <- c("Accession", "Gene_name")

  # Extract first gene name if multiple are listed
  human_proteome_df$Gene_name <- sub(";.*", "", human_proteome_df$Gene_name)

  human_proteome <- human_proteome_df

  # Move everything after first space in Gene_name column to a new column
  human_proteome$Gene_name_full <- human_proteome$Gene_name
  human_proteome$Gene_name <- sub(" .*", "", human_proteome$Gene_name)

  human_proteome$Gene_name <- trimws(human_proteome$Gene_name)

  saveRDS(human_proteome, file = "human_proteome.rds")
}

# Set working directory back to project root
setwd(here::here())

# Run this script to load all lookup tables into environment
source("quantification_complex/00-LookupTables/01-GetLookupTables.R")
