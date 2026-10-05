#!/usr/bin/env Rscript
# Evaluate whether source-specific marker discordances change complete GIFT calls.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: Rscript benchmarks/gift_impact.R RESULTS_DIR GIFTER_SQLITE", call. = FALSE)
}
results_dir <- args[[1]]
db_path <- args[[2]]
connection <- DBI::dbConnect(RSQLite::SQLite(), db_path)

read_calls <- function(path) {
  cols <- c("genome_id", "gene_id", "namespace", "accession")
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(cols %in% names(x))) stop("missing marker columns in ", path)
  unique(x[, cols, drop = FALSE])
}

baseline <- read_calls(file.path(results_dir, "giftag-combined.tsv"))
panel <- sort(unique(baseline$genome_id))
if (length(panel) == 0L) stop("empty giftag marker table")

extract_gifts <- function(calls, scenario) {
  result <- gifter::evaluate_gifts_community(
    calls, db = connection, workers = 1L, progress = FALSE
  )
  rows <- lapply(names(result$results), function(genome) {
    gifts <- result$results[[genome]]$gifts
    data.frame(
      scenario = scenario, genome_id = genome, gift_id = gifts$gift_id,
      gift_type = gifts$gift_type, complete = gifts$complete,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

scenarios <- read.delim(file.path(results_dir, "agreement.tsv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
scenarios <- unique(scenarios[, c("tool", "scope", "source"), drop = FALSE])

select_source <- function(calls, source) {
  switch(source,
    kofam = calls$namespace == "KO",
    dbcan_family = calls$namespace == "CAZY" & !grepl("_e", calls$accession),
    dbcan_sub = calls$namespace == "CAZY" & grepl("_e", calls$accession),
    pfam = calls$namespace == "PFAM",
    ncbifam = calls$namespace %in% c("NCBIFAM", "TIGRFAM")
  )
}

all_gifts <- list(extract_gifts(baseline, "giftag"))
for (i in seq_len(nrow(scenarios))) {
  entry <- scenarios[i, ]
  label <- paste(entry$tool, entry$scope, entry$source, sep = "-")
  reference <- read_calls(file.path(results_dir, paste0(label, "-combined.tsv")))
  hybrid <- unique(rbind(baseline[!select_source(baseline, entry$source), , drop = FALSE],
                         reference))
  if (!identical(sort(unique(hybrid$genome_id)), panel)) {
    stop("hybrid marker table lost a genome: ", label)
  }
  all_gifts[[length(all_gifts) + 1L]] <- extract_gifts(hybrid, label)
}

stacks <- list(
  "stack-targeted" = c(
    "KofamScan-subset-kofam", "run_dbcan-subset-dbcan_family",
    "run_dbcan-subset-dbcan_sub", "HMMER3-subset-pfam", "HMMER3-subset-ncbifam"
  ),
  "stack-full-interpro5" = c(
    "KofamScan-full-kofam", "run_dbcan-full-dbcan_family",
    "run_dbcan-full-dbcan_sub", "InterProScan5-full-pfam",
    "InterProScan5-full-ncbifam"
  ),
  "stack-full-interpro6" = c(
    "KofamScan-full-kofam", "run_dbcan-full-dbcan_family",
    "run_dbcan-full-dbcan_sub", "InterProScan6-full-pfam",
    "InterProScan6-full-ncbifam"
  )
)
all_sources <- c("kofam", "dbcan_family", "dbcan_sub", "pfam", "ncbifam")
covered <- Reduce(`|`, lapply(all_sources, function(source) select_source(baseline, source)))
for (label in names(stacks)) {
  reference <- do.call(rbind, lapply(stacks[[label]], function(component) {
    read_calls(file.path(results_dir, paste0(component, "-combined.tsv")))
  }))
  hybrid <- unique(rbind(baseline[!covered, , drop = FALSE], reference))
  if (!identical(sort(unique(hybrid$genome_id)), panel)) {
    stop("stack marker table lost a genome: ", label)
  }
  all_gifts[[length(all_gifts) + 1L]] <- extract_gifts(hybrid, label)
}

gifts <- do.call(rbind, all_gifts)
write.table(gifts, file.path(results_dir, "gift-calls.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
baseline_gifts <- gifts[gifts$scenario == "giftag", ]
key <- paste(baseline_gifts$genome_id, baseline_gifts$gift_id, sep = "\t")
other <- gifts[gifts$scenario != "giftag", ]
baseline_index <- match(paste(other$genome_id, other$gift_id, sep = "\t"), key)
if (anyNA(baseline_index)) stop("GIFT catalogue changed between evaluations")
other$giftag_complete <- baseline_gifts$complete[baseline_index]
discordant <- other[other$complete != other$giftag_complete, ]
write.table(discordant, file.path(results_dir, "gift-discordant.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

summary <- do.call(rbind, lapply(unique(other$scenario), function(label) {
  scenario_calls <- other[other$scenario == label, , drop = FALSE]
  do.call(rbind, lapply(unique(scenario_calls$gift_type), function(type) {
    x <- scenario_calls[scenario_calls$gift_type == type, , drop = FALSE]
    data.frame(
      scenario = label, gift_type = type, genome_gift_calls = nrow(x),
      giftag_complete = sum(x$giftag_complete),
      scenario_complete = sum(x$complete),
      gained = sum(x$complete & !x$giftag_complete),
      lost = sum(!x$complete & x$giftag_complete),
      flipped_genome_gift_calls = sum(x$complete != x$giftag_complete)
    )
  }))
}))
write.table(summary, file.path(results_dir, "gift-impact-summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(results_dir, "gifter-session-info.txt"))
DBI::dbDisconnect(connection)
cat("GIFT call flips:", nrow(discordant), "\n")
