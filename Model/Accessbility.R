# Step 1: Define the virus pairs to be ignored (this can also be expanded to multiple groups).
excluded_pairs <- list(c("COVID-19", "IBV"))

# helper function
is_excluded <- function(f, t) {
  any(sapply(excluded_pairs, function(p) p[1] == f && p[2] == t))
}

# Step 2: Construct the access_matrix, skipping excluded pairs.
access_matrix <- matrix(NA, nrow = length(virus_nodes), ncol = length(virus_nodes),
                        dimnames = list(virus_nodes, virus_nodes))

for (i in seq_len(nrow(result_table))) {
  from <- result_table$from_virus[i]
  to <- result_table$to_virus[i]
  
  if (!any(sapply(excluded_pairs, function(p) p[1] == from && p[2] == to))) {
    access_matrix[from, to] <- 1 - result_table$final_score[i]
  }
}

# Step 3: Construct the path propagation matrix (excluding excluded pairs in direct and indirect paths)
access_total <- matrix(0, nrow = length(virus_nodes), ncol = length(virus_nodes),
                       dimnames = list(virus_nodes, virus_nodes))

gamma <- 0.25

# Build access_total, skipping all paths (direct or indirect) involving SARS2 → IBV.
for (from in virus_nodes) {
  for (to in virus_nodes) {
    if (from != to) {
      direct <- if (!is_excluded(from, to) && !is.na(access_matrix[from, to])) access_matrix[from, to] else 0
      indirect <- 0
      
      for (mid in virus_nodes) {
        if (mid != from && mid != to &&
            !is_excluded(from, mid) &&
            !is_excluded(mid, to)) {
          
          via1 <- access_matrix[from, mid]
          via2 <- access_matrix[mid, to]
          
          if (!is.na(via1) && !is.na(via2)) {
            indirect <- indirect + gamma * via1 * via2
          }
        }
      }
      
      access_total[from, to] <- min(1, direct + indirect)
    }
  }
}


competition_strength <- 1 - access_total

df_comp <- melt(competition_strength, varnames = c("from", "to"), value.name = "competition")
