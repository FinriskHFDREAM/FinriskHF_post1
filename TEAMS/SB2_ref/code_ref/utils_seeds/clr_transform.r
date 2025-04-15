clr_transform <- function(physeq, pseudocount = 1e-6) {
  # Extract and transpose OTU table to samples x taxa
  otu <- otu_table(physeq)
  if (taxa_are_rows(physeq)) {
    otu <- t(otu)
  } else {
    otu <- as(otu, "matrix")
  }

  # Add pseudocount to avoid log(0)
  otu <- otu + pseudocount

  # Calculate CLR: log(value / geometric mean) per sample
  geometric_mean <- function(x) exp(mean(log(x)))
  otu_clr <- t(apply(otu, 1, function(x) log(x / geometric_mean(x))))

  # Replace OTU table in phyloseq object
  otu_table(physeq) <- otu_table(otu_clr, taxa_are_rows = FALSE)

  return(physeq)
}
