# b_nk_example.R

# Load the ALRA package
# Assuming it's installed, otherwise: remotes::install_github("hk20013106/ALRA")
library(ALRA)

# Load the example data included with the package
data("b_nk_example")
data("labels_example")

# Normalize the data
A_norm <- normalize_data(b_nk_example)

# Choose the rank k
k_choice <- choose_k(A_norm)

# Run ALRA
A_norm_completed <- alra(A_norm, k = k_choice$k)[[3]]

# Print statistics before ALRA
print("Before ALRA:")
print(aggregate(A_norm[, c("NCAM1", "CR2")], by = list(" " = labels_example), FUN = function(x) round(c(percent = 100 * sum(x > 0) / length(x)), 1)))

# Print statistics after ALRA
print("After ALRA:")
print(aggregate(A_norm_completed[, c("NCAM1", "CR2")], by = list(" " = labels_example), FUN = function(x) round(c(percent = 100 * sum(x > 0) / length(x)), 1)))
