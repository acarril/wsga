# Generates the rddsga_synth dataset and saves it to data/rddsga_synth.rda
# Run from the package root: Rscript inst/generate_synth.R

set.seed(42)

n <- 2000

# Running variable centered at 0; cutoff at 0
x <- runif(n, -1, 1)

# Binary subgroup indicator (roughly balanced)
sgroup <- rbinom(n, 1, 0.5)

# Continuous moderator: imbalanced across subgroups at the cutoff
# (sgroup=1 has higher mean m near the cutoff)
m <- rnorm(n, mean = 0.5 * sgroup, sd = 1)

# Sharp treatment assignment
Z <- as.integer(x > 0)

# True RD effects: 2 for G=0, 4 for G=1 (true diff = 2)
y <- 2 * Z * (sgroup == 0) + 4 * Z * (sgroup == 1) + 0.5 * m +
     0.3 * x + rnorm(n, sd = 0.5)

# Fuzzy treatment: compliance ~80% for both groups
u_comply <- runif(n)
D <- integer(n)
D[Z == 1] <- as.integer(u_comply[Z == 1] < 0.80)   # 80% take-up above cutoff
D[Z == 0] <- as.integer(u_comply[Z == 0] < 0.10)   # 10% always-takers

rddsga_synth <- data.frame(
  x      = x,
  sgroup = sgroup,
  m      = m,
  Z      = Z,
  D      = D,
  y      = y
)

save(rddsga_synth, file = "data/rddsga_synth.rda", compress = "bzip2")
message("Saved data/rddsga_synth.rda  (n = ", nrow(rddsga_synth), ")")
