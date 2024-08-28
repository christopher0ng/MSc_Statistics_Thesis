stan_code <- "
data {
  // Data for the multinomial model
  int<lower=1> Y; // Number of years
  int<lower=1> V; // Number of variant combinations
  int<lower=0> Z[V, Y]; // Number of sequenced dengue cases for each variant in each year
  
  // Data for the spatial Poisson model
  int<lower=0> N_spatial; // Number of spatial units
  int<lower=0> N_edges; // Number of edges for spatial structure
  array[N_edges] int<lower=1, upper=N_spatial> node1; // node1[i] adjacent to node2[i]
  array[N_edges] int<lower=1, upper=N_spatial> node2; // and node1[i] < node2[i]
  
  int W[N_spatial,Y]; // total number of dengue cases per spatial unit per year
  int<lower=1> K; // number of covariates
  matrix[N_spatial, K] x; // design matrix for covariates
  
  real<lower=0> scaling_factor; // scales the variance of the spatial effects
}

parameters {
  // Parameters for the multinomial model
  real gamma[V]; // Intercept 
  real phi;  // AR1 Autoregressive coefficient
  real<lower=0> sigma_alpha; // AR1 standard deviation
  matrix[V, Y] alpha;

  // Parameters for the spatial Poisson model
  real beta0; // intercept
  vector[K] betas; // covariates
  real<lower=0> sigma; // overall standard deviation
  real<lower=0, upper=1> rho; // proportion unstructured vs. spatially structured variance
  vector[N_spatial] theta; // heterogeneous effects
  vector[N_spatial] phi_spatial; // spatial effects
}

transformed parameters {
  // Transformed parameters for the multinomial model
  matrix[V, Y] p; // Probabilities in Multinomial distribution
  for (y in 1:Y) {
    p[, y] = softmax(alpha[, y]);
  }

  // Transformed parameters for the spatial Poisson model
  matrix[N_spatial, Y] convolved_re;
  convolved_re = sqrt(1 - rho) * rep_matrix(theta, Y) + sqrt(rho / scaling_factor) * rep_matrix(phi_spatial, Y);
}

model {
  // Priors for the multinomial model
  gamma ~ normal(0, 100);
  phi ~ normal(0, 100);
  sigma_alpha ~ inv_gamma(0.1, 0.1);
  
  // Initial alpha
  for (v in 1:V) {
    alpha[v, 1] ~ normal(gamma[v], sigma_alpha);
  }
  
  // AR1 model for alpha
  for (v in 1:V) {
    for (y in 2:Y) {
      alpha[v, y] ~ normal(gamma[v] + phi * alpha[v, (y - 1)], sigma_alpha);
    }
  }
  
  // Likelihood for the multinomial model
  for (y in 1:Y) {
    Z[, y] ~ multinomial(p[, y]);
  }

  // Model for the spatial Poisson part
  matrix[N_spatial, Y] lambda;
  for (n in 1:N_spatial) {
    for (y in 1:Y) {
      lambda[n, y] = exp(beta0 + dot_product(betas, x[n,]) + convolved_re[n, y] * sigma);
      W[n, y] ~ poisson(lambda[n, y]);
    }
  }

  // Prior for spatial effects (phi_spatial)
  target += -0.5 * dot_self(phi_spatial[node1] - phi_spatial[node2]);
  
  // Soft sum-to-zero constraint on phi_spatial
  sum(phi_spatial) ~ normal(0, 0.001 * N_spatial);

  // Priors for spatial Poisson model
  beta0 ~ normal(0.0, 1.0);
  betas ~ normal(0.0, 1.0);
  theta ~ normal(0.0, 1.0);
  sigma ~ normal(0, 1.0);
  rho ~ beta(0.5, 0.5);
}

generated quantities {
  matrix[N_spatial, Y] eta;
  matrix[N_spatial, Y] lambda_out;
  matrix[N_spatial, Y] W_fitted;
  real <lower=0> W_pred[V, N_spatial, Y];
  
  for (n in 1:N_spatial) {
    for (y in 1:Y) {
      eta[n, y] = beta0 + dot_product(betas, x[n,]) + convolved_re[n, y] * sigma;
      lambda_out[n, y] = exp(eta[n, y]);
      W_fitted[n, y] = poisson_rng(lambda_out[n, y]);
      
      // Calculate W_pred for each variant v
      for (v in 1:V) {
        W_pred[v, n, y] = p[v, y] * W_fitted[n, y];
      }
    }
  }
}
"

# Initialize our data
Y <- 4
V <- 8
Z1 <- c( 1, 1,  0,   7, 0, 0, 16, 0)
Z2 <- c( 6, 0,  2,  92, 9, 0, 42, 0)
Z3 <- c(10, 0,  4, 164, 2, 8, 12, 0)
Z4 <- c( 0, 0, 41, 283, 0, 34 ,5, 8)
Z  <- cbind(Z1, Z2, Z3, Z4)

#import shp file
library(spdep)
map <- st_read("C:/Users/chris/OneDrive - Imperial College London/Desktop/Imperial Masters 2023/MSC Statistics/Summer Research Project/SP_Municipios_2022.shp", quiet = TRUE)
nb <- poly2nb(map)

mungeCARdata4stan <- function(
    adjBUGS,
    numBUGS
) {
  J <- length(numBUGS)
  nn <- numBUGS
  J_edges <- length(adjBUGS) / 2
  node1 <- vector(mode="numeric", length=J_edges)
  node2 <- vector(mode="numeric", length=J_edges)
  iAdj <- 0
  iEdge <- 0
  for (i in 1:J) {
    for (j in 1:nn[i]) {
      iAdj <- iAdj + 1
      if (i < adjBUGS[iAdj]) {
        iEdge <- iEdge + 1
        node1[iEdge] <- i
        node2[iEdge] <- adjBUGS[iAdj]
      }
    }
  }
  return(
    list(
      "J" = J,
      "J_edges" = J_edges,
      "node1" = node1,
      "node2" = node2
    )
  )
}

adjBUGS <- unlist(nb)
numBUGS <- sapply(nb, length)
car_stan_data <- mungeCARdata4stan(adjBUGS, numBUGS)

N_spatial = car_stan_data$J;
node1 = car_stan_data$node1;
node2 = car_stan_data$node2;
#1 edge with no neighbour
N_edges = floor(car_stan_data$J_edges);

library(Matrix)
#Build the adjacency matrix 
adj.matrix = sparseMatrix(i=node1,j=node2,x=1,symmetric=TRUE)

#The ICAR precision matrix (note! This is singular)
Q=  Diagonal(N_spatial, rowSums(adj.matrix)) - adj.matrix
#Add a small jitter to the diagonal for numerical stability (optional but recommended)
Q_pert = Q + Diagonal(N_spatial) * max(diag(Q)) * sqrt(.Machine$double.eps)

# Compute the diagonal elements of the covariance matrix subject to the 
# constraint that the entries of the ICAR sum to zero.

Q_inv <- solve(Q_pert)
A <- matrix(1, 1, ncol(Q_inv))
adjustment <- Q_inv %*% t(A) %*% solve(A %*% Q_inv %*% t(A)) %*% A %*% Q_inv
Q_inv_constrained <- Q_inv - adjustment

#Compute the geometric mean of the variances, which are on the diagonal of Q.inv
scaling_factor = exp(mean(log(diag(Q_inv_constrained))))

#Import data of total cases
df= read.csv("C:/Users/chris/OneDrive - Imperial College London/Desktop/Imperial Masters 2023/MSC Statistics/Summer Research Project/Raw Data/total_cases_with_covariates_scaled.csv")

W <- cbind(df$X2020,df$X2021,df$X2022,df$X2023)
x <- cbind(df$Average.Temperature,df$Average.Rainfall)

# Prepare data for Stan
stan_data <- list(Y = Y, V = V, Z = Z, N_spatial = N_spatial, N_edges = N_edges, node1 = node1, node2=node2, W = W, K = 2, x = x, scaling_factor = scaling_factor)

# Fit our model
stan_fit <- rstan::stan(
  model_code = stan_code,
  data = stan_data,
  model_name = "Hierachical modelling of sequenced cases",
  chains = 2,
  warmup = 2000,
  iter = 4000,
  control = list(max_treedepth = 20) #, adapt_delta = 0.99)
)

samples <- rstan::extract(stan_fit)

#extract predictions for each variant and each year

variant1_2020 <- c()
for (i in 1:4000) {
  variant1_2020 <-c(variant1_2020,sum(samples$W_pred[i,1,,1]))
}
variant1_2021 <- c()
for (i in 1:4000) {
  variant1_2021 <-c(variant1_2021,sum(samples$W_pred[i,1,,2]))
}
variant1_2022 <- c()
for (i in 1:4000) {
  variant1_2022 <-c(variant1_2022,sum(samples$W_pred[i,1,,3]))
}
variant1_2023 <- c()
for (i in 1:4000) {
  variant1_2023 <-c(variant1_2023,sum(samples$W_pred[i,1,,4]))
}
variant2_2020 <- c()
for (i in 1:4000) {
  variant2_2020 <-c(variant2_2020,sum(samples$W_pred[i,2,,1]))
}
variant2_2021 <- c()
for (i in 1:4000) {
  variant2_2021 <-c(variant2_2021,sum(samples$W_pred[i,2,,2]))
}
variant2_2022 <- c()
for (i in 1:4000) {
  variant2_2022 <-c(variant2_2022,sum(samples$W_pred[i,2,,3]))
}
variant2_2023 <- c()
for (i in 1:4000) {
  variant2_2023 <-c(variant2_2023,sum(samples$W_pred[i,2,,4]))
}
variant3_2020 <- c()
for (i in 1:4000) {
  variant3_2020 <-c(variant3_2020,sum(samples$W_pred[i,3,,1]))
}
variant3_2021 <- c()
for (i in 1:4000) {
  variant3_2021 <-c(variant3_2021,sum(samples$W_pred[i,3,,2]))
}
variant3_2022 <- c()
for (i in 1:4000) {
  variant3_2022 <-c(variant3_2022,sum(samples$W_pred[i,3,,3]))
}
variant3_2023 <- c()
for (i in 1:4000) {
  variant3_2023 <-c(variant3_2023,sum(samples$W_pred[i,3,,4]))
}
variant4_2020 <- c()
for (i in 1:4000) {
  variant4_2020 <-c(variant4_2020,sum(samples$W_pred[i,4,,1]))
}
variant4_2021 <- c()
for (i in 1:4000) {
  variant4_2021 <-c(variant4_2021,sum(samples$W_pred[i,4,,2]))
}
variant4_2022 <- c()
for (i in 1:4000) {
  variant4_2022 <-c(variant4_2022,sum(samples$W_pred[i,4,,3]))
}
variant4_2023 <- c()
for (i in 1:4000) {
  variant4_2023 <-c(variant4_2023,sum(samples$W_pred[i,4,,4]))
}
variant5_2020 <- c()
for (i in 1:4000) {
  variant5_2020 <-c(variant5_2020,sum(samples$W_pred[i,5,,1]))
}
variant5_2021 <- c()
for (i in 1:4000) {
  variant5_2021 <-c(variant5_2021,sum(samples$W_pred[i,5,,2]))
}
variant5_2022 <- c()
for (i in 1:4000) {
  variant5_2022 <-c(variant5_2022,sum(samples$W_pred[i,5,,3]))
}
variant5_2023 <- c()
for (i in 1:4000) {
  variant5_2023 <-c(variant5_2023,sum(samples$W_pred[i,5,,4]))
}
variant6_2020 <- c()
for (i in 1:4000) {
  variant6_2020 <-c(variant6_2020,sum(samples$W_pred[i,6,,1]))
}
variant6_2021 <- c()
for (i in 1:4000) {
  variant6_2021 <-c(variant6_2021,sum(samples$W_pred[i,6,,2]))
}
variant6_2022 <- c()
for (i in 1:4000) {
  variant6_2022 <-c(variant6_2022,sum(samples$W_pred[i,6,,3]))
}
variant6_2023 <- c()
for (i in 1:4000) {
  variant6_2023 <-c(variant6_2023,sum(samples$W_pred[i,6,,4]))
}
variant7_2020 <- c()
for (i in 1:4000) {
  variant7_2020 <-c(variant7_2020,sum(samples$W_pred[i,7,,1]))
}
variant7_2021 <- c()
for (i in 1:4000) {
  variant7_2021 <-c(variant7_2021,sum(samples$W_pred[i,7,,2]))
}
variant7_2022 <- c()
for (i in 1:4000) {
  variant7_2022 <-c(variant7_2022,sum(samples$W_pred[i,7,,3]))
}
variant7_2023 <- c()
for (i in 1:4000) {
  variant7_2023 <-c(variant7_2023,sum(samples$W_pred[i,7,,4]))
}
variant8_2020 <- c()
for (i in 1:4000) {
  variant8_2020 <-c(variant8_2020,sum(samples$W_pred[i,8,,1]))
}
variant8_2021 <- c()
for (i in 1:4000) {
  variant8_2021 <-c(variant8_2021,sum(samples$W_pred[i,8,,2]))
}
variant8_2022 <- c()
for (i in 1:4000) {
  variant8_2022 <-c(variant8_2022,sum(samples$W_pred[i,8,,3]))
}
variant8_2023 <- c()
for (i in 1:4000) {
  variant8_2023 <-c(variant8_2023,sum(samples$W_pred[i,8,,4]))
}

variant_array <- array(NA, dim = c(32,4000)) 
variant_array[1,] <- variant1_2020
variant_array[2,] <- variant1_2021
variant_array[3,] <- variant1_2022
variant_array[4,] <- variant1_2023
variant_array[5,] <- variant2_2020
variant_array[6,] <- variant2_2021
variant_array[7,] <- variant2_2022
variant_array[8,] <- variant2_2023
variant_array[9,] <- variant3_2020
variant_array[10,] <- variant3_2021
variant_array[11,] <- variant3_2022
variant_array[12,] <- variant3_2023
variant_array[13,] <- variant4_2020
variant_array[14,] <- variant4_2021
variant_array[15,] <- variant4_2022
variant_array[16,] <- variant4_2023
variant_array[17,] <- variant5_2020
variant_array[18,] <- variant5_2021
variant_array[19,] <- variant5_2022
variant_array[20,] <- variant5_2023
variant_array[21,] <- variant6_2020
variant_array[22,] <- variant6_2021
variant_array[23,] <- variant6_2022
variant_array[24,] <- variant6_2023
variant_array[25,] <- variant7_2020
variant_array[26,] <- variant7_2021
variant_array[27,] <- variant7_2022
variant_array[28,] <- variant7_2023
variant_array[29,] <- variant8_2020
variant_array[30,] <- variant8_2021
variant_array[31,] <- variant8_2022
variant_array[32,] <- variant8_2023

my_array <- array(NA, dim = c(32,3))

for (i in 1:32) {
  quantiles <- quantile(variant_array[i,], probs = c(0.025,0.975))
  my_array[i,1] <- quantiles[1]
  my_array[i,2] <- mean(variant_array[i,])
  my_array[i,3] <- quantiles[2]
}

#quantiles for total cases
write.csv(my_array, "total_cases_quantiles.csv", row.names = FALSE)

#generate traceplots of hyperparameters
png(filename = "traceplots.png", width = 800, height = 600)

library(rstan)
traceplot(stan_fit, pars = c('gamma','phi','sigma_alpha','beta0','betas','sigma','rho'))

# Close the graphics device
dev.off()

#get predictions for each variant at municipality level
info_2 <- as.data.frame(summary(stan_fit))
info_2 <- info_2[ grepl( "W_pred", rownames(info_2)), ]

variant_1 <- info_2[ grepl( "^W_pred\\[1", rownames(info_2)), ]
variant_2 <- info_2[ grepl( "^W_pred\\[2", rownames(info_2)), ]
variant_3 <- info_2[ grepl( "^W_pred\\[3", rownames(info_2)), ]
variant_4 <- info_2[ grepl( "^W_pred\\[4", rownames(info_2)), ] 
variant_5 <- info_2[ grepl( "^W_pred\\[5", rownames(info_2)), ]
variant_6 <- info_2[ grepl( "^W_pred\\[6", rownames(info_2)), ]
variant_7 <- info_2[ grepl( "^W_pred\\[7", rownames(info_2)), ]
variant_8 <- info_2[ grepl( "^W_pred\\[8", rownames(info_2)), ]

write.csv(variant_1, file = "Variant_1.csv")
write.csv(variant_2, file = "Variant_2.csv")
write.csv(variant_3, file = "Variant_3.csv")
write.csv(variant_4, file = "Variant_4.csv")
write.csv(variant_5, file = "Variant_5.csv")
write.csv(variant_6, file = "Variant_6.csv")
write.csv(variant_7, file = "Variant_7.csv")
write.csv(variant_8, file = "Variant_8.csv")

#Load the stan model so you only need to run it once
saveRDS(object = stan_fit, file = "complete_model.rds")