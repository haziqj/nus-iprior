library(iprior)
library(tidyverse)
navyblue <- "#002f5c"

n <- 150
dat <- tibble(
  x = seq(0, 2 * pi, length = n),
  y = sin(x) #+ rnorm(n, sd = 0.01)
)
# mod <- iprior(y ~ x, dat, kernel = "se", control = list(silent = TRUE))
psi <- 1 #get_psi(mod)
lambda <-  1  # get_lambda(mod)
N <- 100
B <- 100

# I-prior sample paths ---------------------------------------------------------
prior.samp <- tibble(x = seq(0, 2 * pi, length = N))
H <- kern_se(prior.samp$x)
for (i in 1:B) {
  w <- rnorm(N, mean = 0, sd = sqrt(psi))
  prior.samp <- cbind(prior.samp, mean(dat$y) + 
                        as.numeric(lambda * H %*% w))
}
colnames(prior.samp) <- c("x", paste0("y", 1:B))
prior.samp <- reshape2::melt(prior.samp, id = "x")
# ggplot(dat, aes(x, y)) +
#   geom_point() +
#   geom_line(data = prior.samp, aes(x, value, group = variable),
#             size = 0.2, alpha = 0.4, col = "gray50") +
#   theme_classic() +
#   coord_cartesian(ylim = c(min(dat$y), max(dat$y)))

# I-prior posterior sample paths -----------------------------------------------
post.samp <- tibble(x = seq(0, 2 * pi, length = N))
H <- kern_se(dat$x)
Vy <- psi * (H * lambda) %*% (H * lambda) + diag(1/psi, n)
w_hat <- as.numeric(psi * lambda * H %*% solve(Vy) %*% (dat$y - mean(dat$y)))
w <- mvtnorm::rmvnorm(B, mean = w_hat, sigma = solve(Vy))
h <- kern_se(dat$x, post.samp$x)
for (i in 1:B) {
  post.samp <- cbind(post.samp, mean(dat$y) + as.numeric(lambda * h %*% w[i, ]))
}
colnames(post.samp) <- c("x", paste0("y", 1:B))
post.samp <- reshape2::melt(post.samp, id = "x")
# ggplot(dat, aes(x, y)) +
#   geom_point() +
#   geom_line(data = post.samp, aes(x, value, group = variable),
#             size = 0.2, alpha = 0.4, col = "gray50") +
#   theme_classic() +
#   coord_cartesian(ylim = c(min(dat$y), max(dat$y)))

bind_rows(
  bind_cols(prior.samp, type = "prior"),
  bind_cols(post.samp, type = "posterior"),
) %>%
  ggplot(aes(x, value, group = variable)) +
  geom_line(size = 0.2, alpha = 0.4, col = "gray50") +

  geom_line(data = as_tibble(prior.samp) %>% 
              group_by(x) %>% 
              summarise(value = mean(value)) %>%
              mutate(type = "prior"),
            aes(x, value), col = "red3", size = 1.5, inherit.aes = FALSE) +
  geom_line(data = as_tibble(post.samp) %>% 
              group_by(x) %>% 
              summarise(value = mean(value)) %>%
              mutate(type = "posterior"),
            aes(x, value), col = "red3", size = 1.5, inherit.aes = FALSE) +
  geom_point(data = dat, aes(x, y), inherit.aes = FALSE,
             shape = 1) +  
  facet_grid(. ~ type) +
  theme_classic() +
  coord_cartesian(ylim = c(-2, 2))

# Looking at effect of scale parameters ----------------------------------------
N <- 100
prior.samp <- tibble(x = seq(0, 2 * pi, length = N))
H <- kern_se(prior.samp$x)
lambda <- c(-2, -1, 0,  0.5, 1, 2, 5, 10)
w <- rnorm(N, mean = 0, sd = sqrt(psi))

res <- NULL
for (i in seq_along(lambda)) {
  tmp <- bind_cols(prior.samp, lambda = lambda[i],
                   y = as.numeric(lambda[i]) * H %*% w)
  res <- bind_rows(res, tmp)
}
res$lambda <- factor(res$lambda)

ggplot(res, aes(x, y, col = lambda, group = lambda)) +
  geom_line() +
  scale_color_viridis_d() +
  theme_classic() 
