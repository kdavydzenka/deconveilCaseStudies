setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/deconveilCaseStudies/")

library(ggplot2)

sigma_CN <- function(CN, C50, h, A, sigma0) {
  result <- sigma0 + (A * CN^h) / (C50^h + CN^h)
  return(result)
}

#logistic_sigmoid <- function(CN, C50, h, A) {
  #return(A / (1 + exp(-h * (CN - C50))))
#}

CN_values <- seq(0, 20, by = 1.0) /2

C50 <- 2.5   # Midpoint for transition 
h <- 1.7   # Hill coefficient (controls steepness)
A <- 6.0    # Saturation effect
sigma0 <- 0.01

sigma_values <- sapply(CN_values, sigma_CN, C50 = C50, h = h, A = A, sigma0)

#sigma_values <- sapply(CN_values, logistic_sigmoid, C50 = C50, h = h, A = A)

df <- data.frame(CN = CN_values, Sigma = sigma_values)

formula_text <- expression(sigma(CN) == frac(A * CN^h, C50^h + CN^h))
CN_max = 10

ggplot(df, aes(x = CN, y = Sigma)) +
  geom_line(color = "blue", size = 1.2) +
  geom_point(data = df[df$CN %% 1 == 0, ], aes(x = CN, y = Sigma), color = "red", size = 2) +
  annotate("text", x = CN_max * 0.7, y = A * 0.5, label = formula_text, parse = TRUE, size = 6) +
  scale_x_continuous(limits = c(0, 10), breaks = 1:10) +
  scale_y_continuous(limits = c(0, 6), breaks = 1:6) +
  labs(title = "Saturation model",
       x = "CN/2 values",
       y = "Sigma(CN)") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, face = "plain"))

