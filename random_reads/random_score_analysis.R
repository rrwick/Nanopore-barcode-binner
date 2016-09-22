
library(ggplot2)
library(scales)

random_scores <- read.delim("~/Programs/nanopore_barcode_binner/random_reads/random_scores")


log_log_breaks <- function(x) {
  b <- 10^10^(seq(log10(log10(min(x))), log10(log10(max(x))), length.out=8))
  return(signif(b,2))
}


double.log.plot.regression <- function(random_scores, col.name) {
  fit <- lm(unlist(random_scores[col.name]) ~ log10(log10(unlist(random_scores['Read.length']))))
  x <- seq(log10(log10(24.0)), log10(log10(1000000.0)), 0.01)
  y <- fit$coefficients[1] + fit$coefficients[2] * x
  x <- 10^(10^x)
  regression = data.frame(x=x, y=y)
  
  print(paste(c("slope:", fit$coefficients[2]), collapse = " "))
  print(paste(c("intercept:", fit$coefficients[1]), collapse = " "))
  print(paste(c("score =", fit$coefficients[2], "* log10(log10(read_length)) -", -fit$coefficients[1]), collapse = " "))

  # ggplot(random_scores, aes_string(x='Read.length', y=col.name)) +
  #   geom_point(size=1, alpha=0.2, color='red') +
  #   geom_line(data=regression, aes(x=x, y=y)) +
  #   theme_bw() +
  #   scale_y_continuous(expand=c(0,0), limits = c(-30, 50)) +
  #   scale_x_continuous(expand=c(0,0), limits = c(24, 1000000))
# 
#   ggplot(random_scores, aes_string(x='Read.length', y=col.name)) +
#     geom_point(size=1, alpha=0.2, color='red') +
#     geom_line(data=regression, aes(x=x, y=y)) +
#     theme_bw() +
#     scale_y_continuous(expand=c(0,0), limits = c(-30, 50)) +
#     scale_x_log10(expand=c(0,0), limits = c(24, 1000000))

  ggplot(random_scores, aes_string(x='Read.length', y=col.name)) +
    geom_point(size=1, alpha=0.2, color='red') +
    geom_line(data=regression, aes(x=x, y=y)) +
    theme_bw() +
    scale_y_continuous(expand=c(0,0), limits = c(-30, 50)) +
    scale_x_continuous(expand=c(0,0), limits = c(24, 1000000), trans=trans_new("log-log", function(x) log10(log10(x)), function(x) 10^10^x), breaks=log_log_breaks)
}





double.log.plot.regression(random_scores, 'Barcode.1.score')
double.log.plot.regression(random_scores, 'Barcode.2.score')
double.log.plot.regression(random_scores, 'Barcode.3.score')
double.log.plot.regression(random_scores, 'Barcode.4.score')
double.log.plot.regression(random_scores, 'Barcode.5.score')
double.log.plot.regression(random_scores, 'Barcode.6.score')
double.log.plot.regression(random_scores, 'Barcode.7.score')
double.log.plot.regression(random_scores, 'Barcode.8.score')
double.log.plot.regression(random_scores, 'Barcode.9.score')
double.log.plot.regression(random_scores, 'Barcode.10.score')
double.log.plot.regression(random_scores, 'Barcode.11.score')
double.log.plot.regression(random_scores, 'Barcode.12.score')

ggplot(random_scores, aes(x=Read.length, y=Barcode.1.score)) +
  geom_point(size=1, alpha=0.2) +
  theme_bw() + 
  scale_y_continuous(expand=c(0,0), limits = c(-30, 50)) + 
  scale_x_log10(expand=c(0,0), limits = c(24, 1000000))


ggplot(random_scores, aes(x=Read.length, y=Barcode.2.score)) +
  geom_point(size=1, alpha=0.2) +
  theme_bw() + 
  scale_y_continuous(expand=c(0,0), limits = c(-30, 50)) + 
  scale_x_log10(expand=c(0,0), limits = c(24, 1000000))
