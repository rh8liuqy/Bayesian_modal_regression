library(tidyverse)

dloglogis <- function(x,scale,shape) {
  alpha <- scale
  beta <- shape
  p1 <- (beta/alpha) * (x / alpha)^(beta - 1)
  p2 <- (1 + (x/alpha)^beta)^2
  output <- p1 / p2
  return(output)
}

xaxis <- seq(0,5,length.out = 1000)
yaxis <- dloglogis(x = xaxis, scale = 1, shape = 2)

df_loglogis <- data.frame(x = xaxis,
                          density = yaxis)

df_loglogis %>%
  ggplot(aes( x = xaxis, y = density)) +
  geom_line()


# theoritical mode: sqrt(1/3)
sqrt(1/3)
