T_Anderson <- function(SAMPLE_NUM, SAMPLE_EV, POP_EV){
  
 return (sqrt((SAMPLE_NUM - 1) / 2) * (SAMPLE_EV - POP_EV) / POP_EV)

}

AndersonHist <- function(POP_NUM, SAMPLE_NUM, mu, sd, p){

  POP <- matrix(rnorm(POP_NUM * p, mu, sd), ncol = p)
  POP_EV <- (prcomp(POP)$sdev[1]) ^ 2

  Anderson <- numeric(100)

  for(i in 1:100){
    SAMPLE <- POP[sample(1:POP_NUM, SAMPLE_NUM), ]
    SAMPLE_EV <- (prcomp(SAMPLE)$sdev[1]) ^ 2
    Anderson[i] <- T_Anderson(SAMPLE_NUM, SAMPLE_EV, POP_EV)
  }

  df <- data.frame(Anderson)

  library(ggplot2)
  g <- ggplot(data = df, aes(x = Anderson)) 
  g <- g + geom_histogram()
  plot(g)
}

AndersonHist(1000, 500, 5, 10, 5)
