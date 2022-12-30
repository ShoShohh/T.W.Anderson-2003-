library(MASS)

T_Anderson <- function(N, L, Lam){
  
 return (sqrt((N - 1) / 2) * (L - Lam) / Lam)

}

AndersonHist <- function(Mu, Sigma, N, Trial){
  
  Lam <- eigen(Sigma)$values[1] # 母固有値

  L <- numeric(100 * Trial) # 標本を抽出する度にその固有値を入れるためのベクトル
  Anderson <- numeric(100 * Trial) # (100 * Trial)の回数だけAnderson()の結果を得て，
                                   # それを入れるベクトル

    for(i in 1:(100 * Trial)){ # (100 * Trial)の回数だけAnderson()の結果を得る．
                               # その結果をAndersonに蓄積する．

      Sam <- mvrnorm(N, Mu, Sigma) #多次元正規分布に従う独立な確率ベクトルをN組発生．
                                   #ただしSigmaは対称行列で正定値行列
      L[i] <- (prcomp(Sam)$sdev[1]) ^ 2 # 標本固有値
      Anderson[i] <- T_Anderson(N, L[i], Lam) # 標本数Nをi回目に抽出して計算したT_Anderson()を
                                              # Anderson[i]に代入する．

    }

  df <- data.frame(Anderson)

  library(ggplot2)
  ggplot(data = df, aes(x = Anderson)) + 
    geom_histogram()
}

Mu0　<-　c(1, 2, 3, 4, 5) #母期待値ベクトル
Sigma0 <- rbind( #母共分散行列，ただし正定値
  c(1, 0, -1, 1, -1),
  c(0, 2, 0, 1, 1),
  c(-1, 0, 3, 0, 2),
  c(1, 1, 0, 4, -1),
  c(-1, 1, 2, -1, 5)
)
AndersonHist(Mu0, Sigma0, 100, 1)