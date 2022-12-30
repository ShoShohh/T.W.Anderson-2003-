library(MASS) # 多変量正規分布に従う乱数を生成するmvrnorm()関数を用いるため，
              # それが定義されているMASSライブラリを読み込む．

# T.W.Andersonで導出された検定統計量を用いて信頼区間を構成し，母固有値に対して漸近的な受容域(4)の検証を行う関数
Anderson <- function(A, N, L, Lam){      # 引数は有意水準A/100のA%，標本数N，標本固有値L，母固有値Lam

  z <- qnorm(1 - (A / 2) / 100, 0, 1)　     # 100(A/2/100)%点zの導出
  T <- abs(sqrt((N - 1) / 2) * (L - Lam) / Lam) # T.W.Andersonで与えられる検定統計量の絶対値
  
 return(if(T <= z) {1} else {0}) # T.W.Andersonの不等式(4)に対応している．
                                 # 与えられた母固有値と標本固有値が(4)を満たせば1，満たさなければ0を返す．

}

# 期待値ベクトルMu，共分散行列Sigmaの多変量正規分布に従う母集団から，無作為抽出する標本の数を増やしていき，
# 各標本数毎に(100 * Trial)回のAnderson()関数を用いた信頼区間の検証を行う．
# (100 * Trial)回はその都度標本を取り直す．
# 信頼区間の検証とはつまり，抽出した標本の集まり毎の標本固有値で構成されるAndersonの信頼区間に母固有値が
# 含まれているかを確認し，さらに構成された信頼区間が母固有値の100(1-A/100)%信頼区間
# (Aは百分率として与えられるものとする)として成り立っているのかを確認し，
# 各標本数毎に構成した信頼区間に母固有値が含まれた回数をプロットする関数
ConfiPlot <- function(A, Mu, Sigma, Max, Trial){ # 引数は有意水準A/100のA%，母集団が従う多変量正規分布の
                                                 # 期待値ベクトルMuと共分散行列Sigma，抽出する標本数の最大値Max，
                                                 # 各標本数で試行する回数(100 * Trial)のTrial
  NVec <- c(2:Max) # 抽出する標本数を並べるベクトル．
                   # NVec[i] = 1のとき検定統計量が0になってしまうことを避ける．
  CountsVec <- numeric(Max - 1) # 構成した信頼区間に母固有値が入った回数(実際はTrialで割る)を並べるベクトル

  Lam <- eigen(Sigma)$values[1] # 母固有値
  
  Guide <- 0 # 考えている有意水準A/100に対して，
             # 上手く100(1-A/100)%信頼区間が構成できたときの標本数を入れるオブジェクト

  for(i in 1:(Max - 1)){ # 抽出する標本数の最大値(Maxで与えられる)-1の回数だけ繰り返す．
                         # つまり各標本数で以下を実行する．

    Counts <- 0  # 各標本数において構成した信頼区間に母固有値が入った回数を入れるオブジェクト

    for(j in 1:(100 * Trial)){ # (100 * Trial)の回数だけ繰り返し各標本数で信頼区間を構成する．
                               # Anderson()の結果がCountsに蓄積される．

      Sam <- mvrnorm(NVec[i], Mu, Sigma) # 多次元正規分布に従う確率ベクトルをNVec[i]組発生
                                         # ただしSigmaは対称行列で正定値行列
      L <- (prcomp(Sam)$sdev[1]) ^ 2     # 標本固有値
      Counts <- Counts + Anderson(A, NVec[i], L, Lam) # 標本数NVec[i]をj回目に抽出して計算したAnderson()が
                                                      # j-1回目までのCountsに足される．

    }

    CountsVec[i] <- Counts / Trial # CountsVec[i]に標本数NVec[i]の時のCountsを代入する．
    if(((CountsVec[i]) >= 100 - A) && (Guide == 0)){ # 満たしたい有意水準A/100に対して，
                                                     # 100(1-A/100)%信頼区間が構成出来た時の
                                                     # 標本数をGuideに代入する．
       Guide <- NVec[i]
    }
  }

  plot(NVec, CountsVec, xlim = c(0, Max), ylim = c(0, 100) # 横軸を標本数NVec，縦軸を各標本数において構成した
                                                           # 信頼区間に母固有値が入った回数CountsVecとしたグラフ
       , xaxt = "n", yaxt = "n", xlab = "抽出した標本の数(NVec)"
       , ylab = "構成した信頼区間に母固有値が入った回数/Trial(CountsVec)", pch = 1)
  abline(h = 100 - A, col = 'red')                         # 横軸に平行で．満たしたい(信頼係数*100)%に線を引く．

  axis(side = 2, at = c(0, 50, 100), labels = c(0, 50, 100), cex.axis=0.6)
  axis(side = 2, at = 100 - A,　labels = 100 - A, col.ticks ='red', col.axis = "red")

  axis(side = 1, at = c(0, Max / 2, Max), labels = c(0, Max / 2, Max), cex.axis=0.6)
  if (Guide != 0){           # Guideの値が標本数の最大値として与えられるMax以下であれば，縦軸に平行で
                             # 満たしたい有意水準A/100に対して100(1-A/100)%信頼区間を
                             # 最初に構成出来た標本数で線を引く．
     abline(v = Guide, col = 'black')
     axis(side = 1, at = Guide,　labels = Guide, col.ticks ='black', col.axis = "black")
  }
}


Mu0　<-　c(1, 2, 3, 4, 5) #母期待値ベクトル
Sigma0 <- rbind(         #母共分散行列，ただし正定値
  c(1, 0, -1, 1, -1),
  c(0, 2, 0, 1, 1),
  c(-1, 0, 3, 0, 2),
  c(1, 1, 0, 4, -1),
  c(-1, 1, 2, -1, 5)
)

# 上記の母期待値ベクトルMu0，母共分散行列Sigma0に従う確率ベクトルから標本を最大でMax = 100まで
# 抽出し，標本数2,...,Maxそれぞれの下で(100 * Trial) = 100の回数だけAndersonを計算する．
# 有意水準はA/100 = 0.1（A=10%）で，100(1-A/100) = 90%信頼区間を構成したい．
ConfiPlot(10, Mu0, Sigma0, 100, 1) 