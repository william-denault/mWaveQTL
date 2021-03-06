#'@title Zooming strategy for mWaveQTl
#'@description   perform zooming strategy on output of mWaveQTL.
#'@param res output of mWaveQTL
#'@param lev_res level of resolution of the mWaveQTL analysis
#'@return Value of the likelihood on a sub tree.
#'@examples \dontrun{
#'
#'set.seed(66)
#'#########################################
#'#Generate a randomly sample signal size=1Mb
#'#########################################
#'
#'#5000 Randomly choosen pos
#'my_pos <- sort(sample(1:1000000, size=5000,replace = FALSE))
#'#############################
#'#Three different bump signals
#'#############################
#'my_functions <-data.frame(f0 = c(rep(0,400000),rep(0,200000),rep(0,400000)),
#'                          f1 = c(rep(0,400000),rep(1,200000),rep(0,400000)) ,
#'                          f2=c(rep(0,400000),rep(2,200000),rep(0,400000)))
#'
#'
#'library(gridExtra)
#'###########################
#'#Minor allele frequency 30%
#'###########################
#'MAF=0.3
#'sampl_schem <- c((1-MAF)^2,2*MAF*(1-MAF),MAF^2)
#'#######################################
#'#sampling at Hardy Weinberg equilibrium
#'#######################################
#'#Assigning class
#'
#'#sample size =4000
#'n_size=4000
#'type_fn <-sample(0:2,replace = TRUE,size=n_size,  prob=  sampl_schem  )
#'
#'
#'signals <-  matrix(my_functions[my_pos,2 ], ncol=1 ) %*%t(matrix(type_fn,ncol=1))
#'#dim(signals)= nSNP, nind
#'
#'###############################################################
#'#Generate a phenotype with variance explained by signals  0.5%
#'###############################################################
#'varexp=0.005
#'var_noise <- (1-varexp)*var(sample(0:2,replace = TRUE,size=10000,
#'                                   prob=sampl_schem ))/varexp
#'Y <-  rnorm(n=n_size,sd=sqrt(var_noise)) +type_fn
#'df <- data.frame(y=Y,signals =factor(type_fn))
#'P1 <- ggplot(df,aes(y=y,x=signals))+
#'  geom_boxplot()+
#'  xlab("Type of signals")+
#'  theme(axis.text=element_text(size=12),
#'        axis.title=element_text(size=14,face="bold"))+
#'  ylab("Simulated Phenotype")+
#'  theme_bw()+
#'  ggtitle("Variation of the phenotype\ndepending of the signals, \nVariance explained =0.5%")
#'
#'df <- data.frame(pos= rep(my_pos,3),y=c(my_functions[my_pos,1],my_functions[my_pos,2],my_functions[my_pos,3]),
#'                 mycol = factor(c(rep("f0",length(my_pos)),rep("f1",length(my_pos)),rep("f2",length(my_pos))) ) )
#'
#'P2 <- ggplot(df,aes(y=y,x=pos,color=mycol))+
#'  geom_point(size=1)+
#'  xlab("Base pair")+
#'  ylab("Number of variants")+
#'  theme_bw()+
#'  theme(legend.title=element_blank())+
#'  ggtitle("Three different kind of signals signal")
#'
#'grid.arrange(P1,P2,ncol=2)
#'
#'##################
#'#Screening
#'##################
#'res <- mWaveQTL( Y,signal=signals,pos=my_pos,
#'                          lev_res=6,sigma_b = 0.2)
#'
#'zooming_strategy(res, 6)
#'
#'}


zooming_strategy <- function(res, lev_res )
  {
    sub <- extract_tree(res,lev_res=lev_res)
    my_pi <- adaptative_EM_Lambda(sub)
    out <-  adaptative_Lambda (my_pi, sub)
    return(out)
  }
