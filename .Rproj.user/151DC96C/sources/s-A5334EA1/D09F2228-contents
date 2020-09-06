#'@title Main function to perform Fast functional association analysis based on wavelet
#'@description  Perform a screening for a trait on a given phenotype and a specified level of resolution
#'@param Y phenotype vector, has to be numeric. For case control code it as 0 and 1.
#'@param genotype genotype matrix (either data.frame or numeric matrix). Lines=SNPs in increasing order in term of base pair, columns=individuals. No missing values allowed.
#'@param pos vector of spatial position of the variables. It has to be in the same order as long as the number of column of the genotype matrix  line order/length.
#'@param confounder the confounding matrix with the same sample order as Y. The intercept should not be included, if missing will generate a intercept matrix.
#'@param lev_res the maximum level of resolution needed
#'@param sigma_b the parameter of the NIG prior used for the Bayes Factor computation. We advised to set this value between 0.1 and 0.2
#'@param nperm the number of permutation peformed to asses the significance of a regions. Default value 10,000
#'@param npermstop Number of permutations allowed to be larger than the observed test statistics for the regions, before stopping permuation procedure. Default value =100
#'@param para logical parameter for parallelisation, if not specified set at FALSE.
#'@details The mWaveQTL function computes the Likelihood ratio used for testing significance of a genetic region. In addition it computes the porportion of wavelets coefficients associated by level of resolution, and the Bayes factor used for this estimation.
#'@return A named vector. First position the estimated value of the Lambda statistics, then the proportion of association per level of resolution then the computed Bayes Factor per wavelet coefficient.
#'@references Shim and Stephens,Wavelet-based genetic association analysis of functional phenotypes arising from high-thoughput sequencing asssays,The Annals of Applied Statistics, 2015, Vol. 9, No. 2, 665–686
#'@export
#'@examples \dontrun{
#'
#'
#'set.seed(66)
#'#########################################
#'#Generate a randomly sample genotype size=1Mb
#'#########################################
#'
#'#5000 Randomly choosen pos
#'my_pos <- sort(sample(1:1000000, size=5000,replace = FALSE))
#'#############################
#'#Three different bump genotype
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
#'genotype <-  matrix(my_functions[my_pos,2 ], ncol=1 ) %*%t(matrix(type_fn,ncol=1))
#'#dim(genotype)= nSNP, nind
#'
#'###############################################################
#'#Generate a phenotype with variance explained by genotype  0.5%
#'###############################################################
#'varexp=0.005
#'var_noise <- (1-varexp)*var(sample(0:2,replace = TRUE,size=10000,
#'                                   prob=sampl_schem ))/varexp
#'Y <-  rnorm(n=n_size,sd=sqrt(var_noise)) +type_fn
#'df <- data.frame(y=Y,genotype =factor(type_fn))
#'P1 <- ggplot(df,aes(y=y,x=genotype))+
#'  geom_boxplot()+
#'  xlab("Type of genotype")+
#'  theme(axis.text=element_text(size=12),
#'        axis.title=element_text(size=14,face="bold"))+
#'  ylab("Simulated Phenotype")+
#'  theme_bw()+
#'  ggtitle("Variation of the phenotype\ndepending of thegenotype, \nVariance explained =0.5%")
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
#'  ggtitle("Three different kind of genotype")
#'
#'grid.arrange(P1,P2,ncol=2)
#'
#'##################
#'#Screening
#'##################
#'res <- mWaveQTL( Y,
#'                 genotype=genotype,
#'                 pos=my_pos,
#'                          lev_res=6,
#'                          sigma_b = 0.2)
#'##############
#'#Visualisation
#'##############
#'pos <- c(min(my_pos),max(my_pos))
#'plot_mWaveQTL(res=res,pos=pos,lev_res=6)
#'
#'
#'}


mWaveQTL <- function(Y,genotype,pos,confounder,lev_res,sigma_b,para=FALSE,
                     nperm=10000,
                     npermstop=100)
{
  #genotype: genotype matrix, line=SNP order in increasing pos, column individual genoype
  #pos: position of the SNP in term of base pair
  #confounder: designed matrix of the confounding effect size = n,c
  #n= n ind, c= number of confounder
  #lev_res: lev of resolution for the wavelet filtering
  #sigma_b= Para of prior, should be <1 advised 0.2


  #To ensure the length not to be 0
  Y <- as.vector(Y)


  # INPUT CHECKS
  print("Input dimensions:")
  if(!is.numeric(Y) || length(Y)==0){
    stop("ERROR: Y is not a numeric vector")
  } else {
    print(sprintf("%i phenotypes detected", length(Y)))
    if(all(Y %in% c(0,1))){
      print("Binary phenotype detected")
    } else if(!is.vector(Y)){
      stop("ERROR: Y is not a vector. Multi-phenotype analysis not implemented yet.")
    } else {
      print("Continuous phenotype detected")
    }
  }

    betas <- FALSE

  # Writing the design matrix
  if(missing(confounder)) {
    print("no covariates provided, using intercept only")
    confounder <- data.frame(confounding=rep(1,length(Y)) )
  } else if(nrow(confounder)!=length(Y)) {
    stop("ERROR: number of samples in Y and confounder does not match")
  } else {
    print(sprintf("%i covariates for %i samples detected", ncol(confounder), nrow(confounder)))
    confounder <- cbind(rep(1,length(Y)),confounder)
  }


  # Check genotype matrix
  if(is.data.frame(genotype)){
    print("Converting genotype data to matrix")
    genotype <- as.matrix(genotype)
  }
  if(missing(genotype) || !is.numeric(genotype)){
    stop("ERROR: genotype matrix missing or not numeric")
  } else if(ncol(genotype)!=length(Y)){
    stop("ERROR: number of samples in Y and genotype does not match")
  } else {
    print(sprintf("%i SNPs for %i samples detected", nrow(genotype), ncol(genotype)))
  }

  # Check position vector
  if(!is.numeric(pos) || !is.vector(pos)){
    stop("ERROR: must provide numeric position vector")
  } else {
    print(sprintf("positions for %i SNPs read", length(pos)))
  }

  # Clean missing samples from all inputs
  keepY <- complete.cases(Y)
  keepC <- complete.cases(confounder)
  keepGT <- complete.cases(t(genotype))
  nonmissing_index <- which(keepGT & keepY & keepC)
  if(length(nonmissing_index) != length(Y)){
    print(sprintf("Warning: %i individuals will be removed due to missingness",
                  length(Y) - length(nonmissing_index)))
  }

  Y <- Y[nonmissing_index]
  confounder <- confounder[nonmissing_index,]
  genotype <- genotype[,nonmissing_index]
  sigma_b <- 10*sigma_b
  print(paste("N individuals analysed = ", dim(genotype)[2],
              ", N SNPs analysed = ",dim(genotype)[1]))

  # workaround for git issue #1 - mysteriously empty slices
  if(is.null(dim(genotype)) || dim(genotype)[1] < 2^lev_res || dim(genotype)[2] < 2){
    print("Warning: not enough genotype remaining, returning empty output")

    # Naming the output
    names_BF <- c("BF_0_0")
    for(i in 1:lev_res){
      for (j in 1:(2^i)){
        names_BF <- c(names_BF,paste("BF",i,j,sep = "_"))
      }
    }
    out = rep(NA, 1+lev_res+1+length(names_BF))
    names(out) <- c("Lambda", paste("pi",0:lev_res, sep = "_"), names_BF)
    return(out)
  }

  ####################################
  #Redefinition of the needed function
  ####################################
  n_coef_wc <- function(lev_res)
  {
    temp <- c()
    for(i in 0:lev_res)
    {
      temp <- c(temp,2^i)
    }
    sum(temp)
  }

  #Quantile transform to prevent for non normaliy distrib WCs
  Quantile_transform  <- function(x)
  {
    .ex.seed <- exists(".Random.seed")
    if(.ex.seed) .oldseed <- .Random.seed
    set.seed(666)
    if(.ex.seed) on.exit(.Random.seed <<- .oldseed)


    x.rank = rank(x, ties.method="random")
    #x.rank = rank(x, ties.method="average")
    return(qqnorm(x.rank,plot.it = F)$x)
  }

  #Estimation of Lambda
  Lambda_stat <- function (my_pi, my_bayes)
  {
    # vector: pi1 pi2 pi2 pi3 pi3 pi3 pi3...
    my_pi_vec = rep(my_pi, 2^(1:length(my_pi)-1))
    coefs = 1-my_pi_vec + my_pi_vec * my_bayes[1:(2^length(my_pi)-1)]
    prod(coefs)
  }

  sumlog <- function (A1, A2)
  {
    if(A1 > A2){
      res = A1 + log(1 + exp(A2 - A1))
    }else{
      res = A2 + log(exp(A1 - A2) + 1)
    }

    return(res)
  }

  max_EM_Lambda <- function(my_bayes)
  {
    niter=10000
    epsilon <- 10^-4
    p_vec <- c()
    for(gi in 0: lev_res)
    {
      # EM algorithm for each group separately

      N_obllikli = 0
      logpi = log(0.5)
      pi <- 0.5
      log1pi = logpi

      pp = 0
      logPiBF = log(my_bayes[(2^gi):(2^(gi+1)-1)]) + logpi
      logden <- c()
      for (i in 1:length(logPiBF))
      {
        logden[i] <- sumlog(logPiBF[i],log1pi)
      }
      pp = pp+sum(exp(logPiBF - logden))
      N_obllikli = sum(logden)
      O_obllikli = N_obllikli


      for(iter in  0:niter){
        pi = pp/(2^(gi))
        logpi  = log(pi)
        log1pi = log(1-pi)
        logPiBF =   log(my_bayes[(2^gi):(2^(gi+1)-1)]) + logpi
        logden <- c()
        for (i in 1:length(logPiBF))
        {
          logden[i] <- sumlog(logPiBF[i],log1pi)
        }
        pp=0
        pp = pp+sum(exp(logPiBF - logden))
        N_obllikli = sum(logden)
        diff = abs(N_obllikli - O_obllikli)

        if(diff < epsilon){
          break
        }else{
          O_obllikli = N_obllikli
        }
      }
      p_vec <-c(p_vec,pi)
    }
    return(p_vec)
  }


  ###############
  #Paralelisation
  ###############
  if(para==TRUE)
  {
    cl <-makeCluster(detectCores(all.tests=TRUE)-1, type = "SOCK")
  }

  ###################
  #Wavelet processing
  ###################
  print("Wavelet processing")

  Time01 <- (pos- min(pos))/(max(pos)-min(pos))
  my_wavproc <- function(y)
  {
    #Kovac and Silvermann 2000
    mygrid <- wavethresh::makegrid(t=Time01,y=y)
    LDIRWD <- irregwd(mygrid,filter.number=1)
    class(LDIRWD) <- "wd"
    #Thresholding here
    LDIRWD <- threshold(LDIRWD,policy = "universal",type="hard",
                        dev = madmad,levels = 1:(LDIRWD$nlevels-1))

    res <- c()
    for(i in 0: lev_res){

      res <- c(res, accessD( LDIRWD,lev = i) )

    }

    return(res)
  }

  if(para==TRUE)
  {
    clusterExport(cl,"irregwd")
    clusterExport(cl,"threshold")
    clusterExport(cl,"madmad")
    clusterExport(cl,"accessD")
    clusterExport(cl,"accessC")
    clusterExport(cl,"my_wavproc")
    Gen_W_trans <- snow::parApply(cl,genotype,2,my_wavproc)
  }
  else{
    Gen_W_trans <- apply(genotype,2,my_wavproc)
  }

  #Quantile transform for non normal WCs for every scale location
  Gen_W_trans = apply(Gen_W_trans, 1, Quantile_transform)

  ##########
  #Modeling
  ##########
  print("Computing Bayes Factors")
  W <- as.matrix(confounder, ncol=ncol(confounder))
  n = nrow(W)
  q = ncol(W)

  # L <- as.matrix(Y , ncol=ncol(Y)) #reversed regression
  L <- as.matrix(Y,ncol=1)

  p = 1
  PW = diag(n) - W %*% solve(t(W) %*% W) %*% t(W)
  X = PW %*% L
  HB = X %*% solve(t(X) %*% X + diag(1/sigma_b/sigma_b,p)) %*% t(X)
  delta = svd(X)$d
  lambda = delta^2 / (delta^2 + 1/sigma_b/sigma_b)
  log.T = sum(log(1-lambda))/2


  my_bf <- function( y ){
    y <-  as.matrix(y,ncol=1)
    log.R = -0.5*n*log(1 - (t(y) %*% HB %*% y) / (t(y) %*% PW %*% y ))

    bf = exp(log.T + log.R)
    return(c(bf))
  }




  if(para==TRUE)
  {
    clusterExport(cl,"log.T")
    clusterExport(cl,"sigma_b")
    clusterExport(cl,"my_bf")
    my_bayes <- snow::parApply(cl,Gen_W_trans, 2, my_bf )
  }
  else{
    my_bayes <- apply(Gen_W_trans, 2, my_bf )
  }




  #################
  #Estimation Lambda
  #################
  print("Post-processing")
  my_pis <- max_EM_Lambda(my_bayes = my_bayes)
  trueLambda <- Lambda_stat(my_pi = my_pis,my_bayes = my_bayes)



  #Permutation procedure
  perm_proc <- function(x)
  {
    my_bayes_perm <-  apply(Gen_W_trans[sample(1:dim(Gen_W_trans)[1]),], 2, my_bf )
    my_pis <- max_EM_Lambda(my_bayes = my_bayes)
    trueLambda <- Lambda_stat(my_pi = my_pis,my_bayes = my_bayes)
    return(trueLambda)
  }
       batchsize   <- 100 #number of batches for permutation
       larger_perm <- 0   #number of permutation larger than the observe test statistics

 print( 'Performing permutation procedure')
  for ( perm in 1:batchsize)
  {
    print( paste ("Performing permutation procedure: batch ",perm, " out of 100" ))
    lambda_perm <-  do.call(c,lapply(floor(nperm/batchsize),perm_proc ) )

    larger_perm <- larger_perm + length(which( lambda_perm > trueLambda))

    if( larger_perm > floor(nperm/batchsize))  #Early stopping Motne Carlo pvalue
    {
      break
    }
  }

   pv <- (larger+1)/(nperm+1)

    out <- c(trueLambda,my_pis,my_bayes)

  #Naming the output
names_BF <- c("BF_0_0")
    for(i in 1:lev_res)
    {
      for (j in 1:(2^i))
      {
        names_BF <- c(names_BF,paste("BF",i,j,sep = "_"))
      }
    }

    names(out) <- c("pvalue","Lambda",
                    paste("pi",0:lev_res, sep = "_"),
                    names_BF)



  if(para==TRUE)
  {
    stopCluster(cl)
  }
  return(out)
}

