#'@title Simulation of Lambda null distribution on a sub tree
#'@description internal function for qimulation of Lambda null distribution.
#'@param nsimu number of simulation to perform.
#'@param lambda the lambda parameter of the null distribution of the Bayes Factor.
#'@param lev_res the maximum level of resolution needed.
#'@param ncp the lambda parameter of the null distribution of the Bayes Factor, if not specified set at 0.
#'@details Use the theoretical null distribution of the generated Bayes factor to perform simulation of null distriution of the test statistics of the Wavelet screaming procedure.
#'@references Quan Zhou and Yongtao Guan, On the Null Distribution of Bayes Factors in linear Regression, Journal of the American Statistical Association, 518, 2017.










adaptative_Simu_Lambda_null <- function(nsimu,lambda,lev_res,ncp,sub,thresh)

{

  if(missing(thresh))
  {
    thresh <-1
  }
  if(length( which( as.numeric(sub)>thresh ) )==1 )
  {
    print( "Only one Bayes Factor above the choosen threshold(set as one per deault), analytical formula available. Please use the function analytical_p")
    break
  }


  if(missing(ncp)){
    ncp=0
  }


  sumlog <- function (A1, A2)
  {

    if(A1 > A2){
      res = A1 + log(1 + exp(A2 - A1))
    }else{
      res = A2 + log(exp(A1 - A2) + 1)
    }

    return (res)

  }


  #Adapted version that match the sub structure
  loc_EM_Lambda <- function(my_bayes)
  {
    niter=10000
    epsilon <- 10^-4
    p_vec <- c()

    BF_class <-  as.numeric(gsub("^[^_]*_|_[^_]*$", "", colnames(sub )))
    for(gi in unique(BF_class))
    {
      # EM algorithm for each group separately

      N_obllikli = 0
      logpi = log(0.5)
      pi <- 0.5
      log1pi = logpi

      pp = 0
      logPiBF = log(my_bayes[which(BF_class==gi)]) + logpi
      logden <- c()
      for (i in 1:length(logPiBF))
      {
        logden[i] <- sumlog(logPiBF[i],log1pi)
      }
      pp = pp+sum(exp(logPiBF - logden))
      N_obllikli = sum(logden)
      O_obllikli = N_obllikli


      for(iter in  0:niter){
        pi = pp/(length(which(BF_class ==gi)))
        logpi  = log(pi)
        log1pi = log(1-pi)
        logPiBF =   log(my_bayes[which(BF_class==gi)]) + logpi
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

  nBF <- length(as.numeric(sub))
  BF <- matrix(NA,ncol=nBF,nrow=nsimu)


  print("Simulation Bayes Factors")
  for (i in  1:nsimu)

  {

    for (j in 1: nBF)
    {
      BF[i,j] <- exp( (lambda*rchisq( n=1 ,ncp=ncp, df=1)+log(1-lambda) )/2 )
    }
  }
  loc_Lambda_stat <- function (my_pi, my_bayes,sub)
  {
    BF_class <-  as.numeric(gsub("^[^_]*_|_[^_]*$", "", colnames(sub )))

    my_pi_vec <- c()
    for ( i in 1: length(unique(BF_class)))
    {
      temp <- rep(my_pi[i], c( length( which(BF_class == unique(BF_class)[i] ) ) ) )
      my_pi_vec <- c(my_pi_vec, temp )
    }
    coefs = 1-my_pi_vec + my_pi_vec * my_bayes

    prod(coefs)
  }

  print("Lambda maximisation")
  my_pis<- t(apply(BF,1,loc_EM_Lambda))

  Simu_Lambda <- c()
  for (i in 1:dim(my_pis)[1])
  {
    Simu_Lambda <-c( Simu_Lambda , Lambda_stat(my_bayes=BF[i,],my_pi=my_pis[i,]))

  }


  return(Simu_Lambda)





}








#'@title Likelihood ratio for sub tree
#'@description   internal function objective function for sub trees.
#'@param my_pi a vector of the proportion of association per level of resolution.
#'@param sub Output of extract_treet.
#'@return Value of the likelihood on a sub tree.
#'@examples \dontrun{
#'
#'#using res for the Wavelet_screaming exemple
#'
#'
#'sub_analysis <- function(res, lev_res )
#'{
#'  sub <- extract_tree(res,lev_res=lev_res)
#'  my_pi <- adaptative_EM_Lambda(sub)
#'  out <-  adaptative_Lambda (my_pi, sub)
#'  return(out)
#'}
#'
#'
#'sub_analysis(res, 6)
#'
#'}



adaptative_Lambda <- function(my_pi,sub)
{

  my_bayes <- as.numeric(sub)
  BF_class <-  as.numeric(gsub("^[^_]*_|_[^_]*$", "", colnames(sub )))

  my_pi_vec <- c()
  for ( i in 1: length(unique(BF_class)))
  {
    temp <- rep(my_pi[i], c( length( which(BF_class == unique(BF_class)[i] ) ) ) )
    my_pi_vec <- c(my_pi_vec, temp )
  }
  coefs = 1-my_pi_vec + my_pi_vec * my_bayes

  prod(coefs)
}



#'@title  EM procedure for sub tree
#'@description internal function for EM procedure for sub tree
#'@param sub a Output of extract_tree
#'@return Estimated proportion for the sub tree
#'@examples \dontrun{
#'
#'#using res for the Wavelet_screaming exemple
#'
#'
#'sub_analysis <- function(res, lev_res )
#'{
#'  sub <- extract_tree(res,lev_res=lev_res)
#'  my_pi <- adaptative_EM_Lambda(sub)
#'  out <-  adaptative_Lambda (my_pi, sub)
#'  return(out)
#'}
#'
#'
#'sub_analysis(res, 6)
#'
#'}




adaptative_EM_Lambda <- function(sub)
{
  niter=10000
  epsilon <- 10^-4
  p_vec <- c()

  my_bayes <- as.numeric(sub)
  BF_class <-  as.numeric(gsub("^[^_]*_|_[^_]*$", "", colnames(sub )))
  for(gi in unique(BF_class))
  {
    # EM algorithm for each group separately

    N_obllikli = 0
    logpi = log(0.5)
    pi <- 0.5
    log1pi = logpi

    pp = 0
    logPiBF = log(my_bayes[which(BF_class==gi)]) + logpi
    logden <- c()
    for (i in 1:length(logPiBF))
    {
      logden[i] <- sumlog(logPiBF[i],log1pi)
    }
    pp = pp+sum(exp(logPiBF - logden))
    N_obllikli = sum(logden)
    O_obllikli = N_obllikli


    for(iter in  0:niter){
      pi = pp/(length(which(BF_class ==gi)))
      logpi  = log(pi)
      log1pi = log(1-pi)
      logPiBF =   log(my_bayes[which(BF_class==gi)]) + logpi
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

#'@title Extract sub result from result of the mWaveQTL
#'@description  internal function forto performed zoomed analysis of the mWaveQTL output
#'@param res Output of Wavelet_screaming.
#'@param lev_res the maximum level of resolution needed, has to be less or equal to the request level of resolution in the Wavelet_screaming.
#'@param thresh Minimal value of the Bayes Factor to  defined the a sub region, if missing set as 1.
#'@return A vector correpsonding of the sub tree for the zoomed analysis.
#'@examples \dontrun{
#'
#'#using res for themWaveQTL exemple
#'
#'
#'sub_analysis <- function(res, lev_res )
#'{
#'  sub <- extract_tree(res,lev_res=lev_res)
#'  my_pi <- adaptative_EM_Lambda(sub)
#'  out <-  adaptative_Lambda (my_pi, sub)
#'  return(out)
#'}
#'
#'
#'sub_analysis(res, 6)
#'
#'}



extract_tree <- function(res,lev_res,thresh)
{

  if(missing(thresh))
  {
    thresh <- 1
  }
  res <- res[-c(1:(lev_res+2))]
  temp <- colnames(res)[which(as.numeric(res)>thresh)]

  #get the Scale
  my_s <- as.numeric(gsub("^[^_]*_|_[^_]*$", "", temp))

  cor <- c()
  for (i in 1:length(my_s))
  {
    temp <- sum(2^(0:(my_s[i]-1)))
    cor <- c(cor, temp)
  }
  #Get the location of the BF over thresh
  my_l <-   which(as.numeric(res)>thresh)-cor

  #check if the lower scale are nested into the upper scale

  temp <- which(my_s >min(my_s))

  #Case where ther is only one Bayes Factor above 1/ all at the same scale
  if(length(temp)>0)
  {
    temps <- my_s[temp] - min(my_s)

    #my_l essemble of the location needed at the min scale
    cors <- c()
    for( i in 1: length(temp))
    {
      temp1 <- my_l[temp[i]]

      #divided enought time to get the corresponding scale location
      for(j in 1: temps[i])
      {
        temp1 <- ceiling(temp1/2)
      }
      cors <- c(cors,temp1)
    }
    #Corresponding WC at min scale
    locations <- c( my_l[which(my_s==min(my_s))],  cors)
  }
  else{
    locations <- c( my_l[which(my_s==min(my_s))])
  }

  #Index to select
  my_index <- c()

  span <- min(locations):max(locations)
  for ( i in min(my_s):lev_res)
  {
    temp <- sum(2^( 0:( min(i) -1) ) )+ span
    span <- (2*min(span)-1):(2*max(span))
    my_index <- c( my_index,temp)

  }

  return(res[my_index])

}


#'@title Defintion of the regions with Bayes factor over a threshold.
#'@description  internal function for fine mapping tool for output of the mWaveQTL function.
#'@param res Output of mWaveQTL.
#'@param lev_res the maximum level of resolution of the previous analysis.
#'@param thresh numeric, Bayes factor threshold to defined the fine mapping. If missing set as 1.
#'@param start numeric, start in base pair of the analyzed regions .
#'@param chr numeric, end in base pair of the analyzed regions .
#'@details return a list of chr, start, end position, that correspond of the sub regions defined by the dyadic decomposition of wavelet that are associated with a Bayes factor over the defined threshold.


fine_map <- function(res,lev_res,thresh,start,end,chr)
{
  if(missing(thresh))
  {
    thresh=1

  }
  k <-1
  temp <- lev_res+2
  l1 <- which(res[- c(1:temp )] >thresh )
  fmap <- list()
  for (i in 1: lev_res)
  {
    lt <- l1[which(l1 < sum(2^(0:i)))]
    lt <-  lt[which(lt> sum(2^(0:(i-1 ) )  ) ) ]



    if(length(lt) >0 )
    {
      ltt <- ( lt- sum(2^(0:(i-1 ) ) ) )
      for ( j in 1: length(ltt))
      {
        startpos <- as.numeric(start)-1 +
          (ltt[j] -1)*(1/(2^i)) *(as.numeric(end) -as.numeric(start) )
        endpos <- as.numeric(start)-1 +
          (ltt[j] )*(1/(2^i)) *(as.numeric(end) -as.numeric(start)  )

        chr <- chr

        fmap[[k]] <-  c( chr, startpos,endpos)

        k <- k+1
      }



    }



  }
  return(fmap)
}





#'@title Compute the relative size of a sub region.
#'@description  internal function to compute the relative size of a sub region.
#'@param sub Output of extract_tree.
#'@return A proportion, used for multiple testing correction for the region






prop_sub <- function(sub)
{
  my_s <- as.numeric(gsub("^[^_]*_|_[^_]*$", "", colnames(sub )))
  out <- length(which(my_s ==min(my_s)))/   2^(min(my_s))
  return(out)

}
