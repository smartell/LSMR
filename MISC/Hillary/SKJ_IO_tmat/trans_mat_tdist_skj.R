# Monte Carlo transition matrices for paper using IO skipjack
# now with t-distributed increment residuals
# ::
# R. Hillary Jun 2009

# call in growth parameters for skipjack

pars <- matrix(scan("params_t.dat"),ncol=3,byrow=T)

# define growth increment function

lvbinc <- function(lrel,tau,k,Linf) {

  return(max((Linf-lrel)*(1-exp(-k*tau)),0))
}

# define the length bins

lbins <- c(30,35,40,45,50,55,60,65,70,75,80)
lbins <- seq(30, 80, by=5)
# define the Monte Carlo transition matrices

nits <- dim(pars)[1]
tmat <- array(dim=c(length(lbins)-1,length(lbins)-1,nits))

# define the time interval

tau <- 1

# calculate the transition matrices

for(k in 1:nits) {

  phi <- pars[k,3]

  for(i in 1:dim(tmat)[1]) {

  lx <- lbins[i]
  ly <- lbins[i+1]
  df <- phi/mean(lbins[i]+lbins[i+1])
  epsl <- rt(1,df)

    for(j in 1:dim(tmat)[2]) {

      llj <- lbins[j]
      luj <- lbins[j+1]
          
      lli <- lx + lvbinc(lx,tau,pars[k,1],pars[k,2]) + epsl
      lui <- ly + lvbinc(ly,tau,pars[k,1],pars[k,2]) + epsl
      
      # need to work out Lebesgue measure of intersection 
      # of image and actual length bin / length bin

      if(lli > llj & lui < luj) {
        ptmp <- 1
      } else {
        tmp <- c(max(llj,lli),min(luj,lui))
        mu <- ifelse(tmp[1] < tmp[2],tmp[2]-tmp[1],0)
        nu <- lui-lli
        ptmp <- mu/nu 
      }
        
      tmat[i,j,k] <- ptmp
    }
  }
}

# calculate expected transition matrix

tmat.mu <- apply(tmat,1:2,mean)

# now using more traditional method

tmat.old <- matrix(nrow=length(lbins)-1,ncol=length(lbins)-1)

for(i in 1:nrow(tmat)) {
  for(j in 1:ncol(tmat)) {
    lref <- mean(c(lbins[i],lbins[i+1]))
    lu <- lbins[j+1]
    ll <- lbins[j]
    mul <- lref + lvbinc(lref,tau,mean(pars[,1]),mean(pars[,2]))
    df.pe <- mean(pars[,3])/lref
    tmat.old[i,j] <- pt(lu-mul,df.pe,lower.tail=T)-pt(ll-mul,df.pe,lower.tail=T)
  }
}

# %age difference in the Frobenius norms as a measure of difference

fn.new <- vector(length=nits)
for(k in 1:nits)
  fn.new[k] <- sqrt(sum(diag(tmat[,,k]%*%t(tmat[,,k]))))
fn.old <- sqrt(sum(diag(tmat.old%*%t(tmat.old))))

# rescale

fn.old <- fn.old/sqrt(ncol(tmat.old))
fn.new <- fn.new[]/sqrt(ncol(tmat.mu))

# probability that the two are different

p.diff <- length(fn.new[fn.new>fn.old])/nits

# now using E(tmat)

fn.newb <- sqrt(sum(diag(tmat.mu%*%t(tmat.mu))))/sqrt(ncol(tmat.mu))

# store MC matrices for later

save(tmat,tmat.old,file='tmat_tdist_skj.RData')
