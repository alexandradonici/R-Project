# GRUPA 244
# Donici Alexandra-Maria
# Muscalu Elena-Claudia


#Problema 1

# 1. Generarea a 1000 de realizari independente 

b <- rbinom(1000,size=10,p=0.7)
p <- rpois(1000, lambda=3)
e <- rexp(1000,3)
n <- rnorm(1000, mean=0, sd=1)

#varianta esantionului
variance <- function(vec)
{
  range(vec)
  counts <- hist(vec + .01, breaks = 7)$counts
  fx <- counts / (sum(counts))
  x <- c(min(vec): max(vec))
  exp <- sum(x * fx)  #media
  exp.square <- sum(x^2 * fx) 
  var <- exp.square - (exp)^2 ; var #varianta
}

#media esantionului
mean <- function(vec)
{
  range(vec)
  counts <- hist(vec + .01, breaks = 7)$counts
  fx <- counts / (sum(counts))
  x <- c(min(vec): max(vec)) 
  exp <- sum(x * fx) ; exp 
}

#media si varianta repartitiei binomiale
mean(b)
variance(b)

#media si varianta repartitiei Poisson
mean(p)
variance(p)


#media si varianta repartitiei exponentiale
mean(e)
variance(e)

#media si varianta repartitiei uniforme
mean(n)
variance(n)

###############################################################

# 2. Functiile de MASA si de DENSITATE

fct2 <- function(a,b,n,p)
{
  lambda = n*p
db <- dbinom(a:b, size = n, prob =p)
dp <- dpois(a:b, lambda)
de <- dexp(a:b,lambda)
dn <- dnorm(a:b, mean=n, sd=p)
plot(c(a:b,a:b,a:b,a:b), c(db,dp,de,dn),type="n", ylim=c(0, 0.5), main="Functia de masa")

points((a:b)+.15, de, type = "h",
       col ="orange", lwd = 8)
points((a:b)-.15, db, type = "h",
       col ="blue", lwd = 8)
points((a:b)+.15, dp, type = "h",
       col ="red", lwd = 8)

points((a:b)+.15, dn, type = "h",
       col ="black", lwd = 8)
  legend(b-b/2, 0.5, legend = c(paste0("Binomiala(",n,",",p,")"),paste0("Poisson(",lambda,")"), paste0("Normala(",lambda,")"),paste0("Binomiala(",n,",",p,")")), fill = c("blue", "red", "orange", "black"), bg="white", bty = "n")
}

fct2(0,20,20,0.5)
fct2(0,7,10,0.1)
fct2(0,7,15,0.3)
fct2(0,7,8,0.8)
fct2(0,10,5,1)

######################################################################################
# 3. Functiile de repartitie
fct3 <- function(a,b,n,p)
{
  lambda = n*p
  db <- pbinom(a:b, size = n, prob =p)
  dp <- ppois(a:b, lambda)
  de <- pexp(a:b,lambda)
  dn <- pnorm(a:b, mean=n, sd=p)
  plot(c(a:b,a:b,a:b,a:b), c(db,dp,de,dn),type="n", ylim=c(0, 2.5), main="Functia de repartitie")
  
  points((a:b)+.15, de, type = "h",
         col ="orange", lwd = 8)
  points((a:b)-.15, db, type = "h",
         col ="blue", lwd = 8)
  points((a:b)+.15, dp, type = "h",
         col ="red", lwd = 8)
  
  points((a:b)+.15, dn, type = "h",
         col ="black", lwd = 8)
  legend(b-b/2, 2.5, legend = c(paste0("Binomiala(",n,",",p,")"),paste0("Poisson(",lambda,")"), paste0("Normala(",lambda,")"),paste0("Binomiala(",n,",",p,")")), fill = c("blue", "red", "orange", "black"), bg="white", bty = "n")
}

fct3(0,20,20,0.5)
fct3(0,7,10,0.1)
fct3(0,7,15,0.3)
fct3(0,7,8,0.8)
fct3(0,10,5,1)

###################################################################################

# 4. Tabel - aproximari functie binomiala

tabel4 = function(n,p,a,b)
{
  lambda = n*p
  x = matrix(numeric((b-a+1)*6),ncol=6,
             
             dimnames = list(a:b, c("k","Binomiala","Poisson","Normala","Normala Corectie", "Camp-Paulson")))
  
  for(k in a:b)
    x[k][1] = k
  
  
  
  x[,2] = pbinom(a:b,n,p)
  x[,3] = ppois(a:b,lambda)
  
  
  
  for( k in a:b)
  {
    ab <- (k-n*p)/sqrt(n*p*(1-p))
    ac <- (k+0.5-n*p)/sqrt(n*p*(1-p))
    x[k,4] = pnorm(ab,0,1)
    x[k,5] = pnorm(ac,0,1)
  }
  
  
  
  for(k in a:b)
  {
    a1 = 1/(9*(n-k))
    b1 = 1/(9*(k+1))
    r = ((k+1)*(1-p))/(p*(n-k))
    c = (1-b1)*r^(1/3)
    miu = 1-a1
    sigma = sqrt(a1+ (b1*r^(2/3)))
    ad <- (c-miu)/sigma
    x[k,6] = pnorm(ad)
    
  }
  
  error = max(abs(x[,6]))
  return(list(x = as.data.frame(x), error = error, param = c(n, p, lambda)))
}

#apeluri pentru n apartine {25,50, 100} si p apartine {0.05, 0.1} si k apartine [a,b]
tabel4(25,0.05,1,10)
tabel4(25,0.1,1,10)
tabel4(50,0.05,1,10)
tabel4(50,0.1,1,10)
tabel4(100,0.05,1,10)
tabel4(100,0.1,1,10)


#####################################################################
#Problema 2

# a) functia fgam

fgam <- function(a)
{
  f2 <- function(x)
  {
    x^(a-1)*exp(-x)
  }
  if(a<=0)
    return(-1)
  
  if(a==1)
    return(1)
  if(a==1/2)
    return(sqrt(pi))
  
  if(is.integer(a)==TRUE && a>0)
    return(factorial(a-1))
  
  if(a>1)
    return((a-1)*fgam(a-1))
  
  return (integrate(f2,0,Inf))
  
}
##############################################################

# b) functia fbet
fbet <- function(a,b)
{
  if(a+b==1 && a>0 && b>0)
    return(pi/sin(a*pi))
  return((fgam(a)*fgam(b))/fgam(a+b))
    
}
#################################################################

fGamma <- function(x,a,b)
{
  if(x>0 && (a>0 && b>0))
    fG=1/(b^a*fgam(a))*x^(a-1)*(exp(1))^(-x/b)
  return (fG)
}

fBeta <-  function(x,a,b)
{
  #densitatea pentru Y ~ Beta
  if((a>0 && b>0) && (x<1 && x>0))
    fB= 1/(fbet(a,b))*x^(a-1)*(1-x)^(b-1)
  return (fB)
  
}



FGamma1 <- function(x,miu,lambda)
{
  f <- function(u)
  {
    u^(miu-1)*exp((-1)*lambda*u)
  }
  z <- integrate(f,0,x)
  y <- ((lambda^miu) / fgam(miu))*z$value
  if(x>0)
    return (y)
  if(x<=0)
    return (0)
}

FBeta1 <- function(x, alfa, beta)
{
  if(x<=0)
    return (0)
  
  if(x>=1)
    return (1)
  
  f <- function(u)
  {
    u^(alfa-1)*(1-u)^(beta-1)
  }
  
  z <- integrate(f,0,x)
  y <- fgam(alfa+beta)/(fgam(alfa)*fgam(beta))*z$value
  return (y)
}
#### 1) P(X < 3)

fprobgamma1 <- function(a,b)
{
  return (FGamma1(3,a,b))
}

### apel pentru 3 seturi de valori
fprobgamma1(2,1)
fprobgamma1(3,2)
fprobgamma1(4,1)

####


#### 2) P(2 < X < 5)
fprobgamma2 <- function(a,b)
{
  return (FGamma1(5,a,b)-FGamma1(3,a,b))
}

fprobgamma2(2,1)
fprobgamma2(3,2)
fprobgamma2(3,1)

### 3) P(3 < X < 4| X>2)
fprobgamma3 <- function(a,b)
{
  m=FGamma1(4,a,b)-FGamma1(3,a,b)
  n= 1-FGamma1(2,a,b)
  return (m/n)
}

fprobgamma3(2,1)
fprobgamma3(3,2)
fprobgamma3(3,1)


#### 4) P(Y > 2)
fprobbeta4 <- function(a,b)
{
  return (1-FBeta1(2,a,b))
}

fprobbeta4(2,1)
fprobbeta4(3,2)
fprobbeta4(3,1)


#### 5) P(4 < X < 6)
fprobgamma5 <- function(a,b)
{
  return (FGamma1(6,a,b)-FGamma1(4,a,b))
}

fprobgamma5(2,1)
fprobgamma5(3,2)
fprobgamma5(3,1)

#### 6) P(0< X < 1 | X < 7)
fprobgamma6 <- function(a,b)
{
  x <- FGamma1(1,a,b)-FGamma1(0,a,b)
  y <- FGamma1(7,a,b)
  return (x/y)
}

fprobgamma6(2,1)
fprobgamma6(3,2)
fprobgamma6(3,1)


#### 7) P(X + Y < 5)
fprob7 <- function(a, b)
{
  
  f <- function(x, y)
  {
    #densitatea pentr X~ Gamma
    fG <- function(x)
    {
      fG=1/(b^a*fgam(a))*x^(a-1)*(exp(1))^(-x/b)
      return (fG)
    }
    
    #densitatea pentru Y ~ Beta
    fB <-  function(x)
    {
      
      fB= 1/(fbet(a,b))*x^(a-1)*(1-x)^(b-1)
      return (fB)
      
    }
    
    #densitatea comuna f(x,y)
    return (fG(x)* fB(y))
  }
  
  
  ymax <-  function(x)
  {
    # stim ca y apartine(0,1) functia Beta
    min <-  min(1,5-x) 
  }
  
  I <- integral2(f,0,5,0, ymax) #pracma
  return (I)
  
}
fprob7(2,1)

####8) P(X-Y>0.5)

#fprob8 <- function(a,b)
#{
  
#InnerFunc <- function(x, y)
# {
#fG=1/(b^a*fgam(a))*x^(a-1)*(exp(1))^(-x/b)
#fB= 1/(fbet(a,b))*y^(a-1)*(1-y)^(b-1)
# #densitatea comuna f(x,y)
#(fG(x)* fB(y))
#} 
# InnerIntegral= Vectorize( function(x) { integral(InnerFunc,x-0.5,1)$value})
#return(integral(InnerIntegral,0.5,Inf))
#}

#fprob8(2,1)
#fprob8(3,2)
#fprob8(3,1)



####9) P(X+Y>3|X-Y>0.5)

#fprob9 <- function(a,b)
#{
  
# InnerFunc <- function(x, y)
# {
    
#   fG=1/(b^a*fgam(a))*x^(a-1)*(exp(1))^(-x/b)
# fB= 1/(fbet(a,b))*y^(a-1)*(1-y)^(b-1)
#densitatea comuna f(x,y)
# (fG(x)* fB(y))
# } 
  
# ymin <- function(x)
# {
#   min <- (1,(x-0.5))
# }
  
# ymax <- function(x)
# {
#   max <- max(3-x,0)
# }
  
  
# InnerIntegral= Vectorize( function(x) { integral(InnerFunc,ymin,ymax)$value})
# return(integral(InnerIntegral,2,Inf))
  
#}

#fprob9(2,1)
#fprob9(3,2)
#fprob9(3,1)

#################################################3
# d) Tabel diferente functii implementate si functii predefinite
tabel<- function()
{
 
  x <-  matrix(nrow=9,ncol=2, dimnames = list(1:9, c("punctul c)","punctul d)")))
  
  
    x[1,1]=fprobgamma1(2,1)
    x[2,1]=fprobgamma2(2,1)
    x[3,1]=fprobgamma3(2,1)
    x[4,1]=fprobbeta4(2,1)
    x[5,1]=fprobgamma5(2,1)
    x[6,1]=fprobgamma6(2,1)
    #x[7,1]=fprob7(2,1)
    #x[8,1]=fprob8(2,1)
    #x[9,1]=fprob9(2,1)
    x[7,1]=0
    x[8,1]=0
    x[9,1]=0
  
    x[1,2]=pgamma(3,2,1)
    x[2,2]=pgamma(5,2,1)-pgamma(3,2,1)
    x[3,2]=(pgamma(4,2,1)-pgamma(3,2,1))/(1-pgamma(2,2,1))
    x[4,2]=1-pbeta(2,2,1)
    x[5,2]=pgamma(6,2,1)-pgamma(4,2,1)
    x[6,2]=(pgamma(1,2,1)-pgamma(0,2,1))/(pgamma(7,2,1))
    x[7,2]=0
    x[8,2]=0
    x[9,2]=0

error <-  0.5
  return(list(x = as.data.frame(x), error = error))
}
tabel()

########################################################################
#Problema 3

# a)

frepcomgen<- function(n,m)
{
  
  x <-  matrix(nrow=n+1,ncol=m+1)
  
  a <- c("x1")
  
  for(k in 2:n )
    a <- c(a, paste("x",k,sep=""))
  a <- c(a,"pi")
  
  b <- c("y1")
  for(k in 2:m)
    b <- c(b, paste("y",k, sep=""))
  b <- c(b,"qi")
  colnames(x, do.NULL = FALSE)
  colnames(x) <- b
  rownames(x, do.NULL = FALSE)
  rownames(x) <- a
  
  x[,]=0

 q <- 0
 p <- 0
 
 for(i in 1:(n-1))
 {
   r <- runif(1,0.00001, 1-q)
   x[i,m+1] <- r
   q <- q+r
 }
 
 x[n,m+1] <-1-q
 
 for(j in 1:(m-1))
 {
   r <- runif(1, 0.00001, 1-p)
   x[n+1,j] <- r
   p <- p+r
 }
 x[n+1,m]=1-p
 
 s <- rep(0,n)
 t <- rep(0,m)
 for(i in 1:(n-1))
   {
   for(j in 1:(m-1))
   {
    
       if(x[i,m+1]-s[i] < x[n+1,j]-t[j])
       r <- runif(1, 0.00001, x[i,m+1]-s[i])
       else
        r <- runif(1,0.00001,x[n+1,j]-t[j])
       s[i] <- s[i]+r
       t[j] <- t[j]+r
     x[i,j] <- r
   }
   
 }
 for(i in 1:n)
 { 
   x[i,m] <- x[i,m+1]-s[i]
   
 }
 
 for(j in 1:m)
   x[n,j] <- x[n+1,j]-t[j]
 
 min1 <- min(n,m)
 for(k in 1:min1)
   x[k,k] <- -1
 x[n+1,m+1] <- 1
 
  error <-  1
  list(x = as.data.frame(x), error = error)
  return (x)
  
}
frepcomgen(3,4)

########################################################################

# b)

fcomplrepcom <- function(n,m)
{
  x <- frepcomgen(n,m)
  
  
  t <- rep(0,m)
   min1 <- min(m,n)
  for(i in 1:n)
  {
    for(j in 1:min1)
    {
      if(x[i,j]!=-1)
      t[j] <- t[j]+x[i,j]
    }
    
  }
 
 
  for(j in 1:min1)
    x[j,j] <- x[n+1,j]-t[j]
   
   error=-3
  list(x = as.data.frame(x), error = error)
  return (x)
    
}
fcomplrepcom(3,4)

##################################################################

frepcomgen_verif<- function(n,m)
{
  
  x <-  matrix(nrow=n+1,ncol=m+1)
  
  a <- c("x1")
  
  for(k in 2:n )
    a <- c(a, paste("x",k,sep=""))
  a <- c(a,"pi")
  
  b <- c("y1")
  for(k in 2:m)
    b <- c(b, paste("y",k, sep=""))
  b <- c(b,"qi")
  colnames(x, do.NULL = FALSE)
  colnames(x) <- b
  rownames(x, do.NULL = FALSE)
  rownames(x) <- a
  
  x[,]=1/(m*n)
  
  x[n+1,]=1/m
  x[,m+1]=1/n
  x[n+1,m+1] <- 1
  
  min1 <- min(n,m)
  for(k in 1:min1)
    x[k,k] <- -1
  
  error <-  1
  list(x = as.data.frame(x), error = error)
  return (x)
  
}
frepcomgen_verif(4,6)

fcomplrepcom_verif <- function(n,m)
{
  x <- frepcomgen_verif(n,m)
  
  
  t <- rep(0,m)
  min1 <- min(m,n)
  for(i in 1:n)
  {
    for(j in 1:min1)
    {
      if(x[i,j]!=-1)
        t[j] <- t[j]+x[i,j]
    }
    
  }
  
  
  for(j in 1:min1)
    x[j,j] <- x[n+1,j]-t[j]
  
  error <- -3
  list(x = as.data.frame(x), error = error)
  return (x)
  
}
fcomplrepcom_verif(4,6)

##################################################

#Covarianta Cov(3X , 4Y)

cov <- function(n,m)
{
  xy <- fcomplrepcom_verif(n,m)
  x <- c(-2)
  for(k in 2:n)
    x <- c(x, (-3+k))
  
  y <- c(1)
  for(k in 2:m)
    y <- c(y, k)
  
  media_x <- 0
  for(i in 1:n)
  media_x <-media_x + x[i]*xy[i,m+1]

  media_y <- 0
  for(j in 1:m)
    media_y <-media_y + y[j]*xy[n+1,j]

  media_xy <- 0
  for(i in 1:n)
    for(j in 1:m)
      media_xy <- media_xy + x[i]*y[j]*xy[i,j]
 
  covarianta <- 12*media_xy - 12*media_x*media_y
  return (covarianta)
}

cov(2,3)

########################################################

#Probabilitatea P( 0<X<5 | Y>4 )

prob1 <- function(n,m)
{
  xy <- fcomplrepcom_verif(n,m)
  x <- c(-2)
  for(k in 2:n)
    x <- c(x, (-3+k))
  
  y <- c(1)
  for(k in 2:m)
    y <- c(y, k)
  
  prob_numarator <- 0
  prob_numitor <- 0
  for(i in 1:n)
    for(j in 1:m)
    {
      if(x[i]>0 && x[i]<5 && y[j]>4)
        prob_numarator <- prob_numarator+xy[i,j]
      if(y[j]>4)
        prob_numitor <- prob_numitor+xy[i,j]
    }
  prob <- prob_numarator/prob_numitor
  return (prob)
}

prob1(7,5)

########################################################

#Probabilitatea P(X>3, Y<7)

prob2 <- function(n,m)
{
  xy <- fcomplrepcom_verif(n,m)
  x <- c(-2)
  for(k in 2:n)
    x <- c(x, (-3+k))
  
  y <- c(1)
  for(k in 2:m)
    y <- c(y, k)
  
  prob <- 0
  for(i in 1:n)
    for(j in 1:m)
    {
      if(x[i]>3 && y[j]<7)
        prob <- prob+xy[i,j]
    }
  return (prob)
}

prob2(7,5)

###########################################################
#Verificare variabile aleatoare necorelate

fvernecor <- function(n,m)
{
  x <- cov(n,m)
  if(x==0)
    mesaj <- "Variabilele aleatoare X si Y sunt necorelate"
  else 
    mesaj <- "Variabilele aleatoare X si Y nu sunt necorelate"
  return (mesaj)
}
fvernecor(2,3)

##########################################################
#Verificare variabile aleatoare independente

fverind <- function(n,m)
{
  xy <- fcomplrepcom_verif(n,m)
  
  verif <- 1
  
  for(i in 1:n)
    for(j in 1:m)
    {
      if(xy[i,j]!=xy[i,m+1]*xy[n+1,j])
        verif <- 0
    }
  
  if(verif==1)
    mesaj <- "Variabilele aleatoare X si Y sunt independente"
  else 
    mesaj <- "Variabilele aleatoare X si Y nu sunt independente"
  return (mesaj)
}

fverind(2,3)