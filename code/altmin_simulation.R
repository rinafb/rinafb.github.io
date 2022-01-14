set.seed(1)

# parameters

n = 200
d = 50
m = 20
r = 3

SNR = 3
Cor = c(0.0,0.4,0.6)

ntrial = 100

# writing functions

dataGen = function(n,d,m,r,SNR,cor){
  utrue = svd(matrix(rnorm(m*r),m,r))$u
  vtrue = svd(matrix(rnorm(d*r),d,r))$u
  xtrue = utrue%*%t(vtrue) # all singular values = 1
  
  phi_Sigma = (0.3)^(abs(matrix(rep(1:d,d),d,d) - t(matrix(rep(1:d,d),d,d))))
  phi = matrix(rnorm(n*d),n,d)%*%chol(phi_Sigma)
  
  sig2 = mean(((phi%*%t(xtrue))^2))/SNR
  Sigma_true = sig2 * (cor)^(abs(matrix(rep(1:m,m),m,m) - t(matrix(rep(1:m,m),m,m))))
  
  y = phi%*%t(xtrue) + matrix(rnorm(n*m),n,m)%*%chol(Sigma_true)
  
  mydata = list(y=y,phi=phi)
  return(mydata)
}


initialization = function(mydata){
  y = mydata$y; phi = mydata$phi
  
  x = t(solve(t(phi)%*%phi, t(phi)%*%y))
  svd_x = svd(x)
  x = svd_x$u[,1:r]%*%diag(svd_x$d[1:r])%*%t(svd_x$v[,1:r]) 
  Sigma = t(y - phi%*%t(x))%*%(y - phi%*%t(x))/n
  
  vars = list(x=x,Sigma=Sigma)
  return(vars)
}

loss = function(vars,mydata){
  x = vars$x; Sigma = vars$Sigma
  y = mydata$y; phi = mydata$phi
  
  emp_cov = t(y - phi%*%t(x))%*%(y - phi%*%t(x))/n
  f = log(det(Sigma)) + sum(diag(solve(Sigma, emp_cov)))
  
  return(f)
}


update = function(vars,mydata,eta,method,inner1,inner2){
  x = vars$x; Sigma = vars$Sigma
  y = mydata$y; phi = mydata$phi
  etax = eta$x; etaSigma = eta$Sigma
  
  thresh = 1e-2
  
  # alternating method
  if(method==1){
    Sigma = t(y - phi%*%t(x))%*%(y - phi%*%t(x))/n
    for(inn in 1:inner1){
      gradx = 2 * solve(Sigma, (t(phi%*%t(x) - y)%*%phi))/n
      x_sub = x - etax * gradx
      svd_x = svd(x_sub)
      x = svd_x$u[,1:r]%*%diag(svd_x$d[1:r])%*%t(svd_x$v[,1:r])
    }
  }
  
  # joint gradient descent
  if(method==2){
    for(inn in 1:inner2){
      gradx = 2 * solve(Sigma, (t(phi%*%t(x) - y)%*%phi))/n
      gradSigma = -Sigma + t(y - phi%*%t(x))%*%(y - phi%*%t(x))/n
      
      Theta = solve(Sigma) - etaSigma * gradSigma
      eig_Theta = eigen(Theta)
      if(min(eig_Theta$values) <= thresh){
        # diverge
        return(vars)
      }
      Sigma = solve(Theta)
      x_sub = x - etax * gradx
      svd_x = svd(x_sub)
      x = svd_x$u[,1:r]%*%diag(svd_x$d[1:r])%*%t(svd_x$v[,1:r])
    }
  }
  return(list(x=x,Sigma=Sigma))
}


# simulation

maxIter = 1200
inner1 = 1; inner2 = 1

etax = 0.001
neta = 30
etaSigma = exp(seq(from=log(5),to=log(400),length.out=neta))

lossErr = array(0,c(2,length(Cor),maxIter+1,ntrial)) # method, cor, max_iteration, trial
lossMin = matrix(0,maxIter+1,neta) # max_iteration, step size

for(trial in 1:ntrial){
  
  print(trial)
  
  for(cor in 1:length(Cor)){
    mydata = dataGen(n,d,m,r,SNR,Cor[cor])
    vars_ini = initialization(mydata)
    
    lossErr[,cor,1,trial] = lossMin[1,] = loss(vars_ini,mydata)
    
    for(ieta in 1:neta){
      
      eta = list(x=etax, Sigma=etaSigma[ieta])
      vars_altmtd = vars_jointgrad = vars_ini
      for(t in 1:maxIter){
        
        if(ieta == 1){
          vars_altmtd = update(vars_altmtd,mydata,eta,method=1,inner1,inner2)
          lossErr[1,cor,t+1,trial] = loss(vars_altmtd,mydata)
        }
        
        vars_jointgrad = update(vars_jointgrad,mydata,eta,method=2,inner1,inner2)
        lossMin[t+1,ieta] = loss(vars_jointgrad,mydata)
        
      }
    }
    ieta_opt = max(which.min(apply(lossMin,2,min)))
    #ieta_opt = max(which.min(apply(lossMin[(maxIter-8):(maxIter+1),],2,sum)))
    lossErr[2,cor,,trial] = lossMin[,ieta_opt]
  }
}

save('lossErr',file='multitask_results.RData')


# plot results

rowQuan = function(mat){
  Q1 = apply(mat,1,function(vec){quantile(vec,0.25)})
  Med = apply(mat,1,function(vec){quantile(vec,0.5)})
  Q3 = apply(mat,1,function(vec){quantile(vec,0.75)})
  return(list(Q1=Q1,Med=Med,Q3=Q3))
}

make_plots = function(s,y1,y2,rgb,rho,name){
  
    pdf(name,4,4)
    
    xmax = length(s)-1 
    ymax = max(log(y1$Q3[1]),log(y2$Q3[1]))
    ymin = min(log(y1$Q1[xmax+1]),log(y2$Q1[xmax+1]))
    
    plot(0:1,0:1,type='n',xlab='Iteration',ylab='Log(excess loss)',xlim=c(0,xmax),ylim=c(ymin,ymax),las=0)
    title(main=bquote(paste('Convergence for ', rho==.(rho))))
    
    polygon(c(s,rev(s)),c(log(y1$Q1[s+1]),rev(log(y1$Q3[s+1]))),col=rgb(rgb[1,1],rgb[2,1],rgb[3,1],0.25),border=NA)
    polygon(c(s,rev(s)),c(log(y2$Q1[s+1]),rev(log(y2$Q3[s+1]))),col=rgb(rgb[1,2],rgb[2,2],rgb[3,2],0.25),border=NA)
    points(s,log(y1$Med[s+1]),col=rgb(rgb[1,1],rgb[2,1],rgb[3,1]),type='l')
    points(s,log(y2$Med[s+1]),col=rgb(rgb[1,2],rgb[2,2],rgb[3,2]),type='l')
    
    legend('bottomleft',c('Alt_method','Joint_grad'),fill=c(rgb(rgb[1,1],rgb[2,1],rgb[3,1],0.25),rgb(rgb[1,2],rgb[2,2],rgb[3,2],0.25)))
    
    
    dev.off()
}

for(cor in 1:length(Cor)){
	losshat = rep(0,ntrial)
	for(trial in 1:ntrial){ losshat[trial] = min(lossErr[,cor,,trial]) }
  for(method in 1:2){
    lossDiff = t(apply(lossErr[method,cor,,],1,function(vec){vec-losshat}))
    name = paste('Diff', cor, method, sep = '')
    assign(name, rowQuan(lossDiff))
  }
}


s = 0:600  
rgb = matrix(c(0.25,0.75,0,0.25,0,0.75),3,2)

make_plots(s,Diff11,Diff12,rgb,Cor[1],'multitask_logplot1.pdf')
make_plots(s,Diff21,Diff22,rgb,Cor[2],'multitask_logplot2.pdf')
make_plots(s,Diff31,Diff32,rgb,Cor[3],'multitask_logplot3.pdf')



