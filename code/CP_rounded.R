lasso = function(XX,YY,lam){ # minimize 1/2n ||y-X*b||^2_2 + lam * ||b||_1
glmnet(XX,YY,standardize=FALSE,intercept=FALSE,lambda=lam)$beta
}

PI_oracle = function(X,Y,xnew,EY_X,alpha,sigma){
	c(EY_X(xnew) - sigma*qnorm(1-alpha/2),EY_X(xnew) + sigma*qnorm(1-alpha/2))
}

PI_parametric = function(X,Y,xnew,lam,alpha,sigma){
	betahat = lasso(X,Y,lam); betahat[abs(betahat)<=1e-8] = 0
	S = which(betahat!=0)
	SE = sqrt(1 + t(xnew[S])%*%solve(t(X[,S])%*%X[,S],xnew[S]))
	c(sum(xnew*betahat) - sigma*qnorm(1-alpha/2)*SE,sum(xnew*betahat) + sigma*qnorm(1-alpha/2)*SE)
}

round_Y = function(YY,grid){
	round_Y_one = function(YYone){
		grid[min(which.min(abs(YYone-grid)))]
	}
	unlist(lapply(YY,round_Y_one))
}

PI_CPDD = function(X,Y,xnew,lam,alpha,sigma,grid){
	n = length(Y)
	M = length(grid)
	Del = grid[2]-grid[1] # assuming an evenly spaced grid
	Y_r = round_Y(Y,grid)
	PI = c(Inf,-Inf)
	for(m in 1:M){
		betahat = lasso(rbind(X,xnew),c(Y_r,grid[m]),lam)
		resids_train = Y_r - X%*%betahat
		resids_Q = sort(abs(resids_train))[ceiling((1-alpha)*(n+1))]
		if(abs(grid[m]-sum(xnew*betahat)) <= resids_Q){
			PI[1] = min(PI[1], grid[m]-Del/2)
			PI[2] = max(PI[2], grid[m]+Del/2)
		}
	}
	PI
}

PI_CPDM = function(X,Y,xnew,lam,alpha,sigma,grid){
	n = length(Y)
	M = length(grid)
	Del = grid[2]-grid[1] # assuming an evenly spaced grid
	Y_r = round_Y(Y,grid)
	PI = c(Inf,-Inf)
	for(m in 1:M){
		betahat = lasso(rbind(X,xnew),c(Y_r,grid[m]),lam)
		resids_train = Y - X%*%betahat
		resids_Q = sort(abs(resids_train))[ceiling((1-alpha)*(n+1))]
		ymin = max(grid[m]-Del/2, sum(xnew*betahat)-resids_Q)
		ymax = min(grid[m]+Del/2, sum(xnew*betahat)+resids_Q)
		if(ymin <= ymax){
			PI[1] = min(PI[1], ymin)
			PI[2] = max(PI[2], ymax)
		}
	}
	PI
}

PI_approxCP = function(X,Y,xnew,lam,alpha,sigma,grid){
	n = length(Y)
	M = length(grid)
	Del = grid[2]-grid[1] # assuming an evenly spaced grid
	PI = c(Inf,-Inf)
	for(m in 1:M){
		betahat = lasso(rbind(X,xnew),c(Y,grid[m]),lam)
		resids_train = Y - X%*%betahat
		resids_Q = sort(abs(resids_train))[ceiling((1-alpha)*(n+1))]
		if(abs(grid[m]-sum(xnew*betahat)) <= resids_Q){
			PI[1] = min(PI[1], grid[m]-Del/2)
			PI[2] = max(PI[2], grid[m]+Del/2)
		}
	}
	PI
}


nlist = c(100,400); p = 200; sigma = 1; k = 10
EY_X = function(x){(sum(x[1:k]) + sum(sign(x[1:k])*sqrt(abs(x[1:k]))))/sqrt(k)}
alpha = 0.1
ntrial = 1000
Mlist = 5*2^(0:5)

library(glmnet)
PIlength = PIcov = array(0,c(5,length(Mlist),length(nlist),ntrial))

for(i in 1:ntrial){
	print(i)
	for(i_n in 1:length(nlist)){
	n = nlist[i_n]
	lam = sigma*sqrt(log(p)/2/n)
	
	set.seed(100 + i)

	X = matrix(rnorm(n*p),n,p); xnew = rnorm(p)
	Y = apply(X,1,EY_X) + sigma*rnorm(n)
	ynew = EY_X(xnew) + sigma*rnorm(1)

	PI = PI_oracle(X,Y,xnew,EY_X,alpha,sigma)
	PIlength[1,,i_n,i] = max(0,PI[2]-PI[1])
	if(PI[1]<=ynew & ynew<=PI[2]){PIcov[1,,i_n,i] = 1}
	
	PI = PI_parametric(X,Y,xnew,lam,alpha,sigma)
	PIlength[2,,i_n,i] = max(0,PI[2]-PI[1])
	if(PI[1]<=ynew & ynew<=PI[2]){PIcov[2,,i_n,i] = 1}

	for(iM in 1:length(Mlist)){
		M = Mlist[iM]
		grid = min(Y) + (max(Y)-min(Y)) * (0.5 + (0:(M-1)))/M

		PI = PI_approxCP(X,Y,xnew,lam,alpha,sigma,grid)
		PIlength[3,iM,i_n,i] = max(0,PI[2]-PI[1])
		if(PI[1]<=ynew & ynew<=PI[2]){PIcov[3,iM,i_n,i] = 1}
		
		PI = PI_CPDD(X,Y,xnew,lam,alpha,sigma,grid)
		PIlength[4,iM,i_n,i] = max(0,PI[2]-PI[1])
		if(PI[1]<=ynew & ynew<=PI[2]){PIcov[4,iM,i_n,i] = 1}

		PI = PI_CPDM(X,Y,xnew,lam,alpha,sigma,grid)
		PIlength[5,iM,i_n,i] = max(0,PI[2]-PI[1])
		if(PI[1]<=ynew & ynew<=PI[2]){PIcov[5,iM,i_n,i] = 1}

	}
}}

PIlength_means = PIcov_means = PIlength_SEs = PIcov_SEs = 
	array(0,c(5,length(Mlist),length(nlist)))
for(i in 1:5){for(iM in 1:length(Mlist)){for(i_n in 1:length(nlist)){
		PIlength_means[i,iM,i_n] = mean(PIlength[i,iM,i_n,])
		PIlength_SEs[i,iM,i_n] = sd(PIlength[i,iM,i_n,])/sqrt(ntrial)
		PIcov_means[i,iM,i_n] = 100*mean(PIcov[i,iM,i_n,])
		PIcov_SEs[i,iM,i_n] = 100*sd(PIcov[i,iM,i_n,])/sqrt(ntrial)
}}}


cols = c('dodgerblue3','goldenrod3','darkorchid2','chartreuse4','firebrick3')
pchs = c(25,22,23,24,21)
ltys = rep(1,5)
names = c('Oracle','Parametric','Approximate CP','CP w/ discr. data','CP w/ discr. model')

save('PIlength','PIcov',file='tmpresults.RData')

for(i_n in 1:length(nlist)){

pdf(paste0('CP_sims_plot_',i_n,'.pdf'),8.5,4.5)
par(mfrow=c(1,2))


plot(1,1,type='n',xlim=range(Mlist),log='x',ylim=c(0.9*min(PIlength_means[,,i_n]),1.1*max(PIlength_means[,,i_n])),
	axes=FALSE,xlab='# grid points M',ylab='Average PI length')
axis(side=1,at=Mlist)
axis(side=2,las=1)
for(i in 1:5){
	points(Mlist,PIlength_means[i,,i_n],col=cols[i],bg=cols[i],pch=pchs[i],lty=ltys[i],type='l')
	points(Mlist,PIlength_means[i,,i_n],col=cols[i],bg=cols[i],pch=pchs[i],lty=ltys[i])
	for(iM in 1:length(Mlist)){
		segments(Mlist[iM],PIlength_means[i,iM,i_n]-PIlength_SEs[i,iM,i_n],Mlist[iM],PIlength_means[i,iM,i_n]+PIlength_SEs[i,iM,i_n],col=cols[i])
	}
}
# legend('topright',legend = names, col=cols, pt.bg=cols, pch=pchs, lty=ltys)

plot(1,1,type='n',xlim=range(Mlist),log='x',ylim=100*c(1-2*alpha,1),
	axes=FALSE,xlab='# grid points M',ylab='Empirical coverage (%)')
axis(side=1,at=Mlist)
axis(side=2,las=1)
for(i in 1:5){
	points(Mlist,PIcov_means[i,,i_n],col=cols[i],bg=cols[i],pch=pchs[i],lty=ltys[i],type='l')
	points(Mlist,PIcov_means[i,,i_n],col=cols[i],bg=cols[i],pch=pchs[i],lty=ltys[i])
	for(iM in 1:length(Mlist)){
		segments(Mlist[iM],PIcov_means[i,iM,i_n]-PIcov_SEs[i,iM,i_n],Mlist[iM],PIcov_means[i,iM,i_n]+PIcov_SEs[i,iM,i_n],col=cols[i])
	}
}
legend('topright',legend = names, col=cols, pt.bg=cols, pch=pchs, lty=ltys)


dev.off()
}

