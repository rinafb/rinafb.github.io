


###################################
# simulation
###################################

set.seed(12345)

n_list = seq(from=10000,to=20000,by=2000)
n_trial = 50#500000
noise_sd = 0.1

n_n = length(n_list)
bias_store_curve = bias_store_hinge = rep(0,n_n)

f_curve = function(t){t+sin(4*pi*t)/16}
f_hinge = function(t){pmax(0.1*t,1.9*t-0.9)}

for(i_n in 1:n_n){
	n = n_list[i_n]
	mu_curve = f_curve((1:n)/n)
	mu_hinge = f_hinge((1:n)/n)
	EisoY_curve = EisoY_hinge = rep(0,n)
	for(i_trial in 1:n_trial){
		if(i_trial %% 10000 == 0){
			print(paste0('n=',n,',trial #',i_trial))
		}
		EisoY_curve = EisoY_curve + isoreg(mu_curve+noise_sd*rnorm(n))$yf/n_trial
		EisoY_hinge = EisoY_hinge + isoreg(mu_hinge+noise_sd*rnorm(n))$yf/n_trial
	}
	midpt = round(n/2)
	inds = round(n*seq(from=0.1,to=0.9,by=0.05))
	bias_store_curve[i_n] = mean(abs(EisoY_curve[inds] - mu_curve[inds]))
	bias_store_hinge[i_n] = abs(EisoY_hinge[midpt] - mu_hinge[midpt])
}

results = cbind(n_list,bias_store_curve,bias_store_hinge)


par(mfrow=c(2,2))
n=10000;midpt = round(n/2);inds = round(n*seq(from=0.1,to=0.9,by=0.05))
mu_curve = f_curve((1:n)/n);mu_hinge = f_hinge((1:n)/n)
plot(1:n,mu_curve,type='l',xlab='Index j=1,...,n',ylab='',main=expression(paste('Smooth signal (',beta==2,')')))
points(inds,mu_curve[inds],pch=20)
plot(1:n,mu_hinge,type='l',xlab='Index j=1,...,n',ylab='',main=expression(paste('Non-smooth signal (',beta==1,')')))
points(midpt,mu_hinge[midpt],pch=20)




lm_curve = lm(log(bias_store_curve)~log(n_list))
lm_hinge = lm(log(bias_store_hinge)~log(n_list))


plot(log(n_list),log(bias_store_curve),xlab='n',ylab='bias',axes=FALSE,pch=20,main=expression(paste('Bias (',beta==2,')')))
axis(side=1,at=log(n_list),lab=n_list)
axis(side=2,at=axTicks(2),lab=formatC(exp(axTicks(2)),format='e',digits=2))
abline(lm_curve,col='red')
text(quantile(log(n_list),.75),mean(lm_curve$fit),paste('Slope = ',round(lm_curve$coef[2],3)),col='red')

plot(log(n_list),log(bias_store_hinge),xlab='n',ylab='bias',axes=FALSE,pch=20,main=expression(paste('Bias (',beta==1,')')))
axis(side=1,at=log(n_list),lab=n_list)
axis(side=2,at=axTicks(2),lab=formatC(exp(axTicks(2)),format='e',digits=2))
abline(lm_hinge,col='red')
text(quantile(log(n_list),.75),mean(lm_hinge$fit),paste('Slope = ',round(lm_hinge$coef[2],3)),col='red')

