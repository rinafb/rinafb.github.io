
####################################################################
#### Settings ######################################################
####################################################################

signal_amp = 10
sigma_noise = 1
delta = 0.1 # confidence level

# function f gives the true signal x
f = function(t){pmin(signal_amp,pmax(-signal_amp,(t-0.5)*5*signal_amp))}

# number of sampled points
nlist = 700:1000

####################################################################
#### PAVA algorithm ################################################
####################################################################

# fit min{||y-x||_2 : x1 <= x2 <= ... <= xn}
PAVA = function(y){
n = length(y); x = y
blockstarts = 1:n; blockends = 1:n; nblock = n; done = FALSE
while(!done){
	done = TRUE; i=1
	while(i<nblock){
		this_block = blockstarts[i]:blockends[i]; next_block = blockstarts[i+1]:blockends[i+1]
		if(mean(y[this_block])>mean(y[next_block])){ # pool these blocks, then test block i again
			x[c(this_block,next_block)] = mean(y[c(this_block,next_block)])
			blockends[i] = blockends[i+1]
			if(i+1<nblock){ # shift over remaining blocks
				blockstarts[(i+1):(nblock-1)] = blockstarts[(i+2):nblock]
				blockends[(i+1):(nblock-1)] = blockends[(i+2):nblock]}
			nblock = nblock - 1; done = FALSE
		}else{i = i+1}}}
x
}


####################################################################
#### Some setup ####################################################
####################################################################

set.seed(1)

results = NULL

nn = length(nlist)
store_x = store_y = store_xh = store_band1 = store_band2 = matrix(0,max(nlist),nn)


####################################################################
#### Run experiment at each sample size ############################
####################################################################


for(ni in 1:nn){
n = nlist[ni]

# generate the signal x & noisy observations y
x = f((1:n)/(n+1))
y = x + sigma_noise*rnorm(n)

# isotonic regression
xh = PAVA(y)

# bound on the sliding window norm ||x-y||_SW, with probability 1-delta
SW = sigma_noise*sqrt(2*log((n^2+n)/delta))

# isotonic projection is contractive with respect to sliding window norm,
# so ||x - xh||_SW <= ||x-y||_SW
# we use this to calculate a confidence band: for i<=k and k<=j, 
# x[k] <= mean(x[k:j]) <= mean(xh[k:j]) + sw(x-y)/sqrt(j-k+1)
# x[k] >= mean(x[i:k]) >= mean(xh[i:k]) - sw(x-y)/sqrt(k-i+1)

# calculate the confidence band by taking min over j & max over i
band = c(-Inf,Inf)%*%t(rep(1,n))
for(k in 1:n){
	for(i in 1:k){
		band[1,k] = pmax(band[1,k],mean(xh[i:k])-SW/sqrt(k-i+1))
	}
	for(j in k:n){
		band[2,k] = pmin(band[2,k],mean(xh[k:j])+SW/sqrt(j-k+1))
	}
}

# store all data & results for this sample size
store_x[1:n,ni] = x
store_y[1:n,ni] = y
store_xh[1:n,ni] = xh
store_band1[1:n,ni] = band[1,]
store_band2[1:n,ni] = band[2,]

}



####################################################################
#### Save results to file ##########################################
####################################################################

save(store_x,store_y,store_xh,store_band1,store_band2,file='iso_sim_storeresults.RData')

####################################################################
#### Plot original signal ##########################################
####################################################################

pdf('iso_signal.pdf',5,5.5)
n = nlist[nn]
plot(0:1,0:1,type='n',ylim=c(-15,15),main='',xlab='t',ylab='f(t)',las=1)
tval=(0:10000)/10000
points(tval,f(tval),type='l',lwd=2,lty='dashed')
rect(0.1,-11,0.2,-9,col=rgb(0.65,0.9,0.9,0.5),border=NA)
rect(0.8,9,0.9,11,col=rgb(0.65,0.9,0.9,0.5),border=NA)
rect(0.4,-5,0.6,5,col=rgb(0.85,0.75,1,0.5),border=NA)
text(0.85,-10.5,'Flat regions')
text(0.45,8,'Increasing region')
arrows(0.85,-9.8,0.85,8.5,length=0.05)
arrows(0.7,-10.5,0.22,-10.5,length=0.05)
arrows(0.45,7,0.45,5.5,length=0.05)
legend('topleft',legend=expression(paste('Signal ',x[i])),lwd=2)
dev.off()

####################################################################
#### Plot confidence band at one particular sample size n ##########
####################################################################

pdf('iso_confidenceband.pdf',5,5.5)
n = nlist[nn]
plot(1:n,store_y[,nn],col='dark turquoise',pch=20,cex=.5,ylim=c(-15,15),main='',xlab='Index i',ylab='',las=1)
points(1:n,store_xh[,nn],type='l',lwd=2)
points(1:n,store_band1[,nn],type='l',col='firebrick',lwd=2)
points(1:n,store_band2[,nn],type='l',col='firebrick',lwd=2)
legend('topleft',legend=c(expression(paste('Data ',y[i])),expression(paste('Estimate ',iso(y)[i])),'Confidence band'),lwd=c(NA,2,2),pch=c(20,NA,NA),col=c('dark turquoise','black','firebrick'))
dev.off()




####################################################################
#### Log-log plots for (n/log(n))^{-1/2} or (..)^{-1/3} scaling ####
####################################################################


band_width = store_band2 - store_band1
means_flat = means_incr = rep(0,nn)
for(ni in 1:nn){
	n = nlist[ni]
	tval = (1:n)/(n+1)
	inds1 = which((0.1<=tval & tval<=0.2)|(0.8<=tval & tval<=0.9))
	inds2 = which(0.4<=tval & tval<=0.6)
	means_flat[ni] = mean(band_width[inds1,ni])
	means_incr[ni] = mean(band_width[inds2,ni])
}

pdf('iso_loglogplot_flat.pdf',5,5.5)
plot(log(nlist/log(nlist)),log(means_flat),xlab='log( n/log(n) )',ylab='log( mean confidence band width )',main='',las=1,pch=20,cex=0.8)
abline(lm(log(means_flat)~log(nlist/log(nlist))),col='red',lty='dashed',lwd=2)
legend('topright',legend=c('Least squares regression line',(paste('(Slope = ',round(lm(log(means_flat)~log(nlist/log(nlist)))$coef[2],4),')'))),col='red',lty='dashed',lwd=c(2,NA))
dev.off()

pdf('iso_loglogplot_incr.pdf',5,5.5)
plot(log(nlist/log(nlist)),log(means_incr),xlab='log( n/log(n) )',ylab='log( mean 5,6confidence band width )',main='',las=1,pch=20,cex=0.8)
abline(lm(log(means_incr)~log(nlist/log(nlist))),col='red',lty='dashed',lwd=2)
legend('topright',legend=c('Least squares regression line',(paste('(Slope = ',round(lm(log(means_incr)~log(nlist/log(nlist)))$coef[2],4),')'))),col='red',lty='dashed',lwd=c(2,NA))
dev.off()
