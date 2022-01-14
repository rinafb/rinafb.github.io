set.seed(1)

n = 100
r = 5
psamp = 0.2

neta = 20
eta_min = 0.5; eta_max = 4
eta_list = seq(from=eta_min,to=eta_max,length.out=neta)
SNR = 5


niter = 50 ; niter_init = 1
ntrial = 50

loss_store = array(0,c(3,neta,niter+1,2,ntrial,2))
# methods 1,2,3 ; stepsize ; iteration ; noiseless or noisy ; trial ; initialize


for(trial in 1:ntrial){
	print(trial)



utrue = svd(matrix(rnorm(n*r),n,r))$u
xtrue = utrue%*%t(utrue) # all singular values = 1

sig2 = mean(xtrue^2)/SNR

Omega = matrix(rbinom(n*n,1,psamp),n,n); Omega = Omega*upper.tri(Omega); Omega = Omega+t(Omega)
Noise = sqrt(sig2)*matrix(rnorm(n^2),n,n); Noise = Noise*upper.tri(Noise); Noise = Noise+t(Noise)
Obs_noiseless = xtrue*Omega
Obs_noisy = (xtrue + Noise)*Omega




loss = function(umat,obs){1/2*sum(((umat%*%t(umat))*Omega-obs)^2)}



update = function(umat,method,eta,obs){
	grad = (umat%*%t(umat)-obs)*Omega
	if(method==1){
		x = umat%*%t(umat) - eta*grad
		svd_x = svd(x)
		umat = svd_x$u[,1:r]%*%diag(sqrt(svd_x$d[1:r]))
	}
	if(method==2){
		basis = svd(umat)$u
		proj_grad = basis%*%t(basis)%*%grad + grad%*%basis%*%t(basis) - basis%*%t(basis)%*%grad%*%basis%*%t(basis)
		x = umat%*%t(umat) - eta*proj_grad
		svd_x = svd(x)
		umat = svd_x$u[,1:r]%*%diag(sqrt(svd_x$d[1:r]))
	}
	if(method==3){
		umat = umat - eta*grad%*%umat
	}
	umat
}


for(noise in 1:2){for(initialize in 1:2){
	if(noise==1){obs = Obs_noiseless}else{obs = Obs_noisy}
	svd_init = svd(obs)
	u0 = svd_init$u[,1:r]%*%diag(sqrt(svd_init$d[1:r]))
	loss_store[,,1,noise,trial,initialize] = loss(u0,obs)
	for(ieta in 1:neta){
		eta = eta_list[ieta]
		for(method in 1:3){
			u = u0
			for(t in 1:niter){
				if(initialize==2 & t<=niter_init){
					u = update(u,1,eta,obs)
				}else{
					u = update(u,method,eta,obs)
				}
				loss_store[method,ieta,t+1,noise,trial,initialize] = loss(u,obs)
			}
		}
	}
}}
		

}


save('loss_store',file='lowrank_sim_results.RData')




# plot results



rowQ1 = function(mat){Q1=function(vec){quantile(vec,0.25)};apply(mat,1,Q1)}
rowQ3 = function(mat){Q3=function(vec){quantile(vec,0.75)};apply(mat,1,Q3)}
rowMed = function(mat){Med=function(vec){quantile(vec,0.5)};apply(mat,1,Med)}

plot_band = function(x,bottom,top,rgb){
	nx = length(x)
	polygon(c(x,x[nx:1],x[1]),c(bottom,top[nx:1],bottom[1]),
		col=rgb(rgb[1],rgb[2],rgb[3],0.15),border=NA)
}

ieta_opt = array(0,c(2,2,3)) # noise, initialize, method
ymax = matrix(0,2,2)
for(noise in 1:2){for(initialize in 1:2){
	for(method in 1:3){
		ieta_opt[noise,initialize,method] = max(which.min(rowMeans(
			loss_store[method,,niter+1,noise,,initialize])))
		ymax[noise,initialize]=max(ymax[noise,initialize],
			max(rowQ3(loss_store[method,ieta_opt[noise,initialize,method],,noise,,initialize])))
	}	
}}
titles = matrix(c('Noiseless','Noisy','Noiseless + initialization','Noisy + initialization'),2,2)
rgb = matrix(c(0,0,0,0.25,0.75,0.25,0.25,0.25,0.75),3,3)
for(noise in 1:2){for(initialize in 1:2){
	pdf(paste0('lowrank_sim_results_noise',noise-1,'_initialize',initialize-1,'.pdf'),4,4)
	plot(0:1,0:1,type='n',xlab='Iteration',ylab='Loss',xlim=c(0,niter),ylim=c(0,ymax[noise,initialize]),las=1)
	title(main=titles[noise,initialize])
	for(method in 1:3){
		plot_band(0:niter,
			rowQ1(loss_store[method,ieta_opt[noise,initialize,method],,noise,,initialize]),
				rowQ3(loss_store[method,ieta_opt[noise,initialize,method],
					,noise,,initialize]),rgb[,method])
	}
	for(method in 1:3){
		points(0:niter,
			rowMed(loss_store[method,ieta_opt[noise,initialize,method],,noise,,initialize]),
				col=rgb(rgb[1,method],rgb[2,method],rgb[3,method]),type='l')
	}
	legend('topright',c('Projected grad. desc.','Approx. proj. grad. desc.',
		'Factored grad. desc.'),
			fill=c(rgb(rgb[1,1],rgb[2,1],rgb[3,1]),rgb(rgb[1,2],rgb[2,2],rgb[3,2]),
				rgb(rgb[1,3],rgb[2,3],rgb[3,3])))
	dev.off()
}}


