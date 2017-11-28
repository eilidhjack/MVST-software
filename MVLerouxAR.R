##Block update for Metropolis Hastings Poisson model##
##Update in blocks of beta parameters (block size specifified by user)##

space.time<-function(formula, samples, burnin, prior.mean.beta=NULL, prior.var.beta=NULL,  
                     progress.bar=T, thin=1, block=5, d,W=W, tp=1, S=14, X.HB, X.phi){
    ##Arguments
    ##formula - Insert formula for the covariate part of the model
    ##samples - number of samples to be drawn
    ##burnin - number of samples to be removed from the chain
    ##prior.mean.beta - a vector for the prior mean of each beta (defaults to a vector of 0's)
    ##prior.var.beta - a vector for the prior variances for each beta (defaults to a vector of 100's)
    ##progress.bar - whether progress bar appears or not, default = T
    ##thin - value to thin by, default = 1
    ##block - number of parameters per block, default=5
    ##d - number of diseases in data
    ##W - W neighborhood matrix
    ##tp - number of time points
    ##S - number of Health Boards
    ##X.HB - Matrix of HB by year
    ##X.phi - Design matrix of factors for HB to allow for phi to be mean centered by HB
    #################################################################################################
    ##Initial calculations
    #################################################################################################
    ##Create progress bar
    if(progress.bar){
        pb<-txtProgressBar(min=0, max=samples,initial=0, char="=", style=3)}
    t0<-Sys.time()
    ## Create design matrix from formula argument
    data<-model.frame(formula)
    y<-model.response(data)
    X.big<-model.matrix(formula)
    
    ##Number of observations for each disease
    n <- length(X.big[,1])/d
    
    ## Number of data points
    n.all <- length(y)
    ##Number of IZ's to be sampled within each disease
    n.phi<-n/tp
    
    ##Create y matrix with a column for each disease
    y.mat<-matrix(y,n,d,byrow=F)
    y.arr <- array(y, c(n.phi, tp, d))
    
    ##The design matrix X is the same for all diseases so just hold once
    X<-X.big[1:n,]
    
    ## Format the arguments   
    HB.p <- S*tp  
    p <- dim(X)[ 2 ] - HB.p + 1
    
    ##Store Design Matrix seperately for Beta parameters and HB parameters
    X.beta <- X[,1:p]
    
    
    #### Offset variable
    ## Create the offset
    offset <- try(model.offset(data), silent=TRUE)
    if(class(offset)=="try-error")   stop("the offset is not numeric.", call.=FALSE)
    if(is.null(offset))  offset <- rep(0,n.all)
    
    ##Create offset matrix with a column for each disease
    offset.mat<-matrix(offset,n,d,byrow=F)
    offset.arr<-  array(offset, c(n.phi, tp, d))
    
    #################################################################################################
    ##Standardise Covariates
    #################################################################################################
    ##Code each covariate as factor or continuous
    ##Intercept = 0, Continuous = 1, Factor = 2
    X.orig <- X.beta
    data.means<-rep(apply(X.beta,2,mean),d)
    data.sd<-rep(apply(X.beta,2,sd),d)
    data.type<-rep(0,p)
    for(i in 1:p){
        data.type[i]<-length(unique(X.beta[,i]))
    }
    data.type[1]<-0
    for (i in 2:p){
        if(data.type[i]!=2){
            data.type[i]<-1
        }
    }
    
    ##Now standardise all continuous covariates
    for(i in 1:p){
        if(data.type[i]==1){
            X.beta[,i]<-scale(X.orig[,i])
        }
    }
    data.type<-rep(data.type,d)
    
    
    #################################################################################################
    ##Create triplet object for W
    #################################################################################################
    W.triplet<-c(0,0,0)
    for(i in 1:n.phi){
        for(j in 1:n.phi){
            if (W[i,j] != 0){
                W.triplet<-rbind(W.triplet, c(i, j, W[i,j]))
            }
        }
    }
    W.triplet<-W.triplet[-1,]
    n.neighbours<-apply(W,1,sum)
    triplet.length<-length(W.triplet[,1])
    
    track<-1
    W.start<-rep(0,n.phi)
    W.end<-rep(0,n.phi)
    for(i in 1:n.phi){
        W.start[i]<-track
        W.end[i]<-track + n.neighbours[i] - 1
        track<-W.end[i]+1
        
    }
    
    # W.sparse<-Matrix(W,sparse=T)
    
    ###################################################
    ##Priors##
    ###################################################
    
    ##Prior mean for beta
    if(is.null(prior.mean.beta)){
        prior.mean.beta<-rep (0 , p)}
    
    ##Prior variance for beta
    if(is.null(prior.var.beta)){  
        prior.var.beta <- rep(10000 , p)}
    ##Save as vector of standard deviations since this is what is taken as an argument by dnorm
    sd.beta<-sqrt(prior.var.beta)
    
    ##Prior for sigma
    prior.sigma.shape<-d
    prior.sigma.scale<-diag(x = 1, nrow=d)
    
    
    ##Prior for sigma2
    prior.sigma2<-c(0.001, 0.001)
    
    
    ###################################################
    ## Store the results
    ###################################################
    
    ##Only create a matrix of size needed minus burnin and thinning
    saved.samples<-(samples - burnin)/thin
    beta.draws <- matrix(0 , nrow=saved.samples , ncol=p*d)
    HB.draws <- matrix(0 , nrow=saved.samples , ncol=HB.p*d)
    phi.draws <- matrix(0, nrow=saved.samples, ncol=n.phi*d)
    sigma.draws <- matrix(0, nrow=saved.samples, ncol=d*d)
    rho.draws <- matrix(0, nrow=saved.samples, ncol=1)
    alpha.draws <- matrix(0, nrow=saved.samples, ncol=d)
    sigma2.draws <- matrix(0, nrow=saved.samples, ncol=d)
    deviance.draws <- matrix(0, nrow=saved.samples, ncol=1)
    like.draws <- matrix(0, nrow=saved.samples, ncol=n.all)
    ###################################################
    ## Acceptance Rates
    ###################################################
    
    ##Beta##
    ## Set up vectors for calculating acceptance rates for beta
    ##Number of blocks
    n.blocks<-ceiling(p/block)
    beta.accepttotal <- matrix(rep(0, n.blocks),nrow=n.blocks,ncol=d)
    beta.acc <- matrix(rep(0, n.blocks),nrow=n.blocks,ncol=d)
    
    ##HB##
    HB.accepttotal <- matrix(rep(0,d*2), ncol=d)
    HB.acc <- matrix(rep(0,d*2), ncol=d)
    
    ##Phi##
    ## Set up vectors for acceptance rates for phi
    phi.acc<-rep(0,4)
    
    ##Rho##
    ## Set up vectors for acceptance rates for rho
    rho.acc<-rep(0,4)
    
    #####################################################
    ##Posterior quantities
    #####################################################
    ##MIGHT NEED TO USE WHOLE XtX INVERSE THEN PICK OUT THE RIGHT VARIANCE VALUES FOR BETA AND HB DRAWS
    
    
    ##Beta##
    ##calculations needed for variance.
    proposal.var.beta<-0.0001
    proposal.corr.beta<-solve(t(X.beta) %*% X.beta ) 
    ##Set up vector to update sigma2 for each block. This will be updated throughout.
    sigma2.update<-matrix(rep(proposal.var.beta,n.blocks),nrow=n.blocks,ncol=d)
    
    ##HB##
    ##Set up vector to update sigma2 for each block. This will be updated throughout.
    sigma2.update.HB<-matrix(rep(0.0001,d), ncol=d)
    
    ##Phi##
    ##Proposal distribution sd
    proposal.sd.phi<-0.1
    
    ##sigma##
    ##scale
    Q.W<-diag(n.neighbours)-W
    
    ##shape
    shape.sigma<-prior.sigma.shape+n.phi
    
    ##Rho##
    ##Proposal distribution sd
    proposal.sd.rho<-0.01
    
    ##Eigenvalues of Q.W
    Q.eigen<-eigen(Q.W)$values
    
    ##Sigma2##
    ##shape
    shape.sigma2<-prior.sigma2[1]+S*(tp)/2
    
    R.1<-matrix(0,nrow=tp, ncol=tp)
    for(k in 1:tp){
        if(k==1){
            R.1[k,k]<-1
            R.1[k,k+1]<-1
        }
        else if(k==tp){
            R.1[k,k]<-1
            R.1[k,k-1]<-1
        }
        else{
            R.1[k,k]<-1
            R.1[k,k+1]<-1
            R.1[k,k-1]<-1
        }
    }
    #################################################################################################
    ## Set up blocking structure for beta
    #################################################################################################
    
    ##Set up start vector
    start<-rep(0, n.blocks)
    
    ##Set up end vector
    end<-rep(0, n.blocks)
    
    ##Now fill these vectors
    ##First for blocks of 5 parameters
    s<-1
    e<-block
    for(i in 1:floor(p/block)){
        start[i]<-s
        s<-s+block
        
        end[i]<-e
        e<-e+block
    }
    
    ##Now for last block when there are less parameters than blocking size
    if(p%%block !=0){
        if(p%%block == 1){
            start[n.blocks]<-p
            end[n.blocks]<-p} else{
                start[n.blocks]<-(p+1) - (p%%block)
                end[n.blocks]<- p
            }
    }
    
    
    #################################################################################################
    ## Select starting values for all quantities
    #################################################################################################
    ## Select starting values for beta and HB to be estimates from GLM
    start.values <- glm(y~X.big[,-1]+offset(offset), family="poisson")
    all.current <- as.numeric(summary(start.values)$coefficients[, 1])
    
    ##Beta##
    beta.current <- matrix(rep(all.current[1:p], d), nrow=p)
    beta.proposal <- beta.current
    
    
    ##Calculate XtransposeBeta for calcuating densities.
    XtB.curr<-X.beta%*%beta.current
    
    ##HB##
    HB.current <-matrix(rep(c(all.current[(p+1)],all.current[(p+1):(p+HB.p-1)]),d),ncol=d)
    HB.proposal <- HB.current
    HB.current.arr <- array(HB.current, c(tp,S,d))
    HB.proposal.arr <- HB.current.arr
    
    
    ## Calculate XtHB for calculating denisty in metropolis draw
    
    XtHB.curr<-X.HB%*%HB.current
    XtHB.prop<-XtHB.curr
    XtHB.arr<-array(XtHB.curr, c(n.phi,tp,  d))
    ##Phi##
    ##Starting values for Phi
    phi.current<-rnorm(n=n.phi, mean=rep(0,n.phi),sd=rep(0.1,n.phi))
    phi.current<-matrix(rep(phi.current,d), nrow=n.phi, ncol=d, byrow=T)
    
    ##Starting value for sigma
    sigma.curr<-diag(runif(1),d, nrow=d)
    sigma.inv<-solve(sigma.curr)
    
    ##Starting value for rho
    rho.current<-runif(1)
    rho.proposal<-rho.current
    
    ##Starting values for alpha
    alpha.current<-rep(runif(1),d)
    
    
    ##Starting value for sigma2
    sigma2.current<-rep(runif(1),d)
    #################################################################################################
    ##Metropolis Algorithm
    #################################################################################################
    for(i in 1:samples){
        ##run progress bar
        if(progress.bar){
            setTxtProgressBar(pb,i)}
        

            
        ###########################################################
        ##Sample Beta
        ###########################################################
        ##Loop over disease
        for(j in 1:d){
            ##Metropolis Algorithm for each block
            for(k in 1: n.blocks){
                
                ##Update variance for each block
                var.prop <-  sigma2.update[k,j]*proposal.corr.beta
                #var.prop <- sigma2.update[k,j]*diag(rep(1,p))
                if(block == 1){
                    beta.proposal[k,j]<-rnorm(1,beta.current[k,j], var.prop[k,k]^(1/2))
                } else{
                    #draw proposals for each of the beta parameters in block k
                    ##If there is only one parameter to be drawn in the last block then draw from rnorm rather than rmvnorm.
                    
                    if((k == n.blocks & p%%block == 1)){
                        beta.proposal[p,j]<-rnorm(1,beta.current[p,j], var.prop[p,p]^(1/2))} else{
                            
                            beta.proposal[start[k]:end[k],j]<-as.numeric(rmvnorm(1,beta.current[start[k]:end[k],j], var.prop[start[k]:end[k],start[k]:end[k]]))
                        }
                }
                
                
                
                #calculate acceptance ratio, r
                top<-sum(dpois(y.mat[,j], exp(X.beta%*%beta.proposal[,j]  + XtHB.curr[,j] + rep(phi.current[ ,j], tp) + offset.mat[,j]) , log=T))+sum(dnorm(beta.proposal[,j], prior.mean.beta, sd.beta, log=T))
                bottom<-sum(dpois(y.mat[,j], exp(X.beta%*%beta.current[,j]+ XtHB.curr[,j] + rep(phi.current[ ,j], tp)   + offset.mat[,j]) , log=T))+sum(dnorm(beta.current[,j], prior.mean.beta, sd.beta, log=T))
                r<-top-bottom
                
                
                #set beta.current = beta.proposal with probability r
                if(log(runif(1))<r){
                    beta.current[,j]<-beta.proposal[,j]
                    #print ("worked")
                    #count the number of accepted draws for this block
                    beta.acc[k,j]<-beta.acc[k,j]+1}else{
                        #print("didn't work")
                        beta.proposal[,j]<-beta.current[,j]}
                
            }}
        #***#
        #beta.current<-beta.sim
        #***#
        ##Calculate XtransposeBeta with current values for Beta
        #for(j in 1:d){
        #  XtB.curr[,j]<-betaX%*%beta.current[,j]
        #}
        
        XtB.curr<-X.beta%*%beta.current
        XtB.arr<-array(XtB.curr, c(n.phi,tp,  d))
        
        

        ###########################################################
        ##Sample HB
        ###########################################################
        ##Metropolis Algorithm for health board parameter
        ########################################################

        O.arr<-XtB.arr+offset.arr+array(rep(phi.current, d*tp),c(n.phi,tp,d))

        HBacc <- rep(0, d)
        HB<-HBUpdate(HB_current_arr=HB.current.arr, alpha_current=alpha.current, sigma2_current=sigma2.current,
                      sigma2_update_HB=sigma2.update.HB, X_phi=X.phi, y_arr=y.arr, O_arr=O.arr,
                      XtHB_arr=XtHB.arr,  d=d, S=S, tp=tp, HBacc=HBacc)

        HB.current.arr <- HB[[1]]
        HB.acc[1, ] <- HB.acc[1, ] + HB[[2]]
        HB.acc[2, ] <- HB.acc[2, ] + rep(S*tp,d)

        HB.current.orig<-matrix(HB.current.arr, ncol=d)
        HB.current<-sweep(HB.current.orig,2,apply(HB.current.orig,2,mean),"-")
        HB.proposal<-HB.current
        
        XtHB.curr<-X.HB%*%HB.current
        
        XtHB.arr<-array(XtHB.curr, c(n.phi,tp,  d))
        HB.current.arr<-array(HB.current, c(tp,S,d))
        HB.proposal.arr<-HB.current.arr
        
        ##############################
        # HB.current<-HB.mat
        # HB.current.arr<-HB.replace
        # HB.proposal.arr<-HB.replace
        #XtHB.curr<-XtHB.mat
        # XtHB.arr<-XtHB.replace
        ##############################
        ###########################################################
        ##Sample Phi
        ###########################################################
        #XtB<-matrix(XtB.curr, nrow=n.phi, ncol=d)
        phiacc<-c(0,0)
        #sigma.inv<-solve(sigma.curr)
        O<-XtB.arr+offset.arr+XtHB.arr
        phi <-PhiUpdate(phi_current=phi.current,rho_current=rho.current,
                        proposal_sd_phi=proposal.sd.phi, n_phi=n.phi, 
                        n_neighbours=n.neighbours,W_triplet=W.triplet, W_start=W.start,
                        y=y.arr, O=O, sigma_inv=sigma.inv, d=d, t=tp, phi_acc=phiacc)
        
        
        
        phi.current <- phi[[1]]
        phi.acc[1:2] <- phi.acc[1:2] + phi[[2]]
        
        #for(j in 1:d){
        #  phi.current[,j] <- phi.current[ ,j] - mean(phi.current[ ,j])
        #}
        
        ##This doesn't seem to be faster
        #phi.current<- sweep(phi.current,1, apply(phi.current,1,mean),"-")
        
        for(j in 1:d){
            phi.HB <- matrix(0, nrow=n.phi, ncol=S)
            ##Create an empty matrix to store the mean centred values of each of the phi's
            phi.HB.centered <- matrix(0, nrow=n.phi, ncol=S)
            for(s in 1:S){
                phi.HB[,s]<-X.phi[,s]*phi.current[,j]
            }
            ##Find the mean phi value seperately for each phi
            mean.phi.HB <- apply(phi.HB, 2, sum) / apply(X.phi,2,sum)
            
            ##Replace 0's with NA's so when the mean is subtracted 
            ##only for elements belonging to the correct HB
            phi.HB[phi.HB == 0] <- NA
            
            ##Mean centre phi for each HB
            ##WHICH IS FASTER LOOP OR SWEEP
            #phi.HB.centered<-sweep(phi.HB,2,mean.phi.HB,"-")
            
            for(s in 1:S){      
                phi.HB.centered[,s]<-phi.HB[,s]-mean.phi.HB[s]
            }
            
            
            ##Change missing values back to 0's as the sum of NA's is NA
            phi.HB.centered[is.na(phi.HB.centered)] <- 0
            
            
            ##We want a vector of phis, not a matrix.
            ##For each row of this matrix, only one element
            ##should be non-zero, so summing across the
            ##rows will create a vector with the mean centred 
            ##values of phi.
            phi.current[,j]<-apply(phi.HB.centered, 1, sum)
        }
        
        
        #phi.current<-phi.sim.mat
        ###########################################################
        ##Sample Sigma
        ###########################################################
        Q.R.W<-rho.current*Q.W+(1-rho.current) * diag(rep(1,n.phi))
        scale.sigma<-prior.sigma.scale+t(phi.current)%*%(Q.R.W)%*%phi.current
        sigma.curr<-riwish(shape.sigma,scale.sigma)
        
        ###########################################################
        ##Sample Rho
        ###########################################################
        ##Draw proposal
        rho.proposal<-rtrunc(1,spec="norm",a=0,b=1,mean=rho.current,sd=proposal.sd.rho)
        
        ##Caluculate proposal quantities
        sigma.inv<-solve(sigma.curr)
        
        
        kr<-KroneckerUpdate(phi_current=phi.current, rho_current=rho.current, rho_proposal=rho.proposal,
                            n_phi=n.phi, triplet_length=triplet.length, n_neighbours=n.neighbours, W_triplet=W.triplet, 
                            sigma_inv=sigma.inv)
        
        current.scale.tau.part<-kr[[1]]
        proposal.scale.tau.part<-kr[[2]]
        
        top<-(d/2)*sum(log(rho.proposal*Q.eigen + (1-rho.proposal))) - proposal.scale.tau.part
        bottom<-(d/2)*sum(log(rho.current*Q.eigen + (1-rho.current))) - current.scale.tau.part
        
        r=top-bottom
        
        
        #set rho.current = rho.proposal with probability r
        if(log(runif(1))<r){
            rho.current<-rho.proposal
            rho.acc[1]<-rho.acc[1]+1}else{
                rho.proposal<-rho.current}
        rho.acc[2]<-rho.acc[2]+1
        
        ###########################################################
        ##Sample sigma2
        ###########################################################
        for(j in 1:d){
            R<-R.1*(-alpha.current[j])
            diag(R)<-c((alpha.current[j]^2),rep((1+alpha.current[j]^2),tp-2),1)
            scale.sigma2<-0
            for(s in 1:S){
                scale.sigma2<-scale.sigma2+t(HB.current.arr[,s,j])%*%R%*%(HB.current.arr[,s,j])
            }
            scale.sigma2<-prior.sigma2[2]+0.5*scale.sigma2
            
            
            
            sigma2.current[j]<-rinvgamma(1,shape.sigma2,scale.sigma2)
            
            ###########################################################
            ##Sample alpha
            ###########################################################
            if(tp==2){
              top.alpha<-sum(HB.current.arr[(2:tp),,j]*HB.current.arr[(1:(tp-1)),,j])
              bottom.alpha <-sum(HB.current.arr[(1:(tp-1)),,j]^2)
            }
            else{
              top.alpha <- sum(HB.current.arr[(2:tp),,j]*HB.current.arr[(1:(tp-1)),,j])
              bottom.alpha <- sum(HB.current.arr[(1:(tp-1)),,j]^2)
            }
            
            mean.alpha<-top.alpha/bottom.alpha
            sd.alpha<-sqrt(sigma2.current[j]/bottom.alpha)
            
            alpha.current[j]<-rtrunc(1, spec="norm", a=0, b=1, mean=mean.alpha, sd=sd.alpha)
        }
        
        #sigma2.current<-sigma2.sim
        #alpha.current<-alpha.sim

        
       
        ###########################################################
        ##Calculate deviance
        ###########################################################
        phi.current.extend <- array(NA, c(tp*n.phi, d))
        for(j in 1:d)
        {
            phi.current.extend[ ,j] <- rep(phi.current[ ,j], tp)    
        }
        
        fitted <- exp(as.vector(XtB.curr) + as.vector(XtHB.curr) + as.vector(phi.current.extend)+ offset)
        deviance.all <- dpois(y, lambda=fitted, log=TRUE)
        deviance <- -2 * sum(deviance.all, na.rm=TRUE)
        like <- exp(deviance.all)
        
        #################################################################################################
        ##Store Results
        #################################################################################################
        ##Store the results from each draw in main beta matrix
        ##Do not store if sample part of burn in
        ##Do not store unless multiple of thinning parameter
        if(i>burnin & i%%thin==0){
            row <- (i-burnin)/thin
            beta.draws[row, ] <- as.numeric(beta.current)
            HB.draws[row, ] <- as.numeric(HB.current)
            phi.draws[row, ] <- as.numeric(phi.current)
            sigma.draws[row ,] <- as.vector(sigma.curr)
            rho.draws[row, ] <- rho.current
            sigma2.draws[row, ] <- sigma2.current
            alpha.draws[row, ] <- alpha.current
            deviance.draws[row ,]<-deviance
            like.draws[row,]<-like
        }
        
        #################################################################################################
        ##Alter Variance
        #################################################################################################
        #If iteration number a multiple of 100, alter the variance if acceptance rates are too high or too low
        #For a multi-dimensional proposal we want the acceptance rate to be around 25%
        #for a 1-dimensional proposal we want the acceptance rate to be around 44%
        t <- i/100
        if(ceiling(t)==floor(t))
        {
            
            ##BETA
            for(j in 1:d){
                beta.accepttotal[1:n.blocks,j] <- beta.accepttotal[1:n.blocks,j] + beta.acc[,j]
                
                for (k in 1:n.blocks){
                    #If acceptance rates far too high (>0.5), double variance
                    if(beta.acc[k,j]/100 >= 0.4){
                        sigma2.update[k,j]<-2*sigma2.update[k,j]}
                    
                    #If acceptance rates between 0.4 and 0.5, add 10% of current value on to variance
                    if(beta.acc[k,j]/100 >= 0.3 & beta.acc[k,j]/100 < 0.4){
                        sigma2.update[k,j]<-sigma2.update[k,j] + 0.1*sigma2.update[k,j]
                    }
                    
                    #If acceptance rates far too low (<0.1), half variance
                    if(beta.acc[k,j]/100 <= 0.1){
                        sigma2.update[k,j]<-sigma2.update[k,j]/2
                    }
                    
                    #If acceptance rates between 0.1 and 0.2, take 10% of current value from variance
                    if(beta.acc[k,j]/100 > 0.1 & beta.acc[k,j]/100 <= 0.2){
                        sigma2.update[k,j]<-sigma2.update[k,j] - 0.1*sigma2.update[k,j]
                    }
                }}
            ## reset count of accepted draws to be zero
            beta.acc <- matrix(rep(0, n.blocks),nrow=n.blocks,ncol=d)
            
            
            ##HB
            for(k in 1:d){
                HB.accepttotal[1,k] <- HB.accepttotal[1,k] + HB.acc[1,k]
                HB.accepttotal[2,k] <- HB.accepttotal[2,k] + HB.acc[2,k]
                
                #If acceptance rates far too high (>0.6), double variance
                if(HB.acc[1,k]/HB.acc[2,k] >= 0.6){
                    sigma2.update.HB[k]<-2*sigma2.update.HB[k]}
                
                #If acceptance rates between 0.5 and 0.6, add 10% of current value on to variance
                if(HB.acc[1,k]/HB.acc[2,k] >= 0.5 & HB.acc[1,k]/HB.acc[2,k] < 0.6){
                    sigma2.update.HB[k]<-sigma2.update.HB[k] + 0.1*sigma2.update.HB[k]
                }
                
                #If acceptance rates far too low (<0.3), half variance
                if(HB.acc[1,k]/HB.acc[2,k] <= 0.3){
                    sigma2.update.HB[k]<-sigma2.update.HB[k]/2
                }
                
                #If acceptance rates between 0.3 and 0.4, take 10% of current value from variance
                if(HB.acc[1,k]/HB.acc[2,k] > 0.3 & HB.acc[1,k]/HB.acc[2,k] <= 0.4){
                    sigma2.update.HB[k]<-sigma2.update.HB[k] - 0.1*sigma2.update.HB[k]
                }
                
                
                
            }
            HB.acc <- matrix(rep(0,d*2), ncol=d)
            #print(sigma2.update.HB)
            ##PHI
            #If acceptance rates far too high (>0.6), double variance
            if(phi.acc[1]/phi.acc[2]>= 0.6){
                proposal.sd.phi<-proposal.sd.phi*2
            }
            #If acceptance rates between 0.5 and 0.6, add 10% of current value on to variance
            
            if(phi.acc[1]/phi.acc[2] >= 0.5 & phi.acc[1]/phi.acc[2] < 0.6){
                proposal.sd.phi<-proposal.sd.phi + 0.1*proposal.sd.phi
            }
            #If acceptance rates far too low (<0.1), half variance
            if(phi.acc[1]/phi.acc[2] <= 0.3){
                proposal.sd.phi<-proposal.sd.phi/2
            }
            #If acceptance rates between 0.1 and 0.2, take 10% of current value from variance
            if(phi.acc[1]/phi.acc[2] > 0.3 & phi.acc[1]/phi.acc[2] <= 0.4){
                proposal.sd.phi<-proposal.sd.phi - 0.1*proposal.sd.phi
            }
            phi.acc[3]<-phi.acc[1]+phi.acc[3]
            phi.acc[4]<-phi.acc[2]+phi.acc[4]
            phi.acc[1:2]<-rep(0,2)
            
            ##RHO
            #If acceptance rates far too high (>0.5), double variance
            if(rho.acc[1]/rho.acc[2]>= 0.5){
                proposal.sd.rho<-proposal.sd.rho*2
            }
            #If acceptance rates between 0.4 and 0.5, add 10% of current value on to variance
            if(rho.acc[1]/rho.acc[2] >= 0.4 & rho.acc[1]/rho.acc[2] < 0.5){
                proposal.sd.rho<-proposal.sd.rho + 0.1*proposal.sd.rho
            }
            #If acceptance rates far too low (<0.1), half variance
            if(rho.acc[1]/rho.acc[2] <= 0.1){
                proposal.sd.rho<-proposal.sd.rho/2
            }
            #If acceptance rates between 0.1 and 0.2, take 10% of current value from variance
            if(rho.acc[1]/rho.acc[2] > 0.1 & rho.acc[1]/rho.acc[2] <= 0.2){
                proposal.sd.rho<-proposal.sd.rho - 0.1*proposal.sd.rho
            }
            rho.acc[3]<-rho.acc[1]+rho.acc[3]
            rho.acc[4]<-rho.acc[2]+rho.acc[4]
            rho.acc[1:2]<-rep(0,2)
            
        }
        
        ##Save results as function is running to avoid losing everything if not enough space R
        # results<-list(list(beta=beta.draws, phi=phi.draws, sigma=sigma.draws, accept.beta=beta.accepttotal/samples, 
        #                     accept.phi=phi.acc[3]/phi.acc[4]))
        
        #save(results, file="E:/results.RData")
    }
    
    #################################################################################################
    ## Deviance information criterion (DIC)
    #################################################################################################
    median.beta <- matrix(apply(beta.draws, 2, median), nrow=p, ncol=d)
    beta.fit<-X.beta%*%median.beta
    
    median.HB <- matrix(apply(HB.draws, 2, median), nrow=HB.p, ncol=d)
    HB.fit <- X.HB%*%median.HB
    
    median.phi <- matrix(apply(phi.draws, 2, median), nrow=n.phi)
    median.phi.extend <- array(NA, c(tp*n.phi, d))
    for(j in 1:d)
    {
        median.phi.extend[ ,j] <- rep(median.phi[ ,j], tp)    
    }
    
    fitted.median <- exp(as.vector(beta.fit) + as.vector(HB.fit) + as.vector(median.phi.extend) + offset)
    deviance.all.fitted<-dpois(y, lambda=fitted.median, log=TRUE)
    deviance.fitted <- -2 * sum(deviance.all.fitted, na.rm=TRUE)
    p.d <- median(deviance.draws) - deviance.fitted
    DIC <- 2 * median(deviance.draws) - deviance.fitted 
    
    
    #################################################################################################
    ## Watanabe-Akaike Information Criterion (WAIC)
    #################################################################################################
    LPPD <- sum(log(apply(like.draws,2,mean)), na.rm=TRUE)
    p.w <- sum(apply(log(like.draws),2,var), na.rm=TRUE)
    WAIC <- -2 * (LPPD - p.w)
    
    #################################################################################################
    ##Back Transform intercept and continuous covariates
    #################################################################################################
    ##Intercept
    ##This is assuming that the same covariates are used for each disease.
    col<-1
    intercept<-1
    for(j in 1:d){
        intercept.terms<-matrix(0 , nrow=(samples-burnin)/thin , ncol=(p))
        for(i in 1:(p)){
            if(data.type[i]==1){
                intercept.terms[,i]<-beta.draws[,col]*data.means[i]/data.sd[i]
                
            }
            col<-col+1
        }
        sum.intercept<-apply(intercept.terms,1,sum)
        beta.draws[,intercept]<-beta.draws[,intercept]-sum.intercept
        intercept<-intercept+p
    }
    
    
    
    
    
    ##Continuous covariates
    for(i in 1:(p*d)){
        if(data.type[i]==1){
            beta.draws[,i]<-beta.draws[,i]/data.sd[i]
        } 
    }
    
    ##calulate time
    t<-c(t,difftime(Sys.time(),t0,u="mins"))
    cat("\n")
    cat("Time to run:",t[2], "minutes")
    #################################################################################################
    ##Return results
    #################################################################################################
    
    ##Return beta vector and acceptance rates
    return(list(beta=beta.draws, HB=HB.draws, phi=phi.draws, sigma=sigma.draws, 
                rho=rho.draws, sigma2=sigma2.draws, alpha=alpha.draws,
                accept.beta=beta.accepttotal/samples, accept.HB=HB.accepttotal[1,]/HB.accepttotal[2,],
                accept.phi=phi.acc[3]/phi.acc[4],DIC=DIC, WAIC=WAIC, 
                p.d=p.d, p.w=p.w, t))
    
    ##Close Progress Bar
    if(progress.bar){
        close(pb)}
    
}
