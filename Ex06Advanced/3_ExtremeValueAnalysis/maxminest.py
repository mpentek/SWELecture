def maxminest (record, dur_ratio = 1):

    import numpy as np
    import scipy.interpolate as interpolate
    import math
    #    if not record:


    n_cdf_pk =1000
    cdf_pk_min = 0.0005
    cdf_pk_max = 0.9995

    cdf_pk = np.linspace(cdf_pk_min,cdf_pk_max,n_cdf_pk)

    rsize = np.array(record).shape
    
    if len(rsize) == 1:
        rec_size = 1
    else:
        rec_size = rsize[0]
    

    max_est = np.zeros((rec_size,1))
    min_est = np.zeros((rec_size,1))
    max_std = np.zeros((rec_size,1))
    min_std = np.zeros((rec_size,1))

    for i in np.arange(rec_size):
        if rec_size == 1:
            x = record
        else:
            x = record[:,i]
        
        n = x.size
        
        mean_x = np.mean(x)

        # std_x = np.std(x,ddof = 1) # Original
        std_x = np.std(x)

        skew_x = np.sum(np.power(x-mean_x,3))/(n*std_x**3)

        X = x*np.sign(skew_x)
        
        sort_X = np.sort(X)

        mean_X = mean_x*np.sign(skew_x)
        std_X = std_x
        CDF_X = np.divide(np.arange(1,n+1),n+1)

        n_coarse = min([n,1000])

        CDF_coarse = np.linspace(1/(n_coarse+1),n_coarse/(n_coarse+1),n_coarse)
        
        f = interpolate.interp1d(CDF_X,sort_X)
        X_coarse = f(CDF_coarse)
        
        
        
        
        mean_X_coarse = np.mean(X_coarse)
        # std_X_coarse = np.std(X_coarse) # Original
        std_X_coarse = np.std(X_coarse, ddof=1)

        gamma_min = 1
        gamma_max = 125
        n_gamma = 19
        n_start = 7

        gamma_list = np.logspace(math.log10(gamma_min),math.log10(gamma_max),n_gamma)
        

        gam_PPCC_list = np.zeros(gamma_list.shape)
        count = 0
        beta_coarse_list = np.zeros((125,1))
        mu_coarse_list =np.zeros((125,1))
        
        for j in np.arange(n_start,-1,-1):
          
            count+=1
            
            
            s_gam_j = stdgaminv(CDF_coarse,gamma_list[j])
            
            mean_s_gam_j = np.mean(s_gam_j)
            print(mean_s_gam_j)
            # linear regression:
            
            beta_coarse_list[j] = (np.sum(np.multiply(s_gam_j,X_coarse))-(n_coarse*mean_s_gam_j*mean_X_coarse))/(np.sum(np.power(s_gam_j,2))-(n_coarse*mean_s_gam_j**2))
            
            mu_coarse_list[j]=(mean_X_coarse - beta_coarse_list[j]*mean_s_gam_j)
            #print('mean = {}'.format(mean_s_gam_j))
            #print('beta_coarse_list[{}] = {}'.format(j,beta_coarse_list[j]))
            #print('mu_coarse_list[{}] = {}'.format(j,mu_coarse_list[j]))
            #Probability Plot Correlation Coefficient:
           
            # gam_PPCC_list[j] = (beta_coarse_list[j]*np.std(s_gam_j)/std_X_coarse) # Original
            gam_PPCC_list[j] = (beta_coarse_list[j]*np.std(s_gam_j, ddof=1)/std_X_coarse)
            #print(gam_PPCC_list[j])
            X_coarse_fit_j = mu_coarse_list[j] + beta_coarse_list[j]*s_gam_j

            if gam_PPCC_list[j] == max(gam_PPCC_list):
                gam = gamma_list[j]
                gam_PPCC_max = gam_PPCC_list[j]
            else:
                break
        
        if gam_PPCC_list[n_start-1] < gam_PPCC_list[n_start]:
            # if the PPCC decreased with decreasing gamda, try increasing gamma: 
            for j in np.arange(n_start+1,n_gamma):
                count += 1
                # Obtain the Gamma Distribution Parameters for current gamma:
                
                s_gam_j = stdgaminv(CDF_coarse,gamma_list[j])   # standard variate
                mean_s_gam_j = np.mean(s_gam_j)
                # linear regression:
                beta_coarse_list[j] = (np.sum(np.multiply(s_gam_j,X_coarse))-(n_coarse*mean_s_gam_j*mean_X_coarse))/(np.sum(np.power(s_gam_j,2))-(n_coarse*mean_s_gam_j**2))
                
                mu_coarse_list[j] = mean_X_coarse - beta_coarse_list[j]*mean_s_gam_j
                #Probability Plot Correlation Coefficient:
                # gam_PPCC_list[j] = beta_coarse_list[j]* np.std(s_gam_j)/std_X_coarse # Original
                gam_PPCC_list[j] = beta_coarse_list[j]* np.std(s_gam_j, ddof=1)/std_X_coarse
                X_coarse_fit_j = mu_coarse_list[j] + beta_coarse_list[j]*s_gam_j

                ##
                # BLOCK needs extra indent not to break out prematurely
                ##
                if gam_PPCC_list[j] == max(gam_PPCC_list):
                    gam = gamma_list[j]
                    gam_PPCC_max = gam_PPCC_list[j]
                else:
                    break
                ##
                # BLOCK needs extra indent
                ##

            ##
            # ORIGINAL
            ##
            # if gam_PPCC_list[j] == max(gam_PPCC_list):
            #     gam = gamma_list[j]
            #     gam_PPCC_max = gam_PPCC_list[j]
            # else:
            #     break
            ##
            # ORIGINAL
            ##

        s_gam = stdgaminv(CDF_X,gam)
        mean_s_gam = np.mean(s_gam)
        print(mean_s_gam)
        beta = (np.sum(np.multiply(s_gam,sort_X))-n*mean_s_gam*mean_X)/(np.sum(np.power(s_gam,2))-n*mean_s_gam**2)
        mu = mean_X - beta*mean_s_gam
        # gam_PPCC = beta*np.std(s_gam)/std_X # Original
        gam_PPCC = beta*np.std(s_gam, ddof=1)/std_X
        
        x_fit = mu +beta*s_gam

        # Obtain the Normal Distribution Parameters for lower portion of CDF

        CDF_split = 0.25
        f = interpolate.interp1d(CDF_X,sort_X)
        X_split = f(CDF_split)

        ind_low = np.where(sort_X<X_split)
        
        X_low = sort_X[ind_low]
        n_low = len(X_low)
        CDF_low = CDF_X[ind_low]

        s_norm_low = stdnorminv(CDF_low)
        mean_s_norm_low = np.mean(s_norm_low)
        mean_X_low = np.mean(X_low)

        # linear regression:
        

                
        sigma_low = (np.sum(np.multiply(s_norm_low,X_low))-n_low*mean_s_norm_low*mean_X_low)/(np.sum(np.power(s_norm_low,2))-n_low*mean_s_norm_low**2)
        
        mu_low=mean_X_low - sigma_low*mean_s_norm_low
        X_low_fit = mu_low +sigma_low*s_norm_low

        # Probability Plot Correlation Coefficient:

        # norm_PPCC = sigma_low*np.std(s_norm_low)/np.std(X_low) # Original
        norm_PPCC = sigma_low*np.std(s_norm_low, ddof=1)/np.std(X_low, ddof=1)

        X_u=np.mean(sort_X[np.where(abs(CDF_X-0.5) == min(abs(CDF_X-0.5)))])

        front = np.where(X[1:]>=X_u)
        back = np.where(X[0:-1]<X_u)
        
        Nupcross = len(set(front[0]) & set(back[0]))
        

        if Nupcross<100:
            print('The number of median upcrossings is low {}'.format(Nupcross))
            print('The record may be too short for accurate peak estimation.')
        
        y_pk = np.sqrt(2.0*np.log(np.divide(-dur_ratio*Nupcross,np.log(cdf_pk))))
        
        CDF_y = stdnormcdf(y_pk)
        
        
        # Perform the mapping procedure to compute the CDF of largest peak for X(t) from y(t)

        X_max = stdgaminv(CDF_y,gam) * beta 
        X_max+= + mu
        
        print(np.mean(beta))
        
        X_min = np.multiply(stdnorminv(1-CDF_y),sigma_low)
        
        X_min+=mu_low
        pdf_pk = np.multiply(np.multiply(-y_pk,cdf_pk),np.log(cdf_pk))
        
        
        
        # Compute the Mean of the Peaks for process X(t)

        if np.sign(skew_x)>0:
            max_est[i] = np.trapz((np.multiply(pdf_pk,X_max)),y_pk)
            min_est[i] = np.trapz((np.multiply(pdf_pk,X_min)),y_pk)
            max_std[i] = np.trapz((np.multiply(np.power((X_max-max_est[i]),2),pdf_pk)),y_pk)
            min_std[i] = np.trapz((np.multiply(np.power((X_min-min_est[i]),2),pdf_pk)),y_pk)
        else:
            ##
            # ORIGINAL
            ##
            # max_est[i] = np.trapz((np.multiply(pdf_pk,X_max)),y_pk)
            # min_est[i] = np.trapz((np.multiply(pdf_pk,X_min)),y_pk)
            ##
            # ORIGINAL
            ##

            ##
            # UPDATE according to initial MATLAB -> seems to be able to robustly handle
            # normal random as well
            ##
            max_est[i] = -np.trapz((np.multiply(pdf_pk,X_min)),y_pk)
            min_est[i] = -np.trapz((np.multiply(pdf_pk,X_max)),y_pk)
            ##
            # UPDATE according to initial MATLAB -> seems to be able to robustly handle
            # normal random as well
            ##
            max_std[i] = np.trapz((np.multiply(np.power((-X_min-max_est[i]),2),pdf_pk)),y_pk)
            min_std[i] = np.trapz((np.multiply(np.power((-X_max-min_est[i]),2),pdf_pk)),y_pk)

    return max_est, min_est, max_std, min_std