def stdgaminv(p,gam):
    import numpy as np
    import scipy.special as special
    import scipy.interpolate as interpolate
    import math
    abs_tol = 10**-3
    rel_tol = 10**-3
    
    if gam<0.1 or gam>150:
        raise ValueError('The shape parameter gamma must be between 0.1 and 150')
    
    
    
    
    
    
    p = np.array(p)
    
    x_max = 10**np.polyval([-0.009486738 ,0.03376901, 0.1151316, 0.2358172, 1.139717],math.log10(gam))
    

    
    max_iter = 200

    current_iter = 0
    
   
    while special.gammainc(gam,x_max)<max(p):
        

        current_iter+=1
        
        if current_iter>max_iter:
            raise ValueError('Maximum specified probability is too high:{}'.format(max(p)))
        else:
            x_max *=1.5
            
    x_min = 10**np.polyval([-0.0854665, 0.866249, -3.25511, 5.14328, -0.90924, -8.09135, 12.3393, -5.89628],math.log10(gam))
    current_iter = 0
    
    while special.gammainc(gam,x_min)>min(p):
        current_iter +=1
        
        if current_iter>max_iter:
            raise ValueError('Minimum specified probability is too low:{}'.format(min(p)))
        else:
            x_min *=0.1
    
    
    n_check = 1000
    x_check = np.linspace(x_min,x_max,n_check)
    

    p_check = special.gammainc(gam,x_check)
    
    
    
    
    

    p_check, ind_u = np.unique(p_check,return_index = True)
    
    
    x_check = x_check[ind_u]

    f = interpolate.interp1d(p_check,x_check,fill_value='extrapolate')

    x_est = f(p)
    
    
    max_iter = 15
    current_iter= 0
    x_step = np.ones(x_est.shape)
    
    
    
    while any(abs(x_step)>abs_tol) and any(abs(np.divide(x_step,x_est)>rel_tol)):
        current_iter+=1
        
        if current_iter>max_iter:
            break
        
        p_est =special.gammainc(gam,x_est)
        

        p_check = np.append(p_check,p_est)
        x_check = np.append(x_check,x_est)
        p_check, ind_u = np.unique(p_check,return_index = True)

        x_check = x_check[ind_u]

        f = interpolate.interp1d(p_check,x_check,fill_value='extrapolate')
        
        x_interp = f(p)

        x_step = x_interp-x_est

        x_est = x_interp
        
    x = x_est.reshape(p.size)

    
    
    return x

def stdnorminv(p):
    import numpy as np
    import scipy.special as special
    x = -1*np.sqrt(2)*special.erfcinv(2*p)

    return x

def stdnormcdf(x):
    import scipy.special as special
    import math
    p = 0.5 *special.erfc(-x/math.sqrt(2))
    return p


def maxminest (record, dur_ratio = 1):

    import numpy as np
    import scipy.interpolate as interpolate
    import math
    #    if not record:


    n_cdf_pk =1000
    cdf_pk_min = 0.025
    cdf_pk_max = 0.975

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

            # linear regression:
            
            beta_coarse_list[j] = (np.sum(np.multiply(s_gam_j,X_coarse))-(n_coarse*mean_s_gam_j*mean_X_coarse))/(np.sum(np.power(s_gam_j,2))-(n_coarse*mean_s_gam_j**2))
            
            mu_coarse_list[j]=(mean_X_coarse - beta_coarse_list[j]*mean_s_gam_j)

            #Probability Plot Correlation Coefficient:
           
            # gam_PPCC_list[j] = (beta_coarse_list[j]*np.std(s_gam_j)/std_X_coarse) # Original
            gam_PPCC_list[j] = (beta_coarse_list[j]*np.std(s_gam_j,ddof=1)/std_X_coarse)

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

        beta = (np.sum(np.multiply(s_gam,sort_X))-n*mean_s_gam*mean_X)/(np.sum(np.power(s_gam,2))-n*mean_s_gam**2) #0.12
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

        