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
        
        p_est =special.gammainc(x_est,gam)

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