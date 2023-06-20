import numpy as np
import math
import warnings
warnings.filterwarnings("ignore")
from numpy import polyfit

# y = kx + b
# RETURN: [k,b],[k_err,b_err],r^2

def polyfit_func(X, Y, Y_ERR=None, degree=1, lims=None):
    if lims==None:
        if Y_ERR==None:
            x = []; y = []; y_err = None
            for k in range(len(X)):
                if not np.isnan(X[k]) and not np.isnan(Y[k]):
                    x.append(X[k])
                    y.append(Y[k])
        else:
            x = []; y = []; y_err = []
            for k in range(len(X)):
                if not np.isnan(X[k]) and not np.isnan(Y[k]) and not np.isnan(Y_ERR[k]):
                    x.append(X[k])
                    y.append(Y[k])
                    y_err.append(Y_ERR[k])
    else:
        if Y_ERR==None:
            x = []; y = []; y_err = None
            for k in range(len(X)):
                if X[k]>lims[0] and X[k]<lims[1] and  Y[k]>lims[2] and Y[k]<lims[3]:
                    x.append(X[k])
                    y.append(Y[k])
        else:
            x = []; y = []; y_err = []
            for k in range(len(X)):
                if X[k]>lims[0] and X[k]<lims[1] and  Y[k]>lims[2] and Y[k]<lims[3] and not np.isnan(Y_ERR[k]):
                    x.append(X[k])
                    y.append(Y[k])
                    y_err.append(Y_ERR[k])    


    results = {}

    if y_err==None:
        p, V = np.polyfit(x, y, degree, cov=True)
    else:
        p, V = np.polyfit(x, y, degree, w = [1.0 / ty for ty in y_err], cov=True)
    
    
    #print 'ord',p
    # Polynomial Coefficients
    results['pol_coeffs'] = p.tolist()

    # r-squared
    P = np.poly1d(p)
    # fit values, and mean
    yhat = P(x)                         # or [P(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    results['R^2'] = ssreg / sstot
    
    errs = []
    print 
    for k in range(len(p)):
       errs.append(np.sqrt(V[k][k])) 
    
    results['pol_coeffs_std_err'] = errs

    return results['pol_coeffs'],results['pol_coeffs_std_err'],results['R^2']


def sm_liniar(x, y, y_err=None):
    import pandas as pd
    import statsmodels.formula.api as sm

    ws = pd.DataFrame({
        'x': x,
        'y': y
    })
    weights = pd.Series(y_err)
    if y_err!=None:
        wls_fit = sm.wls('y ~ x', data=ws, weights=1 / weights**2).fit()
        return [wls_fit.params[1],wls_fit.params[0]],[wls_fit.bse[1],wls_fit.bse[0]],wls_fit.rsquared

    else:
        ols_fit = sm.ols('y ~ x', data=ws).fit()
        return [ols_fit.params[1],ols_fit.params[0]],[ols_fit.bse[1],ols_fit.bse[0]],ols_fit.rsquared
    #print wls_fit.summary()
    #print ols_fit.summary()


'''
# EXAMPLE:
x = [0.3333333333333333, 0.2886751345948129, 0.25, 0.23570226039551587, 0.22360679774997896, 0.20412414523193154, 0.2, 0.16666666666666666]
y = [0.13250359351851854, 0.12098339583333334, 0.12398501145833334, 0.09152715, 0.11167239583333334, 0.10876248333333333, 0.09814170444444444, 0.08560799305555555]
y_err = [0.003306749165349316, 0.003818446389148108, 0.0056036878203831785, 0.0036635292592592595, 0.0037034897788415424, 0.007576672222222223, 0.002981084130692832, 0.0034913019065973983]

res = sm_liniar(x, y, y_err)
print res
res = polyfit_func(x, y, y_err=y_err, degree=1)
print res

res = sm_liniar(x, y)
print res
res = polyfit_func(x, y, y_err=None, degree=1)
print res
'''
