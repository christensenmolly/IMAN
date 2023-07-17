# This functions move instrumental system to standsrt system
from useful_functions import *
from scipy.optimize import curve_fit
from sklearn.linear_model import RANSACRegressor, HuberRegressor
import warnings
from scipy import odr
import stats

def lin(B, x):
    return x + B[0]

def fit(x, y, n=1):
    
    data = odr.Data(x, y)
    results0 = odr.ODR(data=data, model=odr.Model(lin), beta0=[33.0]).run()
    results = results0
    for _ in range(n):
        y_sigma = np.sqrt(np.abs((y - lin(results.beta, x))))
        results = odr.ODR(data=odr.RealData(x=x, y=y, sy=y_sigma),
                              model=odr.Model(lin),
                              beta0=results.beta).run()
    return results.beta, results0.sd_beta


def clean(inst, calib):
    data = calib - inst

    Q1 = np.quantile(data, q=.25)
    Q3 = np.quantile(data, q=.75)
    IQR = Q3 - Q1#data.apply(stats.iqr)

    inds = np.where(~((data < (Q1-1.5*IQR)) | (data > (Q3+1.5*IQR))))

    return inst[inds], calib[inds]

def move_phot_sys(fnameB, fnameV, fnameR, pb0, pb1, pb2, out_dir):

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')

        hdu = fits.open(fnameB)
        dataB = -2.5*np.log10(hdu[0].data)
        headerB = hdu[0].header

         
        hdu = fits.open(fnameV)
        dataV = -2.5*np.log10(hdu[0].data)
        headerV = hdu[0].header
        

        hdu = fits.open(fnameR)
        dataR = -2.5*np.log10(hdu[0].data)
        headerR = hdu[0].header

        c0  = pb0[0]
        cv  = pb0[1]
        c1  = 1/pb1[0]
        cbv = -pb1[1]/c1
        c2  = 1/pb2[0]
        cvr = -pb2[1]/c2

        
        print('c0', c0)
        print('c1', c1)
        print('c2', c2)
        print('cbv', cbv)
        print('cvr', cvr)
        print('cv', cv)

        V = dataV + c0*(dataB-dataV) + cv
        BV = c1*(dataB-dataV) + cbv
        VR = c2*(dataV-dataR) + cvr

        B = BV + V
        R = V - VR

        Iv = 10**(-0.4*V)
        Iv[np.where(np.isnan(Iv))] = 0
        Ib = 10**(-0.4*B)
        Ib[np.where(np.isnan(Ib))] = 0
        Ir = 10**(-0.4*R)
        Ir[np.where(np.isnan(Ir))] = 0
        

        nameB = fnameB.name #.split('/')[-1]
        fits.writeto(Path(out_dir, nameB), Ib, header=headerB, overwrite=True) 

        nameV = fnameV.name #.split('/')[-1]
        fits.writeto(Path(out_dir, nameV), Iv, header=headerV, overwrite=True) 

        nameR = fnameR.name #.split('/')[-1]
        fits.writeto(Path(out_dir, nameR), Ir, header=headerR, overwrite=True) 
        return

def linear(B, x):
    return B[0]*x + B[1]



def get_equals(cat_B, cat_V, cat_R, est_B, est_V, est_R, fnameB, fnameV, fnameR, filters, out_dir, inds0):
    print(filters)

    gal = fnameB.stem[:-1] #fnameB.split('/')[-1].split('-')[0]
    #plot_sky(ima, 0.1, 0.1, xy=list(zip(xB, yB)), xy1=list(zip(xV, yV)), xy2=list(zip(xR, yR)))
    #plt.savefig(gal+'.png')
    
    inds = np.where( ~(np.isnan(est_B) | np.isnan(est_V) | np.isnan(est_R) ))
    est_B = est_B[inds]
    est_V = est_V[inds]
    est_R = est_R[inds]
    
    cat_B = cat_B[inds]
    cat_V = cat_V[inds]
    cat_R = cat_R[inds]
    
    
    def lin(x, *p):
        return p[0]*x + p[1]

    #def lin0(x, *p):
    #    return x + p[0]
    
    #def lin00(x, *p):
    #    return p[0]
    
    bv = est_B - est_V
    BV = cat_B - cat_V
    

    vr = est_V - est_R
    VR = cat_V - cat_R
    
    Vv = cat_V - est_V 
    
    #inds = np.where( (-2 <VR) * (VR < 2))
    #BV = BV[inds]
    #bv = bv[inds]
    #VR = VR[inds]
    #vr = vr[inds]
    #Vv = Vv[inds]
    
    from scipy.optimize import curve_fit
    #pb10, covpb = curve_fit(lin0, xdata=BV, ydata=bv, p0=np.array([0])) #, bounds=np.array([[0.6, -np.inf], [1.5,np.inf]]))
    #sy = abs(BV + pb10[0] - bv)
   
    #for _ in range(1):
    pb1, covpb = curve_fit(lin, xdata=BV, ydata=bv, p0=np.array([1,0])) #, sigma=sy) #, bounds=np.array([[0.6, -np.inf], [1.5,np.inf]]))
    #sy = abs(BV + pb1[1] - bv)
    X = [[a] for a in BV]
    y = bv
    #reg = RANSACRegressor().fit(X, y)
    reg = HuberRegressor().fit(X, y)
    s1, i1 = reg.coef_, reg.intercept_
        
    #s = pb1[0]
    #i = pb1[1]
    s_err = np.sqrt(covpb[0,0])
    i_err = np.sqrt(covpb[0,0])
    
    a, b = np.min(BV), np.max(BV)
    
    plt.figure(dpi=300, figsize=(10,10))
    #plt.title(gal)
    plt.plot(BV, bv, 'ob')
    #for eb, cb, k in zip(BV, bv, inds0):
    #    plt.plot(eb, cb, 'ob')
    plt.plot([a,b], [a*s1+i1, b*s1+i1], '-r' , label=r'$c_1 = %5.3f \pm %5.3f$ $c_{bv} = %5.3f \pm %5.3f$' %(s1, s_err, i1, i_err))
    plt.xlabel('(' + filters[0] + '-' + filters[1] + ') catalog', fontsize=15)
    plt.ylabel('('+filters[0] +' - '+ filters[1]+') inst', fontsize=15)
    #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.legend(fontsize=15)
    plt.gca().invert_yaxis()
    #plt.gca().invert_xaxis()
    plt.tight_layout()
    plt.savefig(Path(out_dir, gal+'%s%s.png' %(filters[0], filters[1])))
    plt.clf()
    plt.close()
    #plt.show()
    
    #pb20, covpb = curve_fit(lin0, xdata=VR, ydata=vr, p0=np.array([0])) #, bounds=np.array([[0.6, -np.inf], [1.5,np.inf]]))
    #sy = abs(VR + pb20[0] - vr)    
    #for _ in range(1):
    pb2, covpb = curve_fit(lin, xdata=VR, ydata=vr, p0=np.array([1,0])) #, sigma=sy)#, bounds=np.array([[0.6, -np.inf], [1.5,np.inf]]))
    #    sy = abs(VR*pb2[0] + pb2[1] - vr)
    #s = pb2[0]
    #i = pb2[1]
    X = [[a] for a in VR]
    y = vr
    reg = HuberRegressor().fit(X, y)
    s2, i2 = reg.coef_, reg.intercept_
    
    s_err = np.sqrt(covpb[0,0])
    i_err = np.sqrt(covpb[1,1])
    a, b = np.min(cat_V - cat_R), np.max(cat_V-cat_R)
    plt.figure(dpi=300, figsize=(10,10))
    #plt.title(gal)
    plt.plot(VR, vr, 'ob')
    #for eb, cb, k in zip(VR, vr, inds0):
    #    plt.plot(eb, cb, 'ob') #, label='%s %s %s' %(k, np.round(eb,2), np.round(cb,2)))
    plt.plot([a,b], [a*s2+i2, b*s2+i2], '-r', label=r'$c_2 = %5.3f \pm %5.3f$ $c_{vr} = %5.3f \pm %5.3f$' %(s2, s_err, i2, i_err))
    plt.ylabel('('+filters[1] +' - '+ filters[2]+') inst', fontsize=15)
    plt.xlabel('(' + filters[1] + '-' + filters[2] + ') catalog', fontsize=15)
    #plt.xlabel(filters[1]+' - '+filters[2])
    #plt.ylabel('v-r')
    plt.legend(fontsize=15)
    #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    #plt.xlim(-2, 2)
    #plt.ylim(-2, 2)
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.tight_layout()
    plt.savefig(Path(out_dir, gal+'%s%s.png' %(filters[1], filters[2])))
    plt.clf()
    plt.close()
    #plt.show()

    #pb00, covpb = curve_fit(lin00, xdata=bv, ydata=Vv, p0=np.array([0])) #, bounds=np.array([[0.6, -np.inf], [1.5,np.inf]]))
    #sy = abs(pb00[0] - Vv) 
    
    #for _ in range(1):
    pb0, covpb = curve_fit(lin, xdata=bv, ydata=Vv, p0=np.array([0,0]))#, sigma=sy)
    #    sy = abs(pb0[0]*bv*0 + pb0[1] - Vv) 
    X = [[a] for a in bv]
    y = Vv
    reg = HuberRegressor().fit(X, y)
    s0, i0 = reg.coef_, reg.intercept_

    #s = pb0[0]
    #i = pb0[1]
    s_err = np.sqrt(covpb[0,0])
    i_err = np.sqrt(covpb[1,1])
    ia = np.argmin(bv)
    ib = np.argmax(bv)
    
    BV = est_B - est_V
    _BV = np.array([])
    _V = np.array([]) 
    for b, a in sorted(zip(est_V, BV)):
        _BV = np.append(_BV, a)
        _V = np.append(_V, b)

    
    plt.figure(dpi=300, figsize=(10,10))
    #plt.title(gal)
    plt.plot(bv, Vv,'ob')
    #for eb, cb, k in zip(bv, Vv, inds0):
    #    plt.plot(eb, cb, 'ob')#, label='%s %s %s' %(k, np.round(eb,2), np.round(cb,2)))
    plt.plot(_BV, _BV*s0 + i0, '-r', label=r'$c_0 = %5.3f \pm %5.3f$ $c_v = %5.3f \pm %5.3f$' %(s0, s_err, i0, i_err))
    plt.xlabel('(' + filters[0] + '-' + filters[1]+ ') catalog', fontsize=15)
    plt.ylabel(filters[1]+ ' catalog' +' - ' +  filters[1] + ' inst', fontsize=15)
    plt.legend(fontsize=15)
    #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.tight_layout()
    plt.savefig(Path(out_dir, gal+'%s.png' %filters[1]))
    plt.clf()
    plt.close()
    #plt.show()

    #plt.figure()
    #plt.plot(cat_B - cat_V, cat_V - cat_R, 'o')
    #plt.show()
    
    #move_phot_sys(fnameB, fnameV, fnameR, pb0, pb1, pb2, out_dir)
    move_phot_sys(fnameB, fnameV, fnameR, [s0, i0], [s1, i1], [s2, i2], out_dir)
#----------------------------------------------------------------------------------------------------------------------------------------


def get_equals_solo(m_calib, m_inst, fname, filt, out_dir, inds0, out_file):

    gal = fname.stem 
   
    inds = np.where(~np.isnan(m_inst))
    m_calib = m_calib[inds]
    m_inst = m_inst[inds]

    m_inst0, m_calib0 = clean(m_inst, m_calib)
    
    
    #def model(x, *p):
    #    return x -  p[0]
    #
    #pb0, covpb = curve_fit(model, xdata=m_calib - m_inst, ydata=0, p0=np.array([30.0])) #, bounds=np.array([[0.6, -np.inf], [1.5,np.inf]]))
    

    #X = [[a] for a in est_B]
    #y = m_calib
    #reg = HuberRegressor().fit(X, y)
    #s, i = reg.coef_, reg.intercept_
        
    #s = pb0[0]
    #i = pb0[1]
    #s_err = np.sqrt(covpb[0,0])
    #i_err = np.sqrt(covpb[0,0])
    
    m0, m0_std = fit(m_inst0, m_calib0)

    a, b = np.min(m_calib), np.max(m_calib)
    #a, b = np.min(m_inst), np.max(m_inst)
    
    plt.figure(dpi=300, figsize=(10,10))
    plt.title(gal)
    plt.plot(m_calib, m_calib-m_inst, 'bo', label='initial')#, label='%s %s %s' %(k, np.round(eb,2), np.round(cb,2)))
    plt.plot(m_calib0, m_calib0-m_inst0, 'go', label='cleaned')#, label='%s %s %s' %(k, np.round(eb,2), np.round(cb,2)))
    plt.plot([a,b], [m0, m0], '-r' , label=r'$m_0 = %5.3f \pm %5.3f$' %(m0, m0_std))
    #plt.plot([a,b], [m0-3*m0_std, m0-3*m0_std], '--r' )
    #plt.plot([a,b], [m0+3*m0_std, m0+3*m0_std], '--r' )
    plt.ylabel('$'+filt+'_{calib}'+'-'+filt+'_{instr}'+'$', fontsize=15)
    plt.xlabel('$' + filt+'_{calib}'+ '$', fontsize=15)
    plt.legend(fontsize=15)
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.tight_layout()
    plt.savefig(Path(out_dir, gal+'.png'))
    plt.clf()
    plt.close()
    
    move_phot_sys_solo(fname, float(m0), float(m0_std), out_dir, out_file)
#--------------------------------------------------------------------------------------------------------------------------------


def move_phot_sys_solo(fnameB, m0, m0_std, out_dir, out_file):

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')

        hdu = fits.open(fnameB)
        header = hdu[0].header
        data = hdu[0].data

        header['m0'] = m0
        header['m0_std'] = m0_std
        #c0  = pb0[0]
        #cv  = pb0[1]


        print('m0', m0)
        print('m0_std', m0_std)

        fits.writeto(Path(out_dir, out_file), data, header=header, overwrite=True) 
    return
