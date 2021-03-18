import numpy as np 
import spify
import x as x_py
from scipy.optimize import curve_fit

if __name__=="__main__":
    # Load the angles
    data = np.load('angles_tin.npz')
    angles = np.array(data['angles'])
    errors = np.array(data['uncertainty'])
    one_d2 = 4*np.sin(np.deg2rad(np.array(angles/2)))**2/1.541838e-10**2

    from itertools import product
    # create a list of all the possible hkl
    coms = list(product(range(5), repeat=3))


    def check_good(x):
        '''
        Checks if the combination of x (hkl) is plausible
        using rules from bilbao crystallographic server
        for No.141 4b (tin)
        '''
        good = False
        [h,k,l] = x
        if sum(x) ==0:
            good =False
        elif (k == 0) and (l ==0) and (h%2 == 0):
            good = True 
        elif (h==0) and (k==0) and (l%4 ==0):
            good = True
        elif (h==k) and ((2*h+l)%4 == 0):
            good  = True 
        elif (h==0) and ((k+l)%2==0):
            good = True 
        elif (h%2 ==0) and (k%2==0) and (l==0):
            good = True 
        elif ((h+k+l)%2 ==0) and (((l-1)%2==0) or ((2*h+1)%4 ==0)):
            good = True
        else:
            pass
    
        return good

    # Refine the list using the rules for tin
    good = list(filter(check_good, coms))

    # Function that orders the choices properly using a and c from litterature
    a,c = 5.83e-10, 3.18e-10
    norm = lambda x: (x[0]**2+x[1]**2)/a**2 + (x[2]**2)/c**2
    # Sorts the good list using the norm function
    good.sort(key = norm)

    from functools import reduce
    # Removes degeneracy in h,k
    good = list(reduce(lambda x,y: x + [y] if norm(x[-1])!= norm(y) else x[:-1] + [max(x[-1],y,key=lambda x: x[0])], good, [good[0]]))
    # print(good)
    # print(norm(np.array(good).T), one_d2)
    # print(good)

    def norm_litt(y):
        for element in one_d2:
            if np.abs(norm(y)-element)<0.1e19:
                return True 
        return False

    final_filter = np.array(list((filter(norm_litt, good))))

    x_fit = [final_filter[:,0],final_filter[:,1],final_filter[:,2],np.ones(len(angles))*0.1541838e-9]
    # print(x_fit)
    popt,pcov = curve_fit(x_py.tin_func, x_fit, angles, p0= [a,c,0], sigma=errors)
    unc = np.sqrt(np.diag(pcov))
    chi_sqd = np.sum((angles - x_py.tin_func(x_fit, *popt))**2/errors**2)/len(angles-2)
    # print(chi_sqd)

    x_plot = 1/4*((x_fit[0]**2+x_fit[1]**2)/popt[0]**2 + x_fit[2]**2/popt[1]**2)*x_fit[3]**2
    spify.residual_plot(x_plot, angles, errors, x_py.tin_func_2, [popt[2]], 
                    r"$\lambda^2(h^2+k^2+l^2) \ (10^{-19} \ \mathrm{m}^2)$",
                        r"$2\theta \ (^{\circ} )$",
                        "Residuals", 
                        'tin' + "_lattice", renorm=False)
    print(popt, unc)