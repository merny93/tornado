import numpy as np
import time

##define a numerical derivative method
def num_deriv(fun,diff_func,pars,dpar):
    #calculate numerical derivatives of 
    #a function for use in e.g. Newton's method or LM
    # derivs=np.zeros([len(x),len(pars)])
    derivs = []
    for i in range(len(pars)):
        pars2=pars.copy()
        pars2[i]=pars2[i]+dpar[i]
        f_right=fun(pars2)
        pars2[i]=pars[i]-dpar[i]
        f_left=fun(pars2)
        derivs.append(diff_func(f_right,f_left)/(2*dpar[i]))
    return np.array(derivs).T

def newton_solve(model_func, resid_func, noise_func, chi_func, diff_func, init_params):

    ##get initial model
    pars = np.array(init_params)
    print(pars)
    init_model = model_func(pars)
    noise_model = noise_func(init_model)
    Ninv= np.diag(1/noise_model) #noise matrix
    
    dpar = pars*1e-4
    dpar = np.maximum(np.ones(dpar.size)*1e-4, dpar)
    # print(dpar)


    print("Now we do a Newthons method solver")

    #initial guess
    chi_org = chi_func(model_func(pars))
    print("original chi", chi_org)
    ##lets keep track of the derivatives:
    derivatives = []
    parameters = []
    for i in range(25):
        #get val where we are
        model=model_func(pars)

        #linearize around that point and save the derivatives 
        derivs=num_deriv(model_func,diff_func,pars,dpar)
        derivatives.append(derivs)

        
        #how far from the truth
        resid= resid_func(model)
        #solve the linear model

        print("sizes", resid.shape, derivs.shape, Ninv.shape)
        lhs=derivs.T@Ninv@derivs #this is curvature
        rhs=derivs.T@Ninv@resid 
        lhs_inv=np.linalg.pinv(lhs)
        step=lhs_inv@rhs #project solution and move over
        # print(step)
        pars=pars+step
        parameters.append(pars)

        ##calculate the next derivative step as 10% of the error
        if (i>3):
            dpar = np.sqrt(np.diag(lhs_inv))*0.1
            dpar = np.maximum(np.ones(dpar.size)*1e-4, np.nan_to_num(dpar))
        # 

        print("On the: " , i+1, " itteration and params are: ", pars)
        #this next step is not really needed but its nice to see how chi^2 changes
        chi_now = chi_func(model_func(pars))
        print("Chi^2 is now: ", chi_now)
        if np.abs(chi_now - chi_org) < 0.1:
            print("we have converged")
            print("with errors given by :", np.sqrt(np.diag(lhs_inv)))
            break
        chi_org = chi_now
    return pars, np.sqrt(np.diag(lhs_inv))