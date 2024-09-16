#----------------------------------------------------------------------------#
# Created on Mon Dec 14 13:12:33 2020
# Last edit September 14, 2024
# James Thompson - University of Texas at Austin, University of Pittsburgh
# For use by the University of West Indies
#----------------------------------------------------------------------------#

# perform T-E separation using TES
#C.G. Hughes, 2012 Caltech/JPL
#

def NEM_planck(emax, surfradin, wave_tes, it):
    #NEM - Normalized Emissivity Method
    #
    #ASSUMPTIONS:
    #wave_tes has the same number of cells as surfradin and skyradin
    #
    #DEPENDENCIES
    #numpy for zeros_like
    #math for pi (mmm, pi.....)
    from numpy import log, exp, zeros_like
    from math import pi


    #define planck constants - nb, c1 includes a pi term. (2hc^2 is 1.19e-16)
    c1=3.7418e-22 # W.m-2
    c2=0.014388 # um K

    #estimate the ground-emitted radiance and brightness temperature
    R=surfradin#-(1-emax)*skyradin/pi
    T=(c2/wave_tes)*(log((emax*c1)/(pi*R*wave_tes**5)+1))**-1

    #Tnem is the maximum brightness temperature
    Tnem=T.max()

    #compute blackbody radiance using Tnem
    B=c1/(wave_tes**5*pi*(exp(c2/(wave_tes*Tnem))-1))
    #compute emissivity from blackbody radiance
    e=R/B

    #iteratively correct for the downwelling radiance
    Rold=R
    QA=0
    Re=R
    k=0
    # for k in range(it):
    #     if k==it-1:
    #         QA=1
    #         break
    #     Re=surfradin-(1-e)*skyradin/pi
    #     Te=(c2/wave_tes)*(log((c1*e)/(pi*Re*wave_tes**5)+1))**-1
    #     diff=Re-Rold

    #     #compute the delta of corrected radiances - for convergence to
    #     #occur, it must occur in all bands

    #     #arrays of T/F of length wave_tes
    #     dcon=(abs(diff)<0.05)
    #     ddiv=(diff>0.05)

    #     #after 2 iterations, check for convergence / divergence
    #     if (dcon.sum() == len(wave_tes)) and k>1:
    #         QA=0
    #         break

    #     if (ddiv.sum() == len(wave_tes)) and k>1:
    #         QA=1
    #         break

    #     #re-assign values for next iteration
    #     Tnem=Te.max()
    #     B=c1/(wave_tes**5*pi*(exp(c2/(wave_tes*Tnem))-1))
    #     e=Re/B
    #     Rold=Re

    return [e, Tnem, Re, QA, k]

def init_TES(surfrad):
    #given surface radiance and sky irradiance, split them up by
    #TES vs non-TES bands for passing to calculation later
    #ASSUMPTIONS:
    #TIR bands are 41-49 (count from 0) - ABI TIR 10-16 
    #TES bands are 42,43,46,47,48 (count from 0) - ABI TES bands are -  
    #surfrad and skyr are bands x lats x lons, with bands=TIR bands
    #
    #DEPENDENCIES:
    #numpy, for zeros
    from numpy import zeros

    #initialize the 4 arrays
    # skyr_tes=zeros((3,skyr.shape[1],skyr.shape[2]),dtype=float)
    # skyr_nontes=zeros((4,skyr.shape[1],skyr.shape[2]),dtype=float)
    surfrad_tes=zeros((3,surfrad.shape[1],surfrad.shape[2]),dtype=float)
    surfrad_nontes=zeros((4,surfrad.shape[1],surfrad.shape[2]),dtype=float)

    #this method of assigning them is inelegant
    #go through the sky irradiance first, and assign tes or nontes
    # skyr_nontes[0,:,:]=skyr[0,:,:]
    # skyr_tes[0,:,:]=skyr[1,:,:]
    # skyr_nontes[1,:,:]=skyr[2,:,:]
    # skyr_tes[1,:,:]=skyr[3,:,:]
    # skyr_tes[2,:,:]=skyr[4,:,:]
    # skyr_nontes[2,:,:]=skyr[5,:,:]
    # skyr_nontes[3,:,:]=skyr[6,:,:]


    #do the same thing now for surface radiance
    surfrad_nontes[0,:,:]=surfrad[0,:,:]
    surfrad_tes[0,:,:]=surfrad[1,:,:]
    surfrad_nontes[1,:,:]=surfrad[2,:,:]
    surfrad_tes[1,:,:]=surfrad[3,:,:]
    surfrad_tes[2,:,:]=surfrad[4,:,:]
    surfrad_nontes[2,:,:]=surfrad[5,:,:]
    surfrad_nontes[3,:,:]=surfrad[6,:,:]


    return [surfrad_tes, surfrad_nontes]

def TES_for_ABI(surfrad_tes, surfrad_nontes, ecw_tes, ecw_tir):
    #given TIR bands, perform T-E separation on them using the TES split
    #window technique
    #
    #ASSUMPTIONS:
    #surfrad_tes through skyr_nontes are all the same size
    #ecw_tes and ecw_tir have band centers in microns
    #
    #DEPENDENCIES
    #numpy - bunches of functions
    #I started with one or two from numpy. I should have just imported numpy.
    import numpy as np
    from numpy import alltrue, any, arange, array_equal, asarray, diff, exp
    from numpy import iscomplex, isnan, log, polyder, polyfit, polyval
    from numpy import setdiff1d, sometrue, unique, var, where, zeros
    from numpy import zeros_like, max, pi

    #wavelengths for this scene in microns
    wave_tes=asarray(ecw_tes)
    wave_nontes=asarray([ecw_tir[0], ecw_tir[2], ecw_tir[5], ecw_tir[6]])

    #set up coefficients and reform the radiance arrays as needed for
    #bad data. also edit set wave_tes and wave_nontes appropriately
    #identify bands that had -999s in them - look for min() below 0.
    #we cannot use these bands. Not applicaable for GOES

    #most common case - check first!
    #bands - 43,44,47,48,49: 0.9921, 0.74329, 0.78522 (all TES bands)
    #coeff=asarray([0.9921, 0.74329, 0.78522])
    #set coefficients
    #bands - 43, 47, 48: 0.98304, 0.73144, 0.74818 ************************************************
    coeff=asarray([0.98304, 0.73144, 0.74818])
    
    #convert wavelengths to meters
    wave_tes=wave_tes*1e-6
    wave_nontes=wave_nontes*1e-6

    #TES algorithm: max emissivity and Planck coefficients
    emax=0.99
    c1=3.7418e-22 # W.m-2
    c2=0.014388 # um K

    #number of iterations for removing reflected downwelling radiance
    it=13

    #initialize my emissivities - all are bands x lats x lons
    emisf=zeros_like(surfrad_tes)
    emisf_nontes=zeros_like(surfrad_nontes)
    emis_NEM=zeros_like(surfrad_tes)

    #initialize surface temp array - lats x lons
    Ts=zeros((surfrad_tes.shape[1], surfrad_tes.shape[2]))
    QAmap=zeros((surfrad_tes.shape[1], surfrad_tes.shape[2]))

    #loop over all pixels
    for i in range(surfrad_tes.shape[1]):
        # if i==0 or ((i+1) % 100)==0:
        #     print("starting TES on row", i+1)

        for j in range(surfrad_tes.shape[2]):
            #check for bad data (NaNs)
            if any(isnan(surfrad_tes[:,i,j])):
                continue

            # #ok, no nans, let's see if all 0s instead
            # surftest=surfrad_tes[:,i,j]
            # if alltrue(surftest == 0):
            #     continue

            #no nans, not all 0s.... proceed

            #NEM module
            surfradin=surfrad_tes[:,i,j]
            #skyradin=skyr_tes[:,i,j]
            surfradin_nontes=surfrad_nontes[:,i,j]
            #skyradin_nontes=skyr_nontes[:,i,j]
            [ef, Tnemf, Re, QA, k]=NEM_planck(emax, surfradin,
                                              wave_tes, it)
            #check if NEM diverged
            if QA==1:
                QAmap[i,j]=1
                continue

            #if radiances or emissivities not real, continue as well
            if sometrue(iscomplex(Re)) or ef.min()<=0:
                continue

            emis_NEM[:,i,j]=ef

            #refinement of Emax
            v=var(ef)

            # if v>1.7e-4:
            #     emax=0.96
            #     #call NEM again
            #     [ef, Tnemf, Reff, QA, k]=NEM_planck(emax, surfradin,
            #                                   wave_tes, it)
            #     if QA==1:
            #         QAmap[i,j]=1
            #         continue

            # else:
            #     #it seems like an emax of 0.99 is already done, so use
            #     #those results rather than calling it again
            #     emaxr=asarray([0.92, 0.95, 0.97, 0.99])
            #     varlist=zeros_like(emaxr)

            #     varlist[3]=v

            #     for val in range(3):
            #         [e, Tnem, Reff, QA, k]=NEM_planck(emaxr[val], surfradin,
            #                                           wave_tes, it)
            #         varlist[val]=var(e)

            #     #fit a parabola to varlist vs emaxr
            #     p=polyfit(emaxr, varlist, 2)
            #     stepval=0.001
            #     x=arange(emaxr[0],emaxr[3]+stepval, stepval)
            #     f=polyval(p,x)

            #     #find Emax(V_min)
            #     xpos=f.argmin()
            #     emax_min=x[xpos]

            #     #Test 1: Find the average slope of the fit
            #     #this is different from Glynn's code - he called the
            #     #MATLAB function gradient with f and x. However, there
            #     #doesn't seem to be a good equivalent in python,so I'll use
            #     #diff to calculate out what I think is the same thing -
            #     #diff just gives the difference between the current val and
            #     #the previous value, or the delta. So, I'm effectively doing
            #     #delta-y / delta-x to calculate slope, and then looking at the
            #     #absolute value of the mean
            #     FX= diff(f) / diff (x)
            #     d1=abs(FX.mean())

            #     #Test 2: Find the 2nd derivitive (minima of fit)
            #     d2=polyder(p, 2)[0]


            #     #if we pass test 1 and test 2, rerun with my the new Emax value
            #     if (d1 < 1.0e-3) and (d2 > 1e-4):
            #         #we have a new acceptable emax calculated
            #         emax=emax_min
            #     else:
            #         #"in all these cases, if the test is failed Emax = 0.983
            #         # is assumed"
            #         emax=0.983

            #     #now that we have a new emax, recall the NEM module
            #     [ef, Tnemf, Reff, QA, k]=NEM_planck(emax, surfradin,
            #                                         wave_tes, it)

            #     if QA==1:
            #         QAmap[i,j]=1
            #         continue
            
            Reff = Re

            #Ratio (RAT) Module

            bm = np.mean(ef);
            # bm = np.max(ef);

            beta=zeros_like(wave_tes)
            beta=ef*bm #emissity divided by the avg

            #MMD Module

            #find the spectral contrast
            MMD=beta.max() - beta.min()

            #calculate Emin - use MASTER coefficients rather than the
            #ASTER ones in the ATBD
            #

            emin=coeff[0]-(coeff[1]*(MMD**coeff[2]))

            #given Emin, calculate our emissivities
            emis=beta*(emin/beta.min())

            #ensure that all emis values are real; if they are not (ie, if
            #beta.min() was 0 somehow...), replace them with emissivity 0.00
            for b in range(wave_tes.shape[0]):
                if not iscomplex(emis[b]):
                    emisf[b,i,j]=emis[b]
                else:
                    emisf[b,i,j]=0.00


            #Final correction for Ts
            bmax=emis.argmax()
            Tf=(c2/wave_tes[bmax])*(log((c1*emis[bmax])/
                                       (pi*Reff[bmax]*wave_tes[bmax]**5)+1))**-1
            if not iscomplex(Tf):
                Ts[i,j]=Tf
            else:
                Ts[i,j]=0.00

            #Calculate emissivity for non-TES bands using retrieved Ts
            if wave_nontes.size > 0:
                #estimate ground-emitted radiance
                R=surfradin_nontes# - (1-emax)*skyradin_nontes/pi

                #compute blackbody radiance using Tf
                B=c1/(wave_nontes**5*pi*(exp(c2/(wave_nontes*Tf))-1))

                #compute emissivity
                emisf_nontes[:,i,j]=R/B
            else:
                emisf_nontes[:,i,j]=0.00

    return [emisf, emisf_nontes, Ts, QAmap, wave_tes, wave_nontes]

def shiftdim(x, n=None):
    r""" Matlab's shiftdim in python.

    Examples
    --------
    >>> import oceans.ff_tools as ff
    >>> a = np.random.rand(1,1,3,1,2)
    >>> print("a shape and dimension: %s, %s" % (a.shape, a.ndim))
    a shape and dimension: (1, 1, 3, 1, 2), 5
    >>> # print(range(a.ndim))
    >>> # print(np.roll(range(a.ndim), -2))
    >>> # print(a.transpose(np.roll(range(a.ndim), -2)))
    >>> b = ff.shiftdim(a)
    >>> print("b shape and dimension: %s, %s" % (b.shape, b.ndim))
    b shape and dimension: (3, 1, 2), 3
    >>> c = ff.shiftdim(b, -2)
    >>> c.shape == a.shape
    True

    Notes
    -----
    http://www.python-it.org/forum/index.php?topic=4688.0
    """

    import numpy as np

    def no_leading_ones(shape):
        shape = np.atleast_1d(shape)
        if shape[0] == 1:
            shape = shape[1:]
            return no_leading_ones(shape)
        else:
            return shape

    if n is None:
        # returns the array B with the same number of
        # elements as X but with any leading singleton
        # dimensions removed.
        return x.reshape(no_leading_ones(x.shape))
    elif n >= 0:
        # When n is positive, shiftdim shifts the dimensions
        # to the left and wraps the n leading dimensions to the end.
        return x.transpose(np.roll(range(x.ndim), -n))
    else:
        # When n is negative, shiftdim shifts the dimensions
        # to the right and pads with singletons.
        return x.reshape((1,) * -n + x.shape)

def MASTER_TLR(surfradtlr,gpMask):
    # Thermal Log Residuals to refine gray-pixel map
    import numpy as np
    import numpy.ma as ma
    import bottleneck as bn

    szt = surfradtlr.shape
    wave = [8.5969,9.0396,10.6386,11.3522,12.1683]

    # Term 4
    t3s = 0
    for i in range(len(wave)):
        t3s += np.nansum( wave[i]*surfradtlr[i,:,:] )

    tlr4 = t3s/( szt[1]*szt[2]*len(wave))
    TLR = np.zeros((len(wave),surfradtlr.shape[1],surfradtlr.shape[2]),dtype=np.float64)

    for b in range(len(wave)):
        # Term 1
        tlr1 = wave[b] * np.log(surfradtlr[b])
        # Term 2
        tlr2 = 0
        for w in range(len(wave)):
            tlr2 += wave[w] * np.log(surfradtlr[w])
        # Term 3
        tlr3 = wave[b] * bn.nanmean(np.log(surfradtlr[b]))
        TLR[b] = tlr1 - tlr2 - tlr3 + tlr4

    # Make masked copies
    # TLR_ng = ma.masked_array(TLR, mask= ( [~gpMask,~gpMask,~gpMask,~gpMask,~gpMask]) )
    # TLR_g  = ma.masked_array(TLR, mask=   [gpMask,gpMask,gpMask,gpMask,gpMask]  )

    #create a nan array
    NaNs = np.empty(TLR.shape)
    NaNs.fill(np.nan)

    multiBandMask = np.array([gpMask.transpose(),gpMask.transpose(),gpMask.transpose(),gpMask.transpose(),gpMask.transpose()])

    # Make NAN copies
    TLR_ng = np.where(multiBandMask, TLR, NaNs)
    TLR_g  = np.where(multiBandMask, NaNs, TLR)

    #we want the nanmean of both the y and z dims (assuming we have a (x,y,z) matrix) so do nanmean twice
    TLR_gm = bn.nanmean( bn.nanmean(TLR_g, axis = 1), axis = 1)

    return [TLR_ng,TLR_gm]

def MASTER_Tg(Ywvs,PWVs,c_wvs,wave):
    # Thermal Log Residuals to refine gray-pixel map
    import numpy as np

    # Constants
    c1 = 3.741775e-22
    c2 = 0.0143877

    # Change negative values to nans
    Ywvs[Ywvs < 0] = np.nan

    # Compute brightness temperature (old way) Wavelengths
    TB = np.zeros((wave.shape[0],Ywvs.shape[1],Ywvs.shape[2]),dtype=float)
    for b in range(TB.shape[0]):
        TB[b] = c2 / ( wave[b] * np.log(c1 / (wave[b]**5*np.pi*Ywvs[b]) + 1) )

    # Compute T_gis
    Tg = np.zeros((wave.shape[0],Ywvs.shape[1],Ywvs.shape[2]),dtype=float)
    for n in range(Tg.shape[0]):

        Tgo = c_wvs[n,0] + c_wvs[n,1]*PWVs + c_wvs[n,2]*PWVs**2

        Tgsum = 0
        c=3
        for k in range(wave.shape[0]):
            Tgsum = Tgsum + ( c_wvs[n,c] + c_wvs[n,c+1] * PWVs + c_wvs[n,c+2]*PWVs**2) * TB[k]
            c+=3

        Tg[n] = Tgo + Tgsum

    return Tg

def master_ginterp_v2(g,gp,sc,r2):
    import numpy as np
    from scipy.misc import imresize

    # Displacement vector
    dispVector = np.array([ [ 1, 0],
                            [-1, 0],
                            [ 1, 1],
                            [ 0, 1],
                            [-1, 1],
                            [ 1,-1],
                            [ 0,-1],
                            [-1,-1]])

    # Set negative g values to 1 in all bands
    # also set non-gray, clear pixels to 1 using gray pixel thresholder, gp
    for badPixel in np.argwhere( (g > 3) | (g < 0.1) | np.isnan(g) | gp):
        g[:,badPixel[1],badPixel[2]] = np.array([ 1,  1,  1,  1,  1])

    # print("Num grey points:",np.argwhere( (g < 3) & (g > 0.1)).shape[0])

    cgp = np.argwhere( g[0,:,:] > 1 )

    if cgp.shape[0] == 0: #if the argwhere returned no results
        gtemp = np.zeros(g.shape,dtype=float)
        gorig = np.copy(gtemp)
    else:
        # resize gammas by 1/sc using nearest neighbor interopolation
        gSize = ( int(np.floor(g.shape[1]/sc)) , int(np.floor(g.shape[2]/sc)))
        gtemp = np.zeros((g.shape[0],gSize[0],gSize[1]),dtype=float)

        for b in range(g.shape[0]):
            gtemp[b] = imresize(g[b,:,:], gSize, interp='nearest', mode='F')

        # Remove gamma values with 2 or less neighbors (e.g. noise)
        gtemp = np.where( gtemp == 1, np.nan, gtemp )

        # [c1 c2] = find(gtemp{1}>0);
        c = np.argwhere(gtemp[0,:,:]>0 )

        for j in range(c.shape[0]):
            neighbors = dispVector + np.tile(c[j],(8,1))
            #remove edge rows & columns below 0
            neighbors = neighbors[~(neighbors < 0).any(1)]
            #then ros and columns above their max values
            neighbors = neighbors[(neighbors[:,1] < gtemp.shape[2])]
            neighbors = neighbors[(neighbors[:,0] < gtemp.shape[1])]

            # get neighbor values
            nbor = gtemp[0,neighbors[:,0],neighbors[:,1]]

            # remove pixels with 2 or less neighbors
            if np.argwhere(nbor>0).shape[0] <= 2 :
                gtemp[:,c[j,0],c[j,1]] = 1

        ## Remove outliers using IQR
        # Find interquartile ranges
        Q1  = np.percentile(gtemp[0],25)
        Q3  = np.percentile(gtemp[0],75)
        iq1 = Q1-Q3
        it  = 1.5

        # Find outliers
        ern = ((gtemp[0] > Q1+it*iq1) | (gtemp[0] < Q3-it*iq1))

        # Set outliers in each band to 1
        for b in range(gtemp.shape[0]):
            gtemp[b,ern] = 1

        # add non-gray 1 values back in
        gtemp[np.isnan(gtemp)] = 1
        gorig = np.copy(gtemp)

        eee = -4 # p value determines amount of weight given to nearest pixels
        testLen1 = 0

        while np.argwhere(gtemp[1] == 1).shape[0] > 0 : #while we still have shape matches

            vc = np.reshape(gtemp, (gtemp.shape[0],gtemp.shape[1]*gtemp.shape[2]))
            # remove bare pixels (=1) from vc
            vcNoOnesList = []
            for b in range(gtemp.shape[0]):
                vcNoOnesList.append([i for i in vc[b].tolist() if (i != 1 and not np.isnan(i)) ])
            vc = np.array(vcNoOnesList)

            grayPix = np.argwhere( gtemp[0,:,:] != 1 )
            xg = grayPix[:,0]
            yg = grayPix[:,1]

            if grayPix.shape[0] == 0:
                print("No Gray Pixels Found. Raise Error.")
                raise IndexError

            barePix = np.argwhere(gtemp[0,:,:] == 1)
            x = barePix[:,0]
            y = barePix[:,1]

            testLen = y.shape[0]
            if testLen == testLen1:
                # print("f")
                # then the remaining bare pixels are not withing r2 of gray pixels
                for j in range(y.shape[0]):
                    D = np.sqrt((x[j]-xg)**2 + (y[j]-yg)**2)
                    D2_select = D<r2
                    D2 = D[D2_select]
                    #Increase distance if D2 is empty
                    if D2.shape[0] < 1:
                        r3 = r2
                        while D2.shape[0] < 1:
                            r3 += (r2*2)
                            D2_select = D<r3
                            D2 = D[D2_select]

                        # print("e")

                        # find D2^eee but skip zero values to avoid div-by-zero errors
                        wV = np.where(D2 < 1e-6 , D2, D2**eee)

                        for b in range(g.shape[0]):
                            vcc = (vc[b,:])[D2_select]

                            Vt1 = vcc * wV
                            gtemp[b,x[j],y[j]] = np.sum(Vt1)/np.sum(wV)
                        #end for b in range(g.shape[0]):
                    del vcc, Vt1, wV
                    #end if D2.shape[0] < 1:
                #end for j in range(testLen):
                break
            #end if testLen == testLen1 :

            for j in range(y.shape[0]):
                # print("g")
                D = np.sqrt( (x[j] - xg)**2 + (y[j] - yg)**2 )
                D2_select = D<r2
                D2 = D[D2_select]

                k = 0
                if D2.shape[0] < 1:
                    continue
                else:
                    # print("h")
                    # find D2^eee but skip zero values to avoid div-by-zero errors
                    wV = np.where(D2 < 1e-5 , D2, D2**eee)
                    for b in range(g.shape[0]):
                        vcc = (vc[b,:])[D2_select]
                        # try:
                        #     vcc = (vc[b,:])[D2_select]
                        # except IndexError as e:
                        #     print()
                        #     print("vc[b,:].shape",vc.shape)
                        #     print("D2_select.shape",D2_select.shape)
                        #     print()
                        #     raise e

                        Vt1 = vcc * wV
                        gtemp[b,x[j],y[j]] = np.sum(Vt1)/np.sum(wV)

                    #end for b in range(g.shape[0]):
                del vcc, Vt1, wV
                #end if D2.shape[0] < 1:
            #end for j in range(y.shape[0]):

            #set our vars up for the next part of the loop
            testLen1 = testLen
        #end while np.argwhere(gtemp[1] == 1).shape[0] > 0 : #while we still have shape matches
    #end if cgp.shape[0] == 0

    return [gtemp, gorig]
