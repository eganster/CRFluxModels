from CRFluxModels import *


class MascarettiBlasiEvoli(PrimaryFlux):
    """
    Model of the primary cosmic ray spectrum based on a paper preprint from C. Mascaretti, P. Blasi and C. Evoli.
    
    Two models, exponential-square cutoff and spectral break were fitted to data from ARGO and KASCADE.
    Each model features 3 componments: protons, Helium and a light extragalactic (carbon) component.
    
    # p/H - 14
    # He - 402
    # light, extragalactic component (proton only) - 14
    
    Args:
      model (str): on of ["cutoff_ARGO", "break_ARGO", "cutoff_KASCADE", "break_KASCADE"]
    """

    def __init__(self, model="cutoff_ARGO"):

        self.name = "MascarettiBlasiEvoli"
        self.model, self.fit = model.split("_")
        self.params = {}

        # check for model and fit string to be valid
        if self.model not in ["cutoff", "break"] or self.fit not in ["ARGO", "KASCADE"]:
            raise Exception(
                    'MascarettiBlasiEvoli(): Unknown model/fit version requested: '+'{}+{}'.format(self.model, self.fit)
                   +'Allowed: model = [\'cutoff\', \'break\'] and fit = [\'ARGO\', \'KASCADE\']')
        
        # set params
        # params[corsika_id] = (normalization @ 10 TeV, slope)
        
        # p and He always the same
        self.params[14]   = (1.50, 2.71) # p
        self.params[402]  = (1.50, 2.63) # He original was (1.50, 2.64)
        self.nucleus_ids = self.params.keys()
        
        # rigiditys different
        self.rigidity = {"cutoff" : {"ARGO"    : 1.30e15,
                                     "KASCADE" : 15.1e15},
                         "break"  : {"ARGO"    : 640e12,
                                     "KASCADE" : 5.8e15}}
        self.delta = 1.0/3.0
        
        self.calculate_params()
        
    # calculate b_i, E_lim and set it (params = (a, \gamma) b, E_lim))
    def calculate_params(self):
        for corsika_id in self.nucleus_ids:
            # get Z, A
            Z, A = self.Z_A(corsika_id)
            E_lim = Z * self.rigidity[self.model][self.fit] / 1e9
            # b = (E/10TeV)**(2 - delta)
            b = self.params[corsika_id][0] * (E_lim/10e3)**(2.0 - self.delta)
            # set it
            self.params[corsika_id] =  self.params[corsika_id][:2] + (b, E_lim)
    
        
    # exponential-square cutoff
    def cutoff(self, corsika_id, E):
        # get Z, A
        Z, A = self.Z_A(corsika_id)
        # a * (E/10TeV)**-gamma * exp(-(E/ZR)**2) * 10**-7
        flux = self.params[corsika_id][0] * (E/10e3)**(-self.params[corsika_id][1]) * np.exp(-1.0*(E/(Z*self.rigidity[self.model][self.fit]/1e9))**2) * 1e-7
        return flux
    
    # change of slope          
    def slope_break(self, corsika_id, E):
        # check for case
        if np.isscalar(E):
            if E <= self.params[corsika_id][3]:
                # a * (E/10TeV)**(-gamma + delta -2) * 20**-7
                flux = self.params[corsika_id][0] * (E/10e3)**(-self.params[corsika_id][1]) * 1e-7
            else:
                # b * (E/10TeV)**(-gamma + delta -2) * 20**-7
                flux = self.params[corsika_id][2] * (E/10e3)**(-self.params[corsika_id][1] + 1.0/3.0 - 2.0) * 1e-7
        else:
            flux = np.where(E <= self.params[corsika_id][3], self.params[corsika_id][0] * (E/10e3)**(-self.params[corsika_id][1]) * 1e-7,
                                                             self.params[corsika_id][2] * (E/10e3)**(-self.params[corsika_id][1] + self.delta - 2.0) * 1e-7)
        
        return flux
        
    # calculate nucleus flux
    def nucleus_flux(self, corsika_id, E):
        # get nearest corsika_id
        corsika_id = self._find_nearby_id(corsika_id)
        
        # check for cutoff model
        if self.model == "cutoff":
            flux = self.cutoff(corsika_id, E)
            if corsika_id == 14:
                flux += 6.00 * (E/100e6)**(-2.70) * 1e-19
        # check for break model
        elif self.model == "break":
            flux = self.slope_break(corsika_id, E)
            if corsika_id == 14:
                flux += 5.00 * (E/100e6)**(-2.70) * 1e-19
            
        return flux




# knee as exponential exp(-(E/ZR)^2) from ARGO
class exp2_argo(PrimaryFlux):
    """data taken from from https://lpsc.in2p3.fr/cosmic-rays-db/#
        """
    
    def __init__(self, model=None):
        self.name = 'exp2_argo'
        self.sname = 'e2ar'
        self.params = {}
        # dictionary[corsika_id] = (normalization @ 10 TeV, slope, charge)
        #self.params[14] = (1.5e-7+0.2e-7,-2.71+0.04,1)   # H
        #self.params[402] = (1.5e-7+0.01e-7,-2.63+0.03,2)  # He
        self.params[14] = (1.5e-7,-2.71,1)   # H
        self.params[402] = (1.5e-7,-2.63,2)  # He
        self.nucleus_ids = self.params.keys()
    
    def nucleus_flux(self, corsika_id, E):
        corsika_id = self._find_nearby_id(corsika_id)
        
        return self.exp2_ar(corsika_id, E)
    
    def exp2_ar(self, corsika_id, E):
        rigidity = 1.3e06
        param = self.params[corsika_id]
        add = 0
        if (corsika_id==14):
            #add = (6.0e-19+0.2e-19)*((E/1e8)**(-2.7))
            add = (6.0e-19)*((E/1e8)**(-2.7))
        Z = param[2]
        return param[0]*((E/1e4)**param[1])*np.exp(-pow(E/(Z*rigidity),2))+add

# knee as exponential exp(-(E/ZR)^2) from KASCADE-Grande
class exp2_kg(PrimaryFlux):
    """data taken from from https://lpsc.in2p3.fr/cosmic-rays-db/#
        """
    
    def __init__(self, model=None):
        self.name = 'exp2_kg'
        self.sname = 'e2kg'
        self.params = {}
        # dictionary[corsika_id] = (normalization @ 10 TeV, slope, charge)
        #self.params[14] = (1.5e-7+0.2e-7,-2.71+0.04,1)   # H
        #self.params[402] = (1.5e-7+0.01e-7,-2.63+0.03,2)  # He
        self.params[14] = (1.5e-7,-2.71,1)   # H
        self.params[402] = (1.5e-7,-2.63,2)  # He
        self.nucleus_ids = self.params.keys()
    
    def nucleus_flux(self, corsika_id, E):
        corsika_id = self._find_nearby_id(corsika_id)
        
        return self.exp2_kagr(corsika_id, E)
    
    def exp2_kagr(self, corsika_id, E):
        rigidity = 15.1e06
        param = self.params[corsika_id]
        Z = param[2]
        add = 0
        if (corsika_id==14):
            #add = (6.0e-19-0.2e-19)*((E/1e8)**(-2.7))
            add = (6.0e-19)*((E/1e8)**(-2.7))
        Z = param[2]
        return param[0]*((E/1e4)**param[1])*np.exp(-pow(E/(Z*rigidity),2))+add

#change slope at rigidity from fit to ARGO data
class dslope_argo(PrimaryFlux):
    """data taken from from https://lpsc.in2p3.fr/cosmic-rays-db/#
        """
    
    def __init__(self, model=None):
        self.name = 'dslope_argo'
        self.sname = 'DS_ARGO'
        self.params = {}
        # dictionary[corsika_id] = (normalization @ 10 TeV, slope, charge)
        #self.params[14] = (1.5e-7-0.2e-7,-2.71-0.04,1)   # H
        #self.params[402] = (1.5e-7-0.01e-7,-2.63-0.03,2)  # He
        self.params[14] = (1.5e-7,-2.71,1)   # H
        self.params[402] = (1.5e-7,-2.63,2)  # He
        self.nucleus_ids = self.params.keys()
    
    def nucleus_flux(self, corsika_id, E):
        corsika_id = self._find_nearby_id(corsika_id)
        
        return self.ds_ar(corsika_id, E)
    
    def ds_ar(self, corsika_id, E):
        param = self.params[corsika_id]
        # rigidity of the knee
        r_knee = 640e3
        # knee in energy
        E0 = param[2]*r_knee
        # add extragalactic light (proton only) component
        add = 0
        ddelta = 1./3
        if (corsika_id==14):
            #add = (6.0e-19+0.2e-19)*((E/1e8)**(-2.7))
            add = 5.0e-19*((E/1e8)**(-2.7))
        if E < E0:
            return param[0]*((E/1e4)**param[1])+add
        else:
            # the normalization has to be changed to grant continuity of the spectrum
            norm = param[0]*((E0/1e4)**(2-ddelta))
            gamm = param[1]+ddelta-2
            return norm*((E/1e4)**gamm)+add

#change slope at rigidity from fit to KASCADE-Grande data
class dslope_kg(PrimaryFlux):
    """data taken from from https://lpsc.in2p3.fr/cosmic-rays-db/#
        """
    
    def __init__(self, model=None):
        self.name = 'dslope_kg'
        self.sname = 'DS_KG'
        self.params = {}
        # dictionary[corsika_id] = (normalization @ 10 TeV, slope, charge)
        #self.params[14] = (1.5e-7-0.2e-7,-2.71-0.04,1)   # H
        #self.params[402] = (1.5e-7-0.01e-7,-2.63-0.03,2)  # He
        self.params[14] = (1.5e-7,-2.71,1)   # H
        self.params[402] = (1.5e-7,-2.63,2)  # He
        self.nucleus_ids = self.params.keys()
    
    def nucleus_flux(self, corsika_id, E):
        corsika_id = self._find_nearby_id(corsika_id)
        
        return self.ds_kg(corsika_id, E)
    
    def ds_kg(self, corsika_id, E):
        param = self.params[corsika_id]
        # rigidity of the knee
        r_knee = 5.8e6
        # knee in energy
        E0 = param[2]*r_knee
        # add extragalactic light (proton only) component
        add = 0
        ddelta = 1./3
        if (corsika_id==14):
            #add = (5.0e-19-0.5e-19)*((E/1e8)**(-2.7))
            add = (5.0e-19)*((E/1e8)**(-2.7))
        if E < E0:
            return param[0]*((E/1e4)**param[1])+add
        else:
            # the normalization has to be changed to grant continuity of the spectrum
            norm = param[0]*((E0/1e4)**(2-ddelta))
            gamm = param[1]+ddelta-2
            return norm*((E/1e4)**gamm)+add
