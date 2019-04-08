from CRFluxModels import *

class MascarettiBlasiEvoli(PrimaryFlux):
    """
    To be updated
    
    # p/H - 14
    # He - 402
    # light component (CNO) - 1206 (C)
        
    """

    def __init__(self, model="cutoff", fit="ARGO"):

        self.name = "MascarettiBlasiEvoli"
        self.model = model
        self.fit = fit
        self.params = {}

        # check for model and fit string to be valid
        if self.model not in ["cutoff", "break"] or self.fit not in ["ARGO", "KASCADE"]:
            raise Exception(
                    'MascarettiBlasiEvoli(): Unknown model/fit version requested.\
                      Allowed: model = [\'cutoff\', \'break\'] and fit = [\'ARGO\', \'KASCADE\']')
        
        # set params
        # params = (a, \gamma, b, E_lim)
        
        # p and He always the same
        self.params[14]   = (1.50, 2.71) # H
        self.params[402]  = (1.50, 2.64) # He
        
        # "light" component different for both models
        if self.model == "cutoff":
            self.params[1206] = (6.00, 2.70) # CNO
        elif self.model == "break":
            self.params[1206] = (5.00, 2.70) # CNO
            
        self.nucleus_ids = self.params.keys()
        
        # rigiditys different
        self.rigidity = {"cutoff" : {"ARGO"    : 1.30e15,
                                     "KASCADE" : 15.1e15},
                         "break"  : {"ARGO"    : 640e12,
                                     "KASCADE" : 5.8e15}}
        
        # calculate b_i, E_lim and set it
        for corsika_id in self.nucleus_ids:
            # get Z, A
            Z, A = self.Z_A(corsika_id)
            E_lim = Z * self.rigidity[self.model][self.fit] / 1e9
            # b = (E/10TeV)**(2 - delta)
            b = self.params[corsika_id][0] * (E_lim/10e3)**(2.0 - 1.0/3.0)
            # set it
            self.params[corsika_id] += (b, E_lim)
    
        
    # exponential-square cutoff
    def cutoff(self, corsika_id, E):
        # get Z, A
        Z, A = self.Z_A(corsika_id)
        # a * (E/10TeV)**-gamma * exp(-(E/ZR)**2) * 10**-7
        flux = self.params[corsika_id][0] * (E/10e3)**(-self.params[corsika_id][1]) * np.exp(-1.0*(E/(Z * self.rigidity[self.model][self.fit] / 1e9))**2) * 1e-7
        return flux
    
    # change of slope          
    def slope_break(self, corsika_id, E):
        # check for case
        
        flux = np.where(E <= self.params[corsika_id][3], self.params[corsika_id][0] * (E/10e3)**(-self.params[corsika_id][1]) * 1e-7,
                                                         self.params[corsika_id][2] * (E/10e3)**(-self.params[corsika_id][1] + 1.0/3.0 - 2.0) * 1e-7)
        
        #if E <= self.params[corsika_id][3]:
        #    # a * (E/10TeV)**(-gamma + delta -2) * 20**-7
        #    flux = self.params[corsika_id][0] * (E/10e3)**(-self.params[corsika_id][1] + 1.0/3.0 - 2.0) * 1e-7
        #else:
        #    # b * (E/10TeV)**(-gamma + delta -2) * 20**-7
        #    flux = self.params[corsika_id][2] * (E/10e3)**(-self.params[corsika_id][1] + 1.0/3.0 - 2.0) * 1e-7
        
        return flux
    
    # extragalactic component
    def eg_component(self, corsika_id, E):
        # a * (E/100PeV)**-gamma
        flux = self.params[corsika_id][0] * (E/100e6)**(-self.params[corsika_id][1]) * 1e-19
        return flux
        
    # calculate nucleus flux
    def nucleus_flux(self, corsika_id, E):
        # get nearest corsika_id
        corsika_id = self._find_nearby_id(corsika_id)
        # set flux to 0.0
        flux = 0.0
        # check for light eg component
        if corsika_id == 1206:
            flux +=  self.eg_component(corsika_id, E)
        else:
            # calculate flux
            if self.model == "cutoff":
                flux += self.cutoff(corsika_id, E)
            else:
                flux += self.slope_break(corsika_id, E)
        return flux