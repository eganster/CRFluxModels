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
        # params = (a, \gamma, b), b only for the change of slope model (break)
        
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
        self.rigidity = {"cutoff" : {"ARGO"     : 1.30e15,
                                     "KASCADE:" : 15.1e15},
                         "break"  : {"ARGO"     : 640e12,
                                     "KASCADE:" : 5.8e15},}
        
    # exponential-square cutoff
    def cutoff(self, E, corsika_id):
        flux = self.params[corsika_id][0] * (E/10e3)**(-self.params[corsika_id][1]) * np.exp(-1.0*(E/(corsika_id%100 * self.rigidity[model][fit]))**2)
        return flux
            
    def slope_break(self, E, corsika_id):
        flux = 
        return flux
    
    def nucleus_flux(self, corsika_id, E):
        corsika_id = self._find_nearby_id(corsika_id)

        flux = 0.0
        for i in range(1, 4):
            p = self.params[corsika_id][i]
            flux += p[0] * E ** (-p[1] - 1.0) * \
                np.exp(-E / p[2] / self.rid_cutoff[i])
        return flux