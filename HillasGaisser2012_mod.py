from CRFluxModels import *

class HillasGaisser2012_mod(HillasGaisser2012):
    """changes 2nd/3rd component to p/He only
    """
    
    def __init__(self, model="H4a", cids_to_remove=None, comp_to_remove=None, dist=None):
        
        """
        Args:
            model (str): H3a or H4a
            cids_to_remove (list): list of corsika_ids that flux should be removed (14, 402, 1206, 2814, 5426)
            comp_to_remove (list): list of components that shoul be removed (1, 2, 3)
            dist (dict): dict with redistribution of removed flux to corsika ids, e.g. {14: 0.5, 402:0.5}
        """
        
        # convert None to default values
        if not cids_to_remove:
            cids_to_remove = []
        if not comp_to_remove:
            comp_to_remove = []
        else:
            comp_to_remove = [v-1 for v in comp_to_remove]
        if not dist:
            dist = {}
        
        # check flux re-distribution
        dist_tot = sum(dist.values())
        if dist_tot != 1.0 and cids_to_remove:
            raise Exception("Make sure that exactly 100% of the flux is redistributed, not {}%".format(dist_tot*100))
            
        # call __init__ from HillasGaisser2012 and setup H3a/H4a model
        super(new_HillasGaisser2012, self).__init__(model)
        
        # set additional parameters
        self.cids_to_remove = cids_to_remove
        self.comp_to_remove = comp_to_remove
        self.dist = dist

        
    def nucleus_flux(self, corsika_id, E):
        # call self.nucleus_flux_components and sum
        return np.sum(self.nucleus_flux_components(corsika_id, E), axis=0)

    
    def nucleus_flux_components(self, corsika_id, E):
        corsika_id = self._find_nearby_id(corsika_id)
        #Z, A = self.Z_A(corsika_id)
        
        # get original flux
        old_flux_comp = np.array(super(new_HillasGaisser2012, self).nucleus_flux_components(corsika_id, E))
        
        # check for corsika_id to be removed
        if corsika_id in self.cids_to_remove:
            for comp in self.comp_to_remove:
                old_flux_comp[comp] = 0.0
            return old_flux_comp
        
        # elif corsika_id's flux gets not modified
        elif corsika_id not in self.dist.keys():
            return old_flux_comp
        
        # else
        else:
            # calculate removed flux
            removed_flux_comp = np.zeros_like(old_flux_comp)
            for cid in self.cids_to_remove:
                removed_flux_comp += np.array(super(new_HillasGaisser2012, self).nucleus_flux_components(cid, E))
            
            new_flux_comp = np.copy(old_flux_comp)
            
            # calculate new flux
            for comp in self.comp_to_remove:
                new_flux_comp[comp] += self.dist[corsika_id]*removed_flux_comp[comp]
                
            # return new flux
            return new_flux_comp