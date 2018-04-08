# class for hold compartment
class Compartment:
    def __init__(self, id=0, sediment_rad=0, remain_fluid=0,
                 sed_solid_frac=0, fluid_solid_frac=0,
                 sed_size=None, fluid_size=None,
                 fluid_out=0, fluid_size_out=None,
                 fluid_out_frac=0,
                 sed_out=0, sed_size_out=None,
                 sed_out_frac=0):
        self.id = id
        # remaining in properties
        self.sediment_radius = sediment_rad
        self.remain_fluid = remain_fluid
        self.sed_solid_frac = sed_solid_frac
        self.fluid_solid_frac = fluid_solid_frac
        self.sed_size = sed_size
        self.fluid_size = fluid_size
        # flow properties
        self.fluid_out = fluid_out
        self.fluid_out_size = fluid_size_out
        self.fluid_out_frac = fluid_out_frac
        self.sed_out = sed_out
        self.sed_out_size = sed_size_out
        self.sed_out_frac = sed_out_frac
        
    def fluid_height(self, width, length):
        return self.remain_fluid/(width*length)
            


