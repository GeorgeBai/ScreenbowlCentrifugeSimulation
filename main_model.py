from math import *
import numpy as np
import class_function as cf


# parameters
feed_rate = 8e-3/60         # L/min to m3/s
feed_solid_fraction = 0.2

solid_density = 1410
liquid_density = 1000
density_diff = solid_density-liquid_density

liquid_viscosity = 1e-3       # Pa.s
max_sed_frac = 0.65              # important

# feed size distribution function ***** important
min_size = 1e-7
max_size = 700e-6
size_class = 50
size_interval = (max_size-min_size)/size_class

init_size = np.zeros([size_class])
for size in range(size_class):
    size_fraction_1 = 100/max_size*(min_size+size*size_interval)         # Accumulative function
    size_fraction_2 = 100/max_size*(min_size+(size+1)*size_interval)     # Accumulative function
    size_fraction = abs(size_fraction_2-size_fraction_1)
    init_size[size] = size_fraction
feed_size = np.divide(init_size, np.sum(init_size))
# print(np.sum(init_size))
# print(init_size)

# centrifuge parameters
rotation_speed = 2400         # rpm
rpm_diff = 10
rad = 2*pi*rotation_speed/60

bowl_radius = 0.075
weir_radius = 0.064
bowl_length = 0.2511
pitch = 0.067
blade_thickness = 0.0025
conical = 8/180*pi
conical_length = 0.1459
blade_to_wall = 0.002
screen_radius = 0.0545
screen_length = 0.203
aperture = 0.0025

bowl_number_winding = bowl_length / (pitch + blade_thickness)
bowl_channel_length = ((2 * pi * bowl_radius) ** 2 + pitch ** 2) ** 0.5 * bowl_number_winding
bowl_volume = bowl_channel_length * pitch * bowl_radius - weir_radius
beta = asin(pitch/bowl_channel_length/bowl_number_winding)
fraction_co = 0.1

# constants for Compartment calc
number_of_compartment = 500
compartment_length = bowl_channel_length/number_of_compartment
slice_thickness = 0.002
time_interval = 0.1

compartment_holder = {}
for comp in range(number_of_compartment):
    compartment_holder["comp_{}".format(comp+1)] = cf.Compartment(id=comp+1,
                                                                  sediment_rad=bowl_radius,
                                                                  sed_size=np.zeros([size_class]),
                                                                  fluid_size=np.zeros([size_class]),
                                                                  fluid_size_out=np.zeros([size_class]),
                                                                  sed_size_out=np.zeros([size_class]))

# print(compartment_holder['comp_2'].id)
# print(compartment_holder['comp_{}'.format(2)].sediment_radius)


# reconstructing size distribution
def re_size_distribution(volume_1, sv_1, size_distr_1,
                         volume_2, sv_2, size_distr_2):
    # inlet/out=[V_fluid, sv_fluid, PSD_fluid, Q_cake, sv_cake, PSD_cake]
    solid_1 = volume_1*sv_1*size_distr_1
    solid_2 = volume_2*sv_2*size_distr_2
    solid_volume = solid_1 + solid_2
    return solid_volume / sum(solid_volume)


# calculating each compartment
def update_compartment(cls, inlet, state):
    outlet = [0., 0., 0., 0., 0., 0.]
    sed_solid_weight_by_size = np.zeros([size_class])
    remain_solid_weight_by_size = np.zeros([size_class])

    fluid_volume = inlet[0] + cls.remain_fluid
    if fluid_volume == 0:
        return None
    mixed_solid_fraction = (inlet[0] * inlet[1] + cls.remain_fluid * cls.fluid_solid_frac) / fluid_volume
    mixed_size = re_size_distribution(inlet[0], inlet[1], inlet[2],
                                      cls.remain_fluid, cls.fluid_solid_frac, cls.fluid_size)
    # state check
    if state == 0:
        fluid_height = fluid_volume / (pitch * compartment_length)
        # print('fluid height {}'.format(fluid_height))
        # v = (2*g*h)^0.5
        mean_velocity = (2 * rad ** 2 * (cls.sediment_radius - fluid_height/2)*fluid_height/2) ** 0.5
        # print(state, mean_velocity)
    else:
        mean_velocity = inlet[0]/time_interval/(cls.sediment_radius - weir_radius)/pitch
        # print(state, mean_velocity)
    number_of_slice = int(ceil((cls.sediment_radius - weir_radius) / slice_thickness))
    loc_slice_thickness = (cls.sediment_radius - weir_radius) / number_of_slice
    for j in range(number_of_slice):
        local_velocity = 1.5 * mean_velocity * \
                         (1 - ((weir_radius + (j + 1 / 2) * loc_slice_thickness) / cls.sediment_radius) ** 2)
        # print(local_velocity)
        local_residence = compartment_length / 2 / local_velocity
        solid_per_slice = fluid_volume/number_of_slice*mixed_solid_fraction
        solid_weight_in_class = solid_per_slice*mixed_size
        # print(solid_weight_in_class)
        if local_residence < time_interval:
            # flow out the current compartment
            outlet[0] = outlet[0] + fluid_volume / number_of_slice
        else:
            # calc separation size in this compartment
            hendered_factor = (1 - mixed_size / max_sed_frac) ** 4.56  # gel_point better change to max sv
            separate_size = (log(cls.sediment_radius / (weir_radius + (j + 1 / 2) * loc_slice_thickness)) *
                             18 * liquid_viscosity /
                             (density_diff * hendered_factor * rad ** 2 * local_residence)) ** 0.5
            separate_size = max(separate_size)
            print(separate_size)
            for k in range(size_class):
                if min_size + (k+1/2)*size_interval >= separate_size:
                    sed_solid_weight_by_size[k] += solid_weight_in_class[k]
                    
                else:
                    remain_solid_weight_by_size[k] += solid_weight_in_class[k]
    # print(cls.id)
    sed_solid = np.sum(sed_solid_weight_by_size)
    # print(sed_solid)
    sed_solid_size = sed_solid_weight_by_size/sed_solid
    # sed_liquid = sed_solid/max_sed_frac*(1-max_sed_frac)
    remained_fluid = fluid_volume - outlet[0] - sed_solid/max_sed_frac
    remain_solid_frac = np.sum(remain_solid_weight_by_size)/remained_fluid
    remain_fluid_size = remain_solid_weight_by_size / np.sum(remain_solid_weight_by_size)

    outlet[1] = mixed_solid_fraction
    outlet[2] = mixed_size

    # calculation for sedimentation part
    remain_sed_pre = (bowl_radius-cls.sediment_radius)*compartment_length*pitch
    sed_volume = inlet[3] + remain_sed_pre
    if sed_volume > 0:
        sed_solid_frac = (inlet[3]*inlet[4]+remain_sed_pre*cls.sed_solid_frac)/sed_volume
        mixed_sed_size = re_size_distribution(inlet[3], inlet[4], inlet[5],
                                              remain_sed_pre, cls.sed_solid_frac, cls.sed_size)
        area_trans = sed_volume/compartment_length
        effective_velocity = pitch*rpm_diff*2*pi/60/((1+tan(beta)*tan(atan(fraction_co)+beta))*sin(beta))
        # cake transfer residence time
        sed_residence = compartment_length/effective_velocity
        if sed_residence > time_interval:
            sed_residence = time_interval
        sed_trans_volume = area_trans*effective_velocity*sed_residence
    else:
        sed_trans_volume = 0
        sed_solid_frac = 0
        mixed_sed_size = np.zeros([size_class])

    remained_sed_volume = sed_volume + sed_solid/max_sed_frac - sed_trans_volume
    remain_sed_radius = remained_sed_volume/(pitch*compartment_length)
    remain_sed_solid_frac = ((sed_volume-sed_trans_volume)*sed_solid_frac+sed_solid)/remained_sed_volume
    remain_sed_size = re_size_distribution(sed_volume-sed_trans_volume, sed_solid_frac, mixed_sed_size,
                                           sed_solid/max_sed_frac, max_sed_frac, sed_solid_size)

    # update compartment
    # remaining in properties
    cls.sediment_radius = remain_sed_radius
    cls.remain_fluid = remained_fluid
    cls.sed_solid_frac = remain_sed_solid_frac
    cls.fluid_solid_frac = remain_solid_frac
    cls.sed_size = remain_sed_size
    cls.fluid_size = remain_fluid_size
    # flow properties
    cls.fluid_out = outlet[0]
    cls.fluid_out_size = outlet[1]
    cls.fluid_out_frac = outlet[2]
    cls.sed_out = sed_trans_volume
    cls.sed_out_size = mixed_sed_size
    cls.sed_out_frac = sed_solid_frac


# calculating filling process
def bowl_section_calculation(state):
    # cls is curent class
    # inlet/out list contains inlet/out fluid rate, in/out solids fraction, in/out size distribution,
    # inlet/out cake rate, in/out solids fraction, in/out size distribution
    # inlet/out=[V_fluid, sv_fluid, PSD_fluid, Q_cake, sv_cake, PSD_cake]
    for cp_ in range(1, number_of_compartment+1):
        comp_calc = compartment_holder['comp_{}'.format(cp_)]
        if cp_ == 1:
            inlet_sed = compartment_holder['comp_2'].sed_out
            inlet_sed_frac = compartment_holder['comp_2'].sed_out_frac
            inlet_sed_size = compartment_holder['comp_2'].sed_out_size

            inlet = [feed_rate*time_interval, feed_solid_fraction, feed_size,
                     inlet_sed, inlet_sed_frac, inlet_sed_size]
            # print(inlet[0])
            update_compartment(comp_calc, inlet, state)

        elif cp_ == number_of_compartment:
            comp_calc_pre = compartment_holder['comp_{}'.format(cp_ - 1)]
            # comp_calc_next = compartment_holder['comp_{}'.format(cp_ + 1)]
            inlet_fluid = comp_calc_pre.fluid_out
            inlet_fluid_frac = comp_calc_pre.fluid_out_frac
            inlet_fluid_size = comp_calc_pre.fluid_out_size
            inlet_sed = 0
            inlet_sed_frac = 0
            inlet_sed_size = np.zeros([size_class])

            inlet = [inlet_fluid, inlet_fluid_frac, inlet_fluid_size,
                     inlet_sed, inlet_sed_frac, inlet_sed_size]

            update_compartment(comp_calc, inlet, state)
        else:
            comp_calc_pre = compartment_holder['comp_{}'.format(cp_-1)]
            comp_calc_next = compartment_holder['comp_{}'.format(cp_+1)]
            inlet_fluid = comp_calc_pre.fluid_out
            inlet_fluid_frac = comp_calc_pre.fluid_out_frac
            inlet_fluid_size = comp_calc_pre.fluid_out_size
            inlet_sed = comp_calc_next.sed_out
            inlet_sed_frac = comp_calc_next.sed_out_frac
            inlet_sed_size = comp_calc_next.sed_out_size
            
            inlet = [inlet_fluid, inlet_fluid_frac, inlet_fluid_size,
                     inlet_sed, inlet_sed_frac, inlet_sed_size]
            if inlet[0] == 0:
                # print('pass no inlet')
                pass
            else:
                update_compartment(comp_calc, inlet, state)
            
            
# feed in initial state
# feed_stream = [feed_rate*time_interval, feed_solid_fraction, feed_size]

# set simulating time
time_step = 1000
total_time = time_interval*time_step/60    # mins
for t_step in range(time_step):
    print('time step {}'.format(t_step+1))
    # check the filling level
    if compartment_holder['comp_{}'.format(number_of_compartment)].sediment_radius - \
            compartment_holder['comp_{}'.format(number_of_compartment)].fluid_height(pitch, compartment_length) \
            > weir_radius:
        # filling up process calculation
        state_indicator = 1
    else:
        # normal calculation
        state_indicator = 1
    bowl_section_calculation(state_indicator)
    print('liquid level in last compartment: {}'.
          format(compartment_holder['comp_{}'.format(number_of_compartment)].fluid_height(pitch, compartment_length)))
    print('liquid level in 4th compartment: {}'.
          format(compartment_holder['comp_{}'.format(2)].fluid_height(pitch, compartment_length)))