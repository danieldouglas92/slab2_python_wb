############################################# MODEL FEATURES FOR THE WB FILE #############################################
################### Functions which write strings for various model features to the world builder file ###################
##########################################################################################################################

# Writes the string to a world builder file that initializes a feature

# model_name = string specifying the name of the feature. 
# feature_name = user specified string naming the feature.
# min_depth = minimum depth extent of the feature
# max_depth = maximum depth extenet of the feautre
# coordinates = array specifying the lateral coordinates of the feature
# is_subducting = boolean, whether the feature is a subducting plate
# dip_point = if is_subducting=True, location of the dip point of the subducting slab
def model_feature_string(model_name, feature_name, min_depth, max_depth, coordinates, is_subducting, dip_point):
    if is_subducting == False:
        string = '"model":"' + str(model_name) + '", "name":"' + str(feature_name) + '", "min depth":' + str(min_depth) +', "max depth":' + \
                 str(max_depth) + ', "coordinates":' + str(coordinates) + ',\n'
    else:
        string = '"model":"' + str(model_name) + '", "name":"' + str(feature_name) + '", "min depth":' + str(min_depth) + ', "max depth":' + \
                 str(max_depth) + ', "coordinates":' + str(coordinates) + ', "dip point":' + str(dip_point) + ',\n'
    return string

# Extracts the location of the trench from the files containing slab profiles

# profile_directory = the pathway to the directory containing the profiles across slab2.0 data
# xshift = amount of shift in the x direction
# yshift = amount of shift in the y direction
# trench_direction = direction of the trench relative to the overriding plate
# coordinate_sys = string specifying whether the worldbuilder file is 'Cartesian' or 'Spherical'
def trench_extractor(profile_directory, xshift, yshift, trench_direction, coordinate_sys):
    import os
    import numpy as np
    trench_x = []
    trench_y = []
    if coordinate_sys == 'Cartesian':
        for file in np.sort(os.listdir(profile_directory)):
            profile_file = np.loadtxt(fname=profile_directory + file)
            if trench_direction == 'East':
                trench_x.append( (np.max(profile_file[:, 0]) + xshift) * 1000)
                trench_y.append( (profile_file[:, 1][np.where(profile_file[:, 0] == np.max(profile_file[:, 0]))][0] + yshift) * 1000)
                
            elif trench_direction == 'West':
                trench_x.append( (np.min(profile_file[:, 0]) + xshift) * 1000)
                trench_y.append( (profile_file[:, 1][np.where(profile_file[:, 0] == np.min(profile_file[:, 0]))][0] + yshift) * 1000)

            elif trench_direction == 'South':
                trench_y.append( (np.min(profile_file[:, 1]) + yshift) * 1000)
                trench_x.append( (profile_file[:, 0][np.where(profile_file[:, 1] == np.min(profile_file[:, 1]))][0] + xshift) * 1000)

            elif trench_direction == 'North':
                trench_y.append( np.max(profile_file[:, 1] + yshift) * 1000)
                trench_x.append( (profile_file[:, 0][np.where(profile_file[:, 1] == np.min(profile_file[:, 1]))][0] + xshift) * 1000)
    else:
        for file in np.sort(os.listdir(profile_directory)):
            profile_file = np.loadtxt(fname=profile_directory + file)
            if trench_direction == 'East':
                trench_x.append( (np.max(profile_file[:, 0]) + xshift))
                trench_y.append( (profile_file[:, 1][np.where(profile_file[:, 0] == np.max(profile_file[:, 0]))][0] + yshift))
                
            elif trench_direction == 'West':
                trench_x.append( (np.min(profile_file[:, 0]) + xshift))
                trench_y.append( (profile_file[:, 1][np.where(profile_file[:, 0] == np.min(profile_file[:, 0]))][0] + yshift))

            elif trench_direction == 'South':
                trench_y.append( (np.min(profile_file[:, 1]) + yshift))
                trench_x.append( (profile_file[:, 0][np.where(profile_file[:, 1] == np.min(profile_file[:, 1]))][0] + xshift))

            elif trench_direction == 'North':
                trench_y.append( np.max(profile_file[:, 1] + yshift))
                trench_x.append( (profile_file[:, 0][np.where(profile_file[:, 1] == np.min(profile_file[:, 1]))][0] + xshift))
        
    return trench_x, trench_y

# Splits the trench values to allow for a better specification of the overriding and subducting plate boundaries

# trench_x = array containing x-points of the trench output by trench_extractor function
# trench_y = array containing y-points of the trench output by trench_extractor function
# x_bounds = array with the max and min x-values of the world
# y_bounds = array with the max and min y-values of the world
# trench_direction = string specifying the direction of the trench relative to the overriding plate

def trench_splitter(trench_x, trench_y, x_bounds, y_bounds, trench_direction):
    import numpy as np
    
    if trench_direction == 'East':
        sub_coords = [[np.max(x_bounds), np.min(y_bounds)]]
        over_coords = [[np.min(x_bounds), np.min(y_bounds)]]
        sub_coords.append([trench_x[0], np.min(y_bounds)])
        over_coords.append([trench_x[0], np.min(y_bounds)])
        for i in range(len(trench_x)):
            sub_coords.append([trench_x[i], trench_y[i]])
            over_coords.append([trench_x[i], trench_y[i]])
        sub_coords.append([trench_x[-1], np.max(y_bounds)])
        over_coords.append([trench_x[-1], np.max(y_bounds)])
        sub_coords.append([np.max(x_bounds), np.max(y_bounds)])
        over_coords.append([np.min(y_bounds), np.max(y_bounds)])
        
    elif trench_direction == 'West':
        sub_coords = [[np.min(x_bounds), np.min(y_bounds)]]
        over_coords = [[np.max(x_bounds), np.min(y_bounds)]]
        sub_coords.append([trench_x[0], np.min(y_bounds)])
        over_coords.append([trench_x[0], np.min(y_bounds)])
        for i in range(len(trench_x)):
            sub_coords.append([trench_x[i], trench_y[i]])
            over_coords.append([trench_x[i], trench_y[i]])
        sub_coords.append([trench_x[-1], np.max(y_bounds)])
        over_coords.append([trench_x[-1], np.max(y_bounds)])
        sub_coords.append([np.min(x_bounds), np.max(y_bounds)])
        over_coords.append([np.max(x_bounds), np.max(y_bounds)])
    
    elif trench_direction == 'North':
        sub_coords = [[np.min(x_bounds), np.min(y_bounds)]]
        over_coords = [[np.min(x_bounds), np.max(y_bounds)]]
        sub_coords.append([np.min(x_bounds), trench_y[0]])
        over_coords.append([np.min(x_bounds), trench_y[0]])
        for i in range(len(trench_y)):
            sub_coords.append([trench_x[i], trench_y[i]])
            over_coords.append([trench_x[i], trench_y[i]])
        sub_coords.append([np.max(x_bounds), trench_y[-1]])
        over_coords.append([np.max(x_bounds), trench_y[-1]])
        sub_coords.append([np.max(x_bounds), np.min(y_bounds)])
        over_coords.append([np.max(x_bounds), np.max(y_bounds)])
        
    elif trench_direction == 'South':
        over_coords = [[np.min(x_bounds), np.min(y_bounds)]]
        sub_coords = [[np.min(x_bounds), np.max(y_bounds)]]
        over_coords.append([np.min(x_bounds), trench_y[0]])
        sub_coords.append([np.min(x_bounds), trench_y[0]])
        for i in range(len(trench_y)):
            over_coords.append([trench_x[i], trench_y[i]])
            sub_coords.append([trench_x[i], trench_y[i]])
        over_coords.append([np.max(x_bounds), trench_y[-1]])
        sub_coords.append([np.max(x_bounds), trench_y[-1]])
        over_coords.append([np.max(x_bounds), np.min(y_bounds)])
        sub_coords.append([np.max(x_bounds), np.max(y_bounds)])
    return sub_coords, over_coords

########################################## TEMPERATURE FEATURES FOR THE WB FILE ##########################################
################ Functions which write strings for various temperature features to the world builder file ################
##########################################################################################################################

# Defines the temperature as either a geotherm generated from the half-space cooling or plate cooling model

# model_name = 'plate model' or 'half space model', specifiying which cooling model you want
# max_depth = the maximum depth to which this temperature feature will be used
# min_depth = the minimum depth to which this temperature feature will be used
# bot_temp = the temperature the maximum depth
# top_temp = the temperature at the minimum depth
# spr_vel = the velocity of the plate
# ridge_coords = the coordinates specifying the axis of the spreading ridge
# first_last = whether this feature is the first or last temperature feature. 'first' is the first feature, 'last' is the last feature, 'both'
# is if this is the only feature

def cooling_model(model_name, max_depth, min_depth, bot_temp, top_temp, spr_vel, ridge_coords, first_last):
    if first_last == 'first' or first_last == 'both':
        string = '"temperature models":['
    else:
        string = ''
        
    if model_name == 'plate model' or model_name == 'half space model':
        string += '{"model":"' + str(model_name) + '", "max depth":' + str(max_depth) + ', "min depth":' + str(min_depth) + ', "top temperature":' + \
                  str(top_temp) + ', "bottom temperature":' + str(bot_temp) + ', "spreading velocity":' + str(spr_vel) + ', "ridge coordinates":' + \
                  str(ridge_coords) + '}'
    if first_last == 'last' or first_last == 'both':
        string += '], \n'
    else:
        string += ',\n'
    return string

# Defines the temperature as a linear conductive geotherm

# model_name = 'linear'
# max_depth = the maximum depth to which this temperature feature will be used
# min_depth = the minimum depth to which this temperature feature will be used
# bot_temp = the temperature the maximum depth
# top_temp = the temperature at the minimum depth
# first_last = whether this feature is the first or last temperature feature. 'first' is the first feature, 'last' is the last feature, 'both'
# is if this is the only feature

def linear_model(model_name, max_depth, min_depth, bot_temp, top_temp, first_last):
    if first_last == 'first' or first_last == 'both':
        string = '"temperature models":['
    string += '{"model":"' + str(model_name) + '", "max depth":' + str(max_depth) + ', "min depth":' + str(min_depth) + ', "top temperature":' + \
               str(top_temp) + ', "bottom temperature":' + str(bot_temp) + '}'
    if first_last == 'last' or first_last == 'both':
        string += '], \n'
    else:
        string += ',\n'
    return string

# Defines the temperature as a constant value

# model_name = 'uniform'
# uni_temp = the constant temperature value
# operation = whether uni_temp overrides the existing temperature field, or adds or subtracts from it. 'replace', 'add', or 'subtract'
# first_last = whether this feature is the first or last temperature feature. 'first' is the first feature, 'last' is the last feature, 'both'
# is if this is the only feature

def uniform_model(model_name, uni_temp, operation, first_last):
    if first_last == 'first' or first_last == 'both':
            string = '"temperature models":['
    else:
        string = ''
            
    if model_name == 'uniform':
        string += '{"model":"' + str(model_name) + '", "temperature":' + str(uni_temp) + ', "operation":"' + str(operation) + '"}'   
        
    if first_last == 'last' or first_last == 'both':
        string += '], \n'
        
    else:
        string += ',\n'
        
    return string
   
# Defines the temperature in a subducting slab

# density = density of the slab
# plate_vel = velocity of the plate entering the trench
# couple_depth = the depth where the subducting slab interfaces the overriding plate
# shallow_dip = the average dip of the slab up to the couple_depth
# ridge_coords = the coordinates specifying the axis of the spreading ridge
# max_slab_top = the distance to the top of the slab from the mid plane of the slab
# min_slab_top = the distance to the bottom of the slab from the mid plane of the slab
# first_last = whether this feature is the first or last temperature feature. 'first' is the first feature, 'last' is the last feature, 'both'
# is if this is the only feature    
    
def mass_conserving_model(density, plate_vel, couple_depth, shallow_dip, ridge_coords, taper, max_slab_top, min_slab_top, first_last):
    if first_last == 'first' or first_last == 'both':
        string = '"temperature models":['
    else:
        string = ''
        
    string += '{"model":"mass conserving", "density":' + str(density) + ', "plate velocity":' + str(plate_vel) + ', "coupling depth":' + str(couple_depth) + \
              ', "shallow dip":' + str(shallow_dip) + ', "ridge coordinates":' + str(ridge_coords) + ', "taper distance":' + str(taper) + \
              ', "max distance slab top":' + str(max_slab_top) + ', "min distance slab top":' + str(min_slab_top) + '}'
    
    if first_last == 'last' or first_last == 'both':
        string += '], \n'
    else:
        string += ', \n'
        
    return string

def subducting_plate_model(model_name, density, plate_vel, adiabatic_heating, first_last):
    if first_last == 'first' or first_last == 'both':
        string = '"temperature models":['
    else:
        string = ''
    string += '{"model":"' + str(model_name) + '", "density":' + str(density) + ', "plate velocity":' + str(plate_vel) + ', "adiabatic heating":' + str(adiabatic_heating) + '}'
    
    if first_last == 'last' or first_last == 'both':
        string += '], \n'
    else:
        string += ', \n'
    return string

########################################## COMPOSITION FEATURES FOR THE WB FILE ##########################################
############################# Functions which write strings segments for composition fields ##############################
##########################################################################################################################

# Writes the string for a composition feature to the world builder file
# model_name = string specifying the composition model, 'uniform'
# comp_index = the index of the current composition field. Indexing starts at 0
# max_depth = maximum depth extent of the compositional field
# min_depth = minimum depth extent of the compositional field
# is_subducting = boolean, whether the current field is a subducting slab
# first_last = whether this feature is the first or last temperature feature. 'first' is the first feature, 'last' is the last feature, 'both'
# is if this is the only feature  

def composition_feature_string(model_name, comp_index, max_depth, min_depth, is_subducting, first_last):
    if first_last == 'first' or first_last == 'both':
        string = '"composition models":['
    else:
        string = ''
    
    if is_subducting == False:
        string += '{"model":"' + str(model_name) + '", "compositions":[' + str(comp_index) + '], "max depth":' + str(max_depth) + \
                  ', "min depth":' + str(min_depth) + '}'
    else:
        string += '{"model":"' + str(model_name) + '", "compositions":[' + str(comp_index) + '], "max distance slab top":' + str(max_depth) + \
                  ', "min distance slab top":' + str(min_depth) + '}'
    if first_last == 'last' or first_last == 'both':
        string += ']\n'
    else:
        string += ',\n'
    return string

######################################## SUBDUCTION ZONE FEATURES FOR THE WB FILE ########################################
############################# Functions which write strings segments of the subducting slab ##############################
##########################################################################################################################

# Writes the string that defines segments of a slab to the world builder file

# When initializing a slab in Worldbuilder, you must specify the same number of segments for each section along strike
# of the slab, which is troublesome since some sections of slab are much longer than others. To get around this,
# the slab is initialized with the same number of segments as the section with the maximum number of segments. 
# Then, sections with less than this number of segments are assigned filler sections with lengths of 1 m until
# they reach the required number of segments
    
# length = the length of a given segment
# thickness = the thickness of a given segment
# angle = the dip of a given segment
# coordinate = the trench coordinate index coupled to the current section
# max_coord_num = the number of sections making up the slab
# segment_num = the number of segments making up each section
# holder = the current segment index of a section

def segment_string(length, thickness, angle, coordinate, max_coord_num, segment_num, holder):
    
    # Here we check if the given section has the maximum number of segments, that way we can avoid
    # adding the filler segments
    if len(thickness) == segment_num:
        string_total = ''
        
        # Here we check if we are at the last segment of the section, if not the string is standard
        if holder != len(thickness) - 1:
            string_total += '{"length":' + str(length[holder]) + ', "thickness":' + str(thickness[holder]) + \
                            ', "angle":' + str(angle[holder]) + '},\n'
            return string_total
        
        # Here we check if we are at the last segment of the section AND we are at the last section, if so
        # the string is modified to close out the feature
        elif holder == (len(thickness) - 1) and coordinate == max_coord_num:
            string_total += '{"length":' + str(length[holder]) + ', "thickness":' + str(thickness[holder]) + \
                        ', "angle":' + str(angle[holder]) + '}]}\n'
            return string_total
        
        # Here is the case where we are at the last segment but not the last section, so a , must be added to 
        # continue onto the next section
        else:
            string_total += '{"length":' + str(length[holder]) + ', "thickness":' + str(thickness[holder]) + \
                        ', "angle":' + str(angle[holder]) + '}]},\n'
            return string_total
        
    # This is the case where the given section has LESS than the maximum number of segments, so we must add filler
    # segments to ensure that the section has the correct amount of segments
    else:
        string_total = ''
        
        # Here we proceed as above, adding the segment strings until we reach the end of the array
        if holder <= len(thickness) - 1:
            if holder != len(thickness) - 1:
                string_total += '{"length":' + str(length[holder]) + ', "thickness":' + str(thickness[holder]) + \
                                ', "angle":' + str(angle[holder]) + '},\n'
                return string_total
            else:
                string_total += '{"length":' + str(length[holder]) + ', "thickness":' + str(thickness[holder]) + \
                            ', "angle":' + str(angle[holder]) + '},\n'
                return string_total

        # Now we have reached the end of the array, but need to fill out the amount of segments to reach the maximum
        # segment number. Add segments with lengths of 1m, thicknesses of 1m, and dips of 1 degree until this is the case.
        else:
            if holder != segment_num - 1:
                string_total += '{"length":' + str(1) + ', "thickness":' + str([1]) + \
                                ', "angle":' + str([1]) + '},\n'
                return string_total
            
            elif holder == (segment_num - 1) and coordinate == max_coord_num:
                string_total += '{"length":' + str(1) + ', "thickness":' + str([1]) + \
                            ', "angle":' + str([1]) + '}]}\n'
                return string_total
            
            else:
                string_total += '{"length":' + str(1) + ', "thickness":' + str([1]) + \
                            ', "angle":' + str([1]) + '}]},\n'
                return string_total


# Calculates dip, thickness, and length of segments then uses segment_string to write to the world builder file.
# Has the option to vary the slab thickness between certain depths. If the varialbe slab_thickness is a scalar,
# a constant thickness is assumed. To vary thickness, slab_thickness must be an array where each entry is an 
# array with length two. The entries of this array are the thickness of the slab and the depth to where 
# that thickness occurs. For example, to have the slab be 100km thick between 0 <= depth <= 200km, and then
# have the thickness increase to 150km for depth > 200km, the variable slab_thickness would need to be set as:
#
# slab_thickness = [ [100e3, 200e3], [150e3, 1e10] ]
#
# Where 1e10 was chosen to be a depth so high it would never be reached. Slab thickness and dip are varied gradually 
# down dip

# world_builder_file = the name of the world builder file
# profile_directory = directory containing the world builder file
# xshift = the amount of shift in the x direction, m
# yshift = the amount of shift in the y direction, m
# slab_thickness = the thickness of the slab, scaler for uniform thickness, array for variable, km
def segment_section(world_builder_file, profile_directory, xshift, yshift, slab_thickness):
    import numpy as np
    import os
    
    length = []
    for file in np.sort(os.listdir(profile_directory)):
        length.append(len(np.loadtxt(fname=profile_directory + file)))
    segment_num = max(length)
    max_coord_num = len(os.listdir(profile_directory)) - 1
    
    # This initializes the slab with the correct number of segments
    world_builder_file.write('"segments":[\n')
    for i in range(segment_num):
        if i != segment_num - 1:
            world_builder_file.write('{"length":1, "thickness":[1.0], "angle":[0]},\n')
        else:
            world_builder_file.write('{"length":1, "thickness":[1.0], "angle":[0]}],\n\n')

    world_builder_file.write('    "sections":[')
    index = 0
    for file in np.sort(os.listdir(profile_directory)):
        index += 1
        filename = os.path.join(profile_directory, file)

        track_x = np.loadtxt(fname=filename, usecols=0) + xshift
        track_y = np.loadtxt(fname=filename, usecols=1) + yshift
        track_z = np.loadtxt(fname=filename, usecols=2)

        slab_length = 0
        segment_length = []
        dip_holder = [0]
        
        # Here we check to see if the thickness is set to vary along the slab by checking if slab_thickness is
        # a scalar or not. Dip is also computed and stored for output to the world builder file here
        if hasattr(slab_thickness, "__len__"):
            thick_holder = [slab_thickness[index - 1][:, 0][0]]
            for i in range(1, len(track_z)):
                thick_index = np.min(np.where( slab_thickness[index - 1][:, 1] > np.abs(track_z[i]) ))
                thick_holder.append(slab_thickness[index - 1][:, 0][thick_index])
                slab_distance = np.sqrt( (track_x[i] - track_x[i - 1])**2 + (track_y[i] - track_y[i - 1])**2 )
                riserun = (track_z[i - 1] - track_z[i]) / (slab_distance)
                segment_length.append(slab_distance * 1000)
                dip_holder.append(np.arctan(abs(riserun)) * 180 / np.pi)
            
            segment_thickness = []
            for k in range(len(thick_holder) - 1):
                if k != int(len(thick_holder) - 1):
                    segment_thickness.append([thick_holder[k], thick_holder[k + 1]])
                else:
                    segment_thickness.append([thick_holder[k]])
         
        else:
            segment_thickness = []
            for i in range(1, len(track_z)):
                segment_thickness.append([slab_thickness])
                slab_distance = np.sqrt( (track_x[i] - track_x[i - 1])**2 + (track_y[i] - track_y[i - 1])**2 )
                riserun = (track_z[i - 1] - track_z[i]) / (slab_distance)
                segment_length.append(slab_distance * 1000)
                dip_holder.append(np.arctan(abs(riserun)) * 180 / np.pi)
                
        dips = []
        for k in range(len(dip_holder) - 1):
            if k != int(len(dip_holder) - 1):
                dips.append([dip_holder[k], dip_holder[k + 1]])
            else:
                dips.append([dip_holder[k]])
        
        # Write to the World Builder File
        world_builder_file.write('    {"coordinate":' + str(index - 1) + ',\n')
        world_builder_file.write('     "segments":[')
        for h in range(segment_num):
            world_builder_file.write('    ' + segment_string(segment_length, segment_thickness, dips, index - 1, max_coord_num, segment_num, h))