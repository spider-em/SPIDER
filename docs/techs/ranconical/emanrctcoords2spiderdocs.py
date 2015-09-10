#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

print "emanrctcoords2spiderdocs.py, Modified 2015 July 14"  
# Syntax:
#   JSON files: http://stackoverflow.com/questions/2835559/parsing-values-from-a-json-file-in-python
#   timestamps: http://stackoverflow.com/questions/12400256/python-converting-epoch-time-into-the-datetime

import sys
import json
import shutil
import os
from time import gmtime, strftime, localtime

def backup(filename):
    if os.path.exists(filename):
        found_vacancy = 0
        tiebreaker = 0
        shortdir = os.path.basename(os.path.dirname(filename))
        while not found_vacancy:
            test_filename = filename + '_' + str(tiebreaker)
            if os.path.exists(test_filename):
                tiebreaker += 1
            else:
                found_vacancy = 1
                short_old = os.path.join(shortdir, os.path.basename(filename))
                short_new = os.path.join(shortdir, os.path.basename(test_filename))
                print short_old, 'renamed to', short_new
                #os.rename(filename, test_filename)
                #shutil.copyfile(filename, test_filename)
                shutil.copy2(filename, test_filename)

def writedoc(list, filename):
    extension = os.path.splitext(filename)[1][1:]
    backup(filename)
    fp = open(filename, mode='w')
    time = strftime("%Y-%m-%d AT %H:%M:%S", localtime())
    list.insert(0," ;rct/%s   %s   %s\n" % (extension, time, filename))
    fp.writelines(list)

###################################################################################################

args  = sys.argv[1:]
nargs = len(args)
#print "%s nargs in %s: %s " % (nargs,type(args),args)

if nargs < 4:
    print 'Usage: '
#   args:                                   0              1              2             3                4                 5                6            7
    print 'emanrctcoords2spiderdocs.py json_file_1 ydim_micrograph_1 json_file_2 ydim_micrograph_2 output_directory output_filenumber {flip,noflip} {optional_extension}'
    sys.exit()

# build list
coords_files = [args[0],args[2]]

output_directory = args[4]
filenum          = int(args[5])

file_counter = 0
pair_coordinates = []
tilt_params = []

for file in coords_files:
    file_counter += 1
    with open(file) as data_file: data = json.load(data_file)
    
    json_coordinates = data["boxes_rct"]  # type 'list'
    
    particle_counter = 0
    coordinate_array = []
    
    if file_counter == 1 : 
        micrograph_ydim = int(args[1])
    else:
        micrograph_ydim = int(args[3])
    
    # loop through coordinates
    for set in json_coordinates:
        particle_counter += 1

        # read coordinates
        xcoord = set[0]
        
        if args[6] == 'noflip' :
            ycoord = set[1]
        elif args[6] == 'flip' :
            ycoord = micrograph_ydim - set[1]  # WAS set[1]
        else: 
            print "flip/noflip parameter not set"
            sys.exit()
            

        which_pair = set[2]
        
        # format SPIDER doc
        data_string = " %4i 3 %12i %15.2f %12.2f " % (particle_counter, particle_counter, xcoord, ycoord)
        coordinate_array.append(data_string + "\n")
    
#    coordinate_array.append(" 0001 6            1          501          103          501          103            1\n")
        
    pair_coordinates.append(coordinate_array)
#    print "which_pair:", which_pair
    #if which_pair == 'untilted' :
        ##
    #else:  # tilted micrograph
        ##
    
    # remove coordinates from dictionary
    del data["boxes_rct"]
    if len(data) == 1:
        
        # read tilt parameters
        for key in data: 
            
            # check whether remaing key is tilt parameters
            if key[:11] != 'tiltparams_':
                print "ERROR: json file format has changed"
                sys.exit()
            else:  # continue
                mic_name = key[11:]
                angle_set = data[key]

                if file_counter == 1 : 
                    first_file_angles = []
                    first_boxfile = mic_name
                    
                    for angle in angle_set:
                        first_file_angles.append(angle)
                        
                    tilt_params.append(" 0124 3     %7.3f      %7.3f       %7.3f\n" % (first_file_angles[0], first_file_angles[1], first_file_angles[2]))
                else :  # if second file read
                    second_boxfile = mic_name
                    notOK = False  # initialize
                    
                    for index in range(len(angle_set)):
                        current_angle = angle_set[index]
                        if current_angle != first_file_angles[index] : notOK = True
                        
                    if notOK : 
                        print "ERROR: angles in tilt pair %s and %s not identical" % (first_boxfile, second_boxfile)
                
                stem = os.path.splitext(mic_name)[0]
                if nargs >= 6: 
                    extension = args[7]
                else:  # extract from filename
                    extension = os.path.splitext(mic_name)[1][1:]
                
                if which_pair == 'tilted' :
                    prefix = 'dct'
                else:  # untilted
                    prefix = 'dcu'
                
                name_wo_dir = prefix + "%03d" % filenum + '.' + extension
                spider_coords = os.path.join(output_directory,name_wo_dir)
                print "spider_coords:", spider_coords
                coordinate_array.insert(0," ;              PARTNUM      X_COORD      Y_COORD\n")
                writedoc(coordinate_array, spider_coords)
    else: 
        print "ERROR: json file format has changed"
    
name_wo_dir = 'dcb' + "%03d" % filenum + '.' + extension
tilt_params_doc = os.path.join(output_directory,name_wo_dir)
print "tilt_params_doc:", tilt_params_doc
tilt_params.insert(0, " ;            THETA       GAMMA          PHI\n")
writedoc(tilt_params, tilt_params_doc)


