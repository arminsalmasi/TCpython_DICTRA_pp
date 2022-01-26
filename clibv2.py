import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq as find_root
from tc_python import *
import pandas as pd
from tkinter import Tk, filedialog, font
import tkinter as ttk
import os
from IPython.display import clear_output
from copy import deepcopy
import copy
from scipy.ndimage.filters import gaussian_filter1d
import glob
#from tkinter import simpledialog
import pdb
import random
import sys
from collections import defaultdict
#********************************************************************************************
def get_values_from_textfiles():
    '''Reads all values from a given path, 
    input:
        path(optional)
    output:dict{path,all_mfs,elnames,all_points,times,n_points,DICTRA_all_npms, 
                 DICTRA_phnames,DICTRA_all_mus, T, interstitials, substitutionals}
    '''
    data = {}
    data['path'] = get_path()
    data['all_mfs'] = np.loadtxt('MOLE_FRACTIONS.TXT', dtype = float)
    data['elnames'] = np.loadtxt('EL_NAMES.TXT', dtype = str)
    data['all_points'] = np.loadtxt('VOLUME_MIDPOINTS.TXT', dtype = float)
    data['times'] = np.loadtxt('TIME.TXT', dtype = float)
    data['n_points'] = np.loadtxt('VOLUMES_PER_REGION.TXT', dtype = float)
    data['DICTRA_all_npms'] = np.loadtxt('PHASE_FRACTIONS.TXT', dtype = float)
    phnames = []
    with open('PH_NAMES.TXT', 'r')as f:
        for line in f:
            phnames.append(line.strip())
    data['DICTRA_phnames'] = np.array(phnames)
    data['DICTRA_all_mus'] = np.loadtxt('CHEMICAL_POTENTIALS.TXT', dtype = float)
    #Loading the data from txt files'''
    data['T'] = np.loadtxt('T.DAT', dtype = int)+273
    int_idx = []
    if 'N' in data['elnames']:
        Nidx = np.where(data['elnames'] == 'N')[0]
        int_idx.append(Nidx[0])
    if 'C' in data['elnames']:
        Cidx = np.where(data['elnames'] == 'C')[0]
        int_idx.append(Cidx[0])
    if 'H' in data['elnames']:
        Hidx = np.where(data['elnames'] == 'H')[0]
        int_idx.append(Hidx[0])
    if 'O' in data['elnames']:
        Oidx = np.where(data['elnames'] == 'O')[0]
        int_idx.append(Oidx[0])
    if 'VA' in data['elnames']:
        VAidx = np.where(data['elnames'] == 'VA')[0]
        int_idx.append(VAidx[0])
    sub_idx = list(np.arange(len(data['elnames'])))
    data['interstitials'] = [data['elnames'][int_idx],int_idx]
    for idx in data['interstitials'][1]: 
        sub_idx.pop(idx)
    data['substitutionals'] = [data['elnames'][sub_idx],sub_idx]
    return data
#********************************************************************************************
def get_path():
    '''change path to the designated folder, root is users home
    output:
        path
    '''
    flag = True
    while flag:
        root = Tk() # pointing root to Tk() to use it as Tk() in program.
        root.withdraw() # Hides small tkinter window.
        root.attributes('-topmost', True) # Opened windows will be active. above all windows despite of selection.
        #help(root)
        open_file = filedialog.askdirectory(initialdir = '~') # Returns opened path as str
        if open_file:
            flag = False
    os.chdir(open_file)
    return os.getcwd()
#********************************************************************************************
def get_timeStamp(times):
    time1 = float(input('Time to plot/{}/timesteps/{}/:'.format(times[-1],len(times))))
    if time1>=0 and time1 <=times[-1]:
        tStamp, nearestTime = find_nearest(times,time1)
    return tStamp, nearestTime
#********************************************************************************************                      
def find_nearest(array, value):
    '''find index with vlaue closest to input value
    output:
        index, nearest value
    '''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx]                      
#********************************************************************************************                      
def calculate_u_fractions(*argv):
    mf = argv[0]
    sub_idx = argv[1]
    elnames =  argv[2]
    sub_sum = np.sum(mf[:,sub_idx],axis=1)
    #uf_dict = {}
    uf = []
    for nel,el in enumerate(elnames):
        #uf_dict[el] = mf[:,nel]/sub_sum[:]
        uf.append(mf[:,nel]/sub_sum[:]) 
    return np.array(uf)
#********************************************************************************************
#def get_values_of_timeStamp(*argv):
#    '''get value of a given timestamp
#    input: Stamp,all_mfs,elnames,all_points,n_points,all_npms,DICTRA_phnames,substitutionals,
#    output: dict{'tStamp_mfs','tStamp_DICTRA_npms','tStamp_DICTRA_mus','tStamp_points','tStamp_DICTRA_phnames'}
#    '''
#    tStamp, mfs, elnames, points, n_points, npms, phnames, mus, subs = argv   
#    dict = {}
#    '''gets molefractions, chemicalpotentials, phasefractions, volumemidpoints and time of the timeindex
#    removes ZZDICTRA from phasenames, calculate u-fractions'''
#    dict['tStamp_points'] = points[int(np.sum(n_points[:tStamp])):int(np.sum(n_points[:tStamp+1]))]*1e6    
#    idx1 = int(np.sum(n_points[:tStamp]))*int(len(phnames))
#    idx2 = int(np.sum(n_points[:tStamp+1]))*int(len(phnames))
#    dict['tStamp_DICTRA_npms'] = npms[idx1:idx2].reshape((-1, int(len(phnames))))
#    if 'ZZDICTRA_GHOST' in phnames[-1]:
#        dict['tStamp_DICTRA_phnames'] = deepcopy(phnames[:-1])
#    else:
#        dict['tStamp_DICTRA_phnames'] = deepcopy(phnames)
#    idx1 = int(np.sum(n_points[:tStamp]))*len(elnames)
#    idx2 = int(np.sum(n_points[:tStamp+1]))*len(elnames)
#    dict['tStamp_mfs'] = mfs[idx1:idx2].reshape((-1, len(elnames)))
#    dict['tStamp_mus'] = mus[idx1:idx2].reshape((-1, len(elnames)))
#    dict['tStamp_ufs'] = calculate_u_fractions(dict['tStamp_mfs'],subs[1],elnames) 
#    return dict
#********************************************************************************************
def get_tStamp_vals(dict_input):
    '''get value of a given timestamp
    input: dict{all_values} 
    output: dict{'tStamp_mfs','tStamp_DICTRA_npms','tStamp_DICTRA_mus','tStamp_points','tStamp_DICTRA_phnames'}
    '''
    dict = copy.deepcopy(dict_input)
    tStamp = dict['tStamp'][0]
    mfs =dict['all_mfs']
    elnames = dict['elnames']
    points = dict['all_points']
    n_points =dict['n_points']
    npms =dict['DICTRA_all_npms']
    phnames = dict['DICTRA_phnames']
    mus = dict['DICTRA_all_mus']
    subs =dict['substitutionals']   
    '''gets molefractions, chemicalpotentials, phasefractions, volumemidpoints and time of the timeindex
    removes ZZDICTRA from phasenames, calculate u-fractions'''
    dict['tStamp_points'] = points[int(np.sum(n_points[:tStamp])):int(np.sum(n_points[:tStamp+1]))]*1e6    
    idx1 = int(np.sum(n_points[:tStamp]))*int(len(phnames))
    idx2 = int(np.sum(n_points[:tStamp+1]))*int(len(phnames))
    dict['tStamp_DICTRA_npms'] = npms[idx1:idx2].reshape((-1, int(len(phnames))))
    if 'ZZDICTRA_GHOST' in phnames[-1]:
        dict['tStamp_DICTRA_phnames'] = deepcopy(phnames[:-1])
    else:
        dict['tStamp_DICTRA_phnames'] = deepcopy(phnames)
    idx1 = int(np.sum(n_points[:tStamp]))*len(elnames)
    idx2 = int(np.sum(n_points[:tStamp+1]))*len(elnames)
    dict['tStamp_mfs'] = mfs[idx1:idx2].reshape((-1, len(elnames)))
    dict['tStamp_mus'] = mus[idx1:idx2].reshape((-1, len(elnames)))
    dict['tStamp_ufs'] = calculate_u_fractions(dict['tStamp_mfs'],subs[1],elnames) 
    #dict.pop('all_mfs')
    #dict.pop('all_points')
    #dict.pop('times')
    #dict.pop('n_points')
    #dict.pop('DICTRA_all_npms') 
    #dict.pop('DICTRA_phnames')
    #dict.pop('DICTRA_all_mus') 
    for key in ['all_mfs','all_points','times','n_points','DICTRA_all_npms' ,'DICTRA_phnames','DICTRA_all_mus']:
        dict.pop(key) 
    return dict
#********************************************************************************************
def tccalc(*argv):
    '''input: Calculate single equilibrium point by point with input condition
            0-timeStamp_values(
                0-mfs
                1-npms
                2-mus
                3-points
                4-phnames_of_timeStamp)
            1-databas
            2-phases_to_suspend
            3-refrence_states_of_activity(only_Carbon, don't put SER here)
            4-T
            5-elnames
            6-pth
            7-poly3file
        output:
            0-stablePhaeses_in_points,
            1-npms_in_points,
            2-vpvs_in_points,
            3-phX_in_points,
            4-activities_with_ref,
            5-activities,
            6-tcCalculated_mus,
            7-tcCalculated_ws
    '''
    try:
        mfs = argv[0][0] #shape(npoints,nelemets)
        database = argv[1]
        number_of_points = len(argv[0][3])
        phases_to_suspend = argv[2]
        refrence_states_of_activity = argv[3]
        T = argv[4]
        elnames = argv[5]
        pth = argv[6]
        poly3_file=argv[7]
    except:
        print('error in calcualtion input data')
    else:
        print('calculating thermodynaics')
    with TCPython() as start:
        if poly3_file and os.path.isfile('{}/p.POLY3'.format(pth)):
            poly = start.select_database_and_elements(database, elnames).get_system().with_single_equilibrium_calculation()
            poly.run_poly_command('read {}/p.POLY3'.format(pth))
            #poly.run_poly_command('list-status ,cps,')
            poly.remove_all_conditions()
            poly.set_condition('N',1)
            poly.set_condition('P',1e5)
        else:
            print('calculation without poly3 file')
            poly = start.set_cache_folder( "./_cache2").select_database_and_elements(database, elnames).get_system().with_single_equilibrium_calculation()
        if len(phases_to_suspend)>0:
            for phase in phases_to_suspend:
                poly.set_phase_to_suspended(phase)
        poly.set_condition("T", T)
        stablePhaeses_in_points = {}
        npms_in_points = {}
        vpvs_in_points = {}
        phX_in_points = {}
        activities_with_ref = defaultdict(list)
        activities = defaultdict(list)
        tcCalculated_mus = defaultdict(list)
        tcCalculated_ws = defaultdict(list)
        for point in range(number_of_points):
            for nel, el in enumerate(elnames[:-1]):
                poly.set_condition("X({})".format(el), mfs[point, nel])
            singlePointEq = poly.calculate()
            stablePhases = singlePointEq.get_stable_phases()
            stablePhaeses_in_points[point] = stablePhases        
            if 'C' in elnames:
                for reference in refrence_states_of_activity:
                    activities_with_ref['ac(C,'+reference+')'].append(singlePointEq.get_value_of('ac(C,{})'.format(reference)))  
            for nel, el in enumerate(elnames[:]):
                activities[el].append(singlePointEq.get_value_of('ac({})'.format(el)))
                tcCalculated_mus[el].append(singlePointEq.get_value_of('mu({})'.format(el)))
                tcCalculated_ws[el].append(singlePointEq.get_value_of('w({})'.format(el)))
            for ph in stablePhases:      
                npms_in_points['{}, {}'.format(point, ph)] = singlePointEq.get_value_of('npm({})'.format(ph))
                vpvs_in_points['{}, {}'.format(point, ph)] = singlePointEq.get_value_of('vpv({})'.format(ph))
                temp1 = []
                temp2 = []
                for nel2, el2 in enumerate(elnames[:]):
                    temp1.append(singlePointEq.get_value_of('X({}, {})'.format(ph, el2)))
                phX_in_points['{}, {}'.format(point, ph)] = np.array(temp1)
    clear_output(wait = False)
    return (stablePhaeses_in_points, npms_in_points, vpvs_in_points, phX_in_points,
            activities_with_ref, activities, tcCalculated_mus, tcCalculated_ws)
#********************************************************************************************

        
#********************************************************************************************

def trim_tcCalculated_values(*argv):
    '''input:
            0-timeStamp_values(
                0-mfs
                1-npms
                2-mus
                3-points
                4-phnames_of_timeStamp )
            1-timStamp_tc_postProcessed_values(
                0-stablePhaeses_in_points
                1-npms_in_points
                2-vpvs_in_points
                3-phX_in_points
                4-activities_with_ref
                5-activities
                6-tcCalculated_mus 
                7- Ws)
            2-elnames 
            3-subidex 
        output:
            0-points
            1-mfs_dict
            2-npms_dict
            3-vpvs_dict
            4-activities_with_ref_dict
            5-activities_dict
            6-int_idx
            7-mus_dict
            8-tcCalculated_mus_dict
            9-ws_dict
            10-phX_in_points_dict
            11-elnames
            12-phase_names_tcCalculated'''
    try:
        mfs = argv[0][0]
        mus_dictra = argv[0][2]
        points = argv[0][3]
        stablePhaeses_in_points = argv[1][0]
        npms_in_points = argv[1][1]
        vpvs_in_points = argv[1][2]
        phX_in_points = argv[1][3]
        activities_with_ref = argv[1][4]
        activities = argv[1][5]
        tcCalculated_mus = argv[1][6]
        ws = argv[1][7]
        elnames = argv[2]
        int_idx = argv[3]
    except:
        print('input error')
    else:
        print('input OK') 
    
    mfs_dict= {}
    for nel, el in enumerate(elnames):
        mfs_dict[el] = mfs[:,nel]
    phase_names = []
    npoints = len(points)
    for i in range(npoints):
        for j in stablePhaeses_in_points[i]:
               if j not in phase_names:
                    phase_names.append(j)
    npms_dict = {}
    phX_in_points_dict = {}
    for ph in phase_names:
        npms_dict[ph] = np.zeros([npoints])
        phX_in_points_dict[ph] = np.zeros([npoints,len(elnames)])
    for nph,ph in enumerate(phase_names):
        for point in range(npoints):
            if '{}, {}'.format(str(point),ph) in npms_in_points:
                npms_dict[ph][point] = npms_in_points['{}, {}'.format(str(point),ph)]
                for nel,el in enumerate(elnames):
                    phX_in_points_dict[ph][point,nel] = phX_in_points['{}, {}'.format(str(point),ph)][nel]
    vpvs_dict = {}
    for ph in phase_names:
        vpvs_dict[ph] = np.zeros([npoints,len(phase_names)])
    for nph,ph in enumerate(phase_names):
        for point in range(npoints):
            if '{}, {}'.format(str(point),ph) in vpvs_in_points:
                vpvs_dict[ph][point,nph] = vpvs_in_points['{}, {}'.format(str(point),ph)]    
    
    return (points, mfs_dict, npms_dict, vpvs_dict, 
            activities_with_ref, activities,subidx,
            mus_dictra, tcCalculated_mus, ws, phX_in_points_dict,
            elnames,phase_names)
# #********************************************************************************************
def correct_phase_indices(*argv):
    '''
    input:
        0-timmed_tcCalculated_vals(
            0-points
            1-mfs_dict
            2-npms_dict
            3-vpvs_dict
            4-activities_with_ref_dict
            5-activities_dict
            6-subidx
            7-mus_dict
            8-tcCalculated_mus_dict
            9-ws_dict
            10-phX_in_points
            11-elnames
            12-phase_names_tcCalculated)
        1-phase_change([(phase_name,major_el,new_name),...])    
    '''
    try:
        points = argv[0][0]
        npoints = len(points)
        npms = argv[0][2]
        phX_in_points = argv[0][10]
        elnames = argv[0][11]
        nelements = len(elnames)
        phase_names = argv[0][12]
        changes = argv[1]
    except:
        print('input error')
    else:
        print('correcting indices')
    
    new_phX_in_points = copy.deepcopy(phX_in_points)
    new_npms = copy.deepcopy(npms)
    
    for change in changes:
        phase_to_change, search_element, cutoff = change
        for phase in phase_names:
            if phase_to_change in phase:
                for point in range(npoints):
                    xs_phase_in_point = phX_in_points[phase][point,:]
                    sorted_indices = np.flip(sorted(range(len(phX_in_points[phase][point,:])), key=lambda k: phX_in_points[phase][point,k]))                               
                    sorted_xs_phase_in_point = np.flip(sorted(phX_in_points[phase][point,:]))
                    sorted_elnames = elnames[sorted_indices]
                    search_element_index_sorted = np.where( sorted_elnames== search_element )[0][0]
                    if cutoff > 0 and cutoff<1:
                        if sorted_xs_phase_in_point[search_element_index_sorted] >= cutoff:
                                st=phase_to_change+'-'
                                for nel,el in enumerate(sorted_elnames[:2]):
                                    st += el
                                if st not in new_npms:
                                    new_phX_in_points[st] = np.zeros((npoints,nelements))
                                    new_npms[st] = np.zeros(npoints)
                                    new_phX_in_points[st][point] = xs_phase_in_point
                                    new_npms[st][point] =  new_npms[phase][point]
                                    new_phX_in_points[phase][point] = [0]*len(elnames) 
                                    new_npms[phase][point] = 0
                    elif cutoff == 1:
                        if search_element_index_sorted==0:
                            if sorted_xs_phase_in_point[search_element_index_sorted] >0:
                                st=phase_to_change+'-'
                                for nel,el in enumerate(sorted_elnames[:1]):
                                    st += el
                                if st not in new_npms:
                                    new_phX_in_points[st] = np.zeros((npoints,nelements))
                                    new_npms[st] = np.zeros(npoints)
                                    new_phX_in_points[st][point,:] = xs_phase_in_point
                                    new_npms[st][point] =  new_npms[phase][point]
                                    new_phX_in_points[phase][point] = [0]*len(elnames) 
                                    new_npms[phase][point] = 0
    return(new_phX_in_points,new_npms)
#********************************************************************************************
def ph_namechange(dictionary,name_pairs):
    for name in name_pairs:
        if name[0] in dictionary.keys():
            dictionary[name[1]] = dictionary[name[0]]
            del dictionary[name[0]]
    return dictionary
#********************************************************************************************
def add_phs_compositionDets_DICTRA(*argv):
    npms = argv[0][1]
    phs_names = argv[0][4]
    npms_dict = {}
    for nph,key in enumerate(phs_names):
        npms_dict[key]=npms[:,nph]
    phs_without_Csets =[]
    keys = [key for key in npms_dict.keys()]
    for key in keys:
        if '#1' in key:
            tmp = key.split('#')
            phs_without_Csets.append(tmp[0])
    for ph in phs_without_Csets:
        for i in np.arange(2,10):
            for key in keys:
                 if key == ph+'#{}'.format(i):
                    npms_dict[ph+'#1'] += npms_dict[ph+'#{}'.format(i)]
                    npms_dict.pop(ph+'#{}'.format(i))
    keys = [key for key in npms_dict.keys()]
    for i in np.arange(len(keys)):
        for j in np.arange(i,len(keys)):
             if sorted(keys[i]) == sorted(keys[j]):
                if keys[i] is not keys[j]:
                    print(keys[i],keys[j])
                    npms_dict[keys[i]] += npms_dict[keys[j]]
                    npms_dict.pop(keys[j])
    return npms_dict  
#********************************************************************************************
def add_composition_sets(*argv): 
    npms_dict = copy.deepcopy(argv[0])
    keys = [key for key in npms_dict.keys()]
    phs_without_Csets =[]
    for key in keys:
        if '#1' in key:
            tmp = key.split('#')
            phs_without_Csets.append(tmp[0])
    for ph in phs_without_Csets:
        for i in np.arange(2,10):
            for key in keys:
                 if key == ph+'#{}'.format(i):
                    npms_dict[ph+'#1'] += npms_dict[ph+'#{}'.format(i)]
                    npms_dict.pop(ph+'#{}'.format(i))
    keys = [key for key in npms_dict.keys()]
    for i in np.arange(len(keys)):
        for j in np.arange(i,len(keys)):
             if sorted(keys[i]) == sorted(keys[j]):
                if keys[i] is not keys[j]:
                    print(keys[i],keys[j])
                    npms_dict[keys[i]] += npms_dict[keys[j]]
                    npms_dict.pop(keys[j])
    print(npms_dict.keys())
    return npms_dict
#********************************************************************************************
def plot_dict(*argv):
    '''input: x,y,legend,title,xlims '''
    x = argv[0]
    y = argv[1]
    legend = argv[2]
    title = argv[3]
    xlims = argv[4]
    plt.figure(figsize=[8,6])
    for key in legend:
        plt.plot(x,y[key])
    plt.legend(legend)
    plt.xlim([xlims[0],xlims[1]])
    plt.title(title)
#********************************************************************************************
def plot_list(*argv):
    '''input: x,y,legend,title,xlims '''
    x = argv[0]
    y = argv[1]
    legend = argv[2]
    title = argv[3]
    xlims = argv[4]
    plt.figure(figsize=[8,6])
    for key in legend:
        plt.plot(x,y)
    plt.legend(legend)
    plt.xlim([xlims[0],xlims[1]])
    plt.title(title)
#********************************************************************************************

# def get_lims_gui(points):
#     master = ttk.Tk()
# #     def get_entry_fields():
# #         print("xlim1: %s\nxlim2: %s" % (e1.get(), e2.get()))    
#     ttk.Label(master, text = "xlim1 min {:3.1f} or -1".format(points[0])).grid(row = 0)
#     ttk.Label(master, text = "xlim2 max {:3.1f} or -1".format(points[-1])).grid(row = 1)
#     e1 = ttk.Entry(master)
#     e2 = ttk.Entry(master)
#     e1.grid(row = 0, column = 1)
#     e2.grid(row = 1, column = 1)
#     ttk.Button(master, text = 'OK', command = master.quit).grid(row = 3, column = 0, sticky = ttk.W, pady = 4)
#     #tk.Button(master, text = 'Get', command = get_entry_fields).grid(row = 3, column = 1, sticky = tk.W, pady = 4)
#     ttk.mainloop()
#     master.withdraw()
#     lim1 = int(e1.get())
#     lim2 = int(e2.get())
#     return lim1, lim2
# #********************************************************************************************
# def set_xlim(points1, points3, lim1 = -2, lim2 = -2):
#     '''Setting plot xlimits, L < 0 : domain xlimits'''
#     if lim1 != -2 and lim2 != -2:
#         xlim1, xlim2 = lim1, lim2
#     else:
#         if   max(points1[-1], points3[-1]) > 700: xlim1, xlim2 = 350, 450
#         elif max(points1[-1], points3[-1]) <= 400 and max(points1[-1], points3[-1]) >= 200 : xlim1, xlim2 = 150, 250
#         elif max(points1[-1], points3[-1]) >= 250 and max(points1[-1], points3[-1]) <= 360 : xlim1, xlim2 = 280, max(points1[-1], points3[-1])
#         elif max(points1[-1], points3[-1]) <= 150 and max(points1[-1], points3[-1]) >= 100 : xlim1, xlim2 = 80, max(points1[-1], points3[-1])
#         elif max(points1[-1], points3[-1]) <= 100 : xlim1, xlim2 = 40, max(points1[-1], points3[-1])
#         else : xlim1, xlim2 = -1, -1
#     return xlim1, xlim2
# x1,x2 = get_lims_gui(tStamp_vals[3])
#********************************************************************************************
def del_pngs():
    '''Deleting png files'''
    fileList = glob.glob('*.png', recursive = True)
    # Iterate over the list of filepaths & remove each file.
    for filePath in fileList:
        try:
            os.remove(filePath)
        except OSError:
            print("Error while deleting file")  
#********************************************************************************************
def all_x_plotter(names, L, tt, xs, lim1 = -10, lim2 = -10, title = '', file_name = 'all_values', ylab = "fractions"\
                  , fsize = 20, lsize = 20, tsize = 20, bins = 6, tksize = 20, stsize=20, lgsize=20, wlsize=3):
    new_names=[]
    for name in names:
        if '_' in name:
            name = name.replace('_','\_')
        if '#' in name:
            name = name.replace('#','\#')
        name = '$'+name+'$'
        new_names.append(name)
    if lim1 < 0: lim1 = L[0]
    if lim2 < 0: lim2 = L[-1]
    plt.figure(figsize = [15, 10])
    for i, x in enumerate(names):
        #ysmoothed = gaussian_filter1d(xs[:, i], sigma = 0.1)
        plt.plot(L, xs[:, i], linewidth=wlsize)
        plt.legend(new_names, fontsize = lgsize)
    plt.title(title, fontsize = tsize)
    plt.xlabel('Distance [um]', fontsize = lsize)
    plt.ylabel(ylab, fontsize = lsize)
    plt.xticks(fontsize = tksize)
    plt.yticks(fontsize = tksize)
    plt.xlim([lim1, lim2])
    plt.savefig(file_name+'.png')
#********************************************************************************************
def plot_x_3ts(x1, x2, x3, L1, L2, L3, tt1, tt2, tt3, names, lim1 = -10, lim2 = -10, title = '', 
               file_name = '3ts', ylab = "fractions",  
               fsize = 20, lsize = 20, tsize = 20, bins = 6, tksize = 20, stsize=20, lgsize=20):  
    if lim1 < 0: lim1 = min(L1[0], L2[0])
    if lim2 < 0: lim2 = max(L1[-1], L2[-1])
    if len(names)//2 > 1:
        if len(names)%2 == 1:
            fig1, ax1 = plt.subplots(len(names)//2+1, 2, figsize = [15, (len(names)//2+1)*25/5])
        else:
            fig1, ax1 = plt.subplots(len(names)//2, 2, figsize = [15, len(names)//2*25/5])
        ie = 0
        for el in names:
            ysmoothed = gaussian_filter1d(x1[:, ie], sigma = 0.1)
            ax1[ie//2, ie%2].plot(L1, ysmoothed, '.-', label = 't = {:2.1e} sec'.format(tt1) )
            ysmoothed = gaussian_filter1d(x2[:, ie], sigma = 0.1)
            ax1[ie//2, ie%2].plot(L2, ysmoothed, '.-', label = 't = {:2.1e} sec'.format(tt2) )
            ysmoothed = gaussian_filter1d(x3[:, ie], sigma = 0.1)
            ax1[ie//2, ie%2].plot(L3, ysmoothed, '.-', label = 't = {:2.1e} sec'.format(tt3) )
            ax1[ie//2, ie%2].legend(fontsize = lgsize)
            ax1[ie//2, ie%2].set_title(el, fontsize = tsize)
            ax1[ie//2, ie%2].set_xlabel('Distance [um]', fontsize = lsize)
            ax1[ie//2, ie%2].set_ylabel(ylab, fontsize = lsize)
            ax1[ie//2, ie%2].tick_params(axis = "x", labelsize = tksize)
            ax1[ie//2, ie%2].tick_params(axis = "y", labelsize = tksize)
            ax1[ie//2, ie%2].set_xlim([lim1, lim2])
            ax1[ie//2, ie%2].locator_params(axis = 'y', nbins = bins)
            ax1[ie//2, ie%2].locator_params(axis = 'x', nbins = bins)
            ie += 1
        fig1.suptitle(title, fontsize = stsize)
        fig1.tight_layout()
        plt.savefig(file_name+'.png')
    else:
        fig1, ax1 = plt.subplots(len(names), 1, figsize = [15, len(names)*25/5])
        ie = 0
        for el in names:
            ysmoothed = gaussian_filter1d(x1[:, ie], sigma = 0.1)
            ax1[ie].plot(L1, ysmoothed, '.-', label = 't = {:2.1e} sec'.format(tt1))
            ysmoothed = gaussian_filter1d(x2[:, ie], sigma = 0.1)
            ax1[ie].plot(L2, ysmoothed, '.-', label = 't = {:2.1e} sec'.format(tt2))
            ysmoothed = gaussian_filter1d(x3[:, ie], sigma = 0.1)
            ax1[ie//2, ie%2].plot(L3, ysmoothed, '.-', label = 't = {:2.1e} sec'.format(tt3))
            ax1[ie].legend(fontsize = lgsize)
            ax1[ie].set_title(el, fontsize = tsize)
            ax1[ie].set_xlabel('Distance [um]', fontsize = lsize)
            ax1[ie].set_ylabel(' ', fontsize = lsize)
            ax1[ie].tick_params(axis = "x", labelsize = tksize) 
            ax1[ie].tick_params(axis = "y", labelsize = tksize)
            ax1[ie].locator_params(axis = 'y', nbins = bins)
            ax1[ie].locator_params(axis = 'x', nbins = bins)       
            ax1[ie].set_xlim([lim1, lim2])
            ie += 1   
        fig1.suptitle(title, fontsize = stsize)
        fig1.tight_layout()
        plt.savefig(title+'.png')
#********************************************************************************************
def plot_tc_corrected_results_phf_all_phases(npm, xph, t, TIM, L, lim1 = -10, lim2 = -10
                                             , fsize = 20, lsize = 20, tsize = 20, bins = 6, tksize = 20, stsize=20, lgsize=20, lwsize=3, filename='NPM-corrected-postprocessed_'):
    if lim1 <= 0: lim1 = L[0]
    if lim2 <= 0: lim2 = L[-1]
    plt.figure(figsize = [15, 15])
    for ph in npm.keys():
        if len(npm[ph]) != 0:
            #ysmoothed = gaussian_filter1d(npm[ph], sigma = 0.1)
            plt.plot(xph[ph][0], npm[ph], linewidth=lwsize)
    plt.legend(npm.keys(), fontsize = lgsize)
    #plt.title(' Corrected NPM t = {:4.0f}sec, Post processed with all phases in the database'.format(TIM[t]), fontsize = 20)
    plt.xlabel(r'Distance ($\mu m$)', fontsize = lsize)
    plt.ylabel('Phase fraction', fontsize = lsize)
    plt.xticks(fontsize = tksize)
    plt.yticks(fontsize = tksize)
    plt.xlim([lim1, lim2])
    plt.savefig(filename+'-{:4.0f}-sec.png'.format(TIM[t]))
#********************************************************************************************
def plot_zero_padded_npm_corrected(Lx2, npm_pp_corrected, ttx2, lim1, lim2, L
                                   , fsize = 20, lsize = 20, tsize = 20, bins = 6, tksize = 20, stsize=20, lgsize=20, lwsize=3, filename='NPM-zeroPadded-Corrected'):
    if lim1 <= 0: lim1 = L[0]
    if lim2 <= 0: lim2 = L[-1]
    r = np.zeros(len(Lx2)+1)
    for i in range(1, len(Lx2)+1):
        r[i] = (Lx2[i-1]-r[i-1])+Lx2[i-1]
    A = []
    for ph in npm_pp_corrected.keys():
        A.append(ph)
    plt.figure(figsize = [15, 10])
    for ph in A:
        vertices_npm = np.zeros(len(r))
        vertices_npm[0] = npm_pp_corrected[ph][0]
        for i in range(1, len(vertices_npm)-1):
                vertices_npm[i] = (npm_pp_corrected[ph][i]+npm_pp_corrected[ph][i-1])/2
        vertices_npm[-1] = npm_pp_corrected[ph][-1]
        ysmoothed = vertices_npm#gaussian_filter1d(vertices_npm, sigma = 0.1)
        plt.plot(r, ysmoothed, linewidth=lwsize)
    plt.legend(A, fontsize = lgsize)
    plt.xlabel(r'Distance ($\mu m$)', fontsize = lsize)
    plt.ylabel('Phase fraction', fontsize = lsize)
    plt.xticks(fontsize = tksize)
    plt.yticks(fontsize = tksize)
    plt.xlim([lim1, lim2])
    #plt.title('zero padded corrected Phase fraction t = {:2.1e} sec'.format(ttx2), fontsize = 20)
    plt.savefig(filename+'-{:4.0f}-sec.png'.format(ttx2))
#********************************************************************************************
def plot_x_3ts_log10(x1, x2, x3, L1, L2, L3, tt1, tt2, tt3, names, lim1 = -10, lim2 = -10, title = '', 
               file_name = '3ts', ylab = "fractions",  fsize = 20, lsize = 20, tsize = 20, bins = 6, 
               tksize = 20, stsize=20, lgsize=20, lwsize=3):
    params = {'text.usetex': True}
    plt.rcParams.update(params)
    
    if lim1 < 0: lim1 = min(L1[0], L2[0])
    if lim2 < 0: lim2 = max(L1[-1], L2[-1])
    if len(names)//2 > 1:
        if len(names)%2 == 1:
            fig1, ax1 = plt.subplots(len(names)//2+1, 2, figsize = [15, (len(names)//2+1)*25/5])
        else:
            fig1, ax1 = plt.subplots(len(names)//2, 2, figsize = [15, len(names)//2*25/5])
        ie = 0
        for el in names:
            ax1[ie//2, ie%2].plot(L1, -1*np.log10(np.abs(x1[:, ie])), label = 't = {:4.0f} sec'.format(tt1) , linewidth=lwsize)
            ax1[ie//2, ie%2].plot(L2, -1*np.log10(np.abs(x2[:, ie])), label = 't = {:4.0f} sec'.format(tt2) , linewidth=lwsize)
            ax1[ie//2, ie%2].plot(L3, -1*np.log10(np.abs(x3[:, ie])), label = 't = {:4.0f} sec'.format(tt3) , linewidth=lwsize)
            ax1[ie//2, ie%2].legend(fontsize = lgsize)
            ax1[ie//2, ie%2].set_title(el, fontsize = tsize)
            ax1[ie//2, ie%2].set_xlabel(r'Distance ($\mu m$)', fontsize = lsize)
            ax1[ie//2, ie%2].set_ylabel(ylab, fontsize = lsize)
            ax1[ie//2, ie%2].tick_params(axis = "x", labelsize = tksize)
            ax1[ie//2, ie%2].tick_params(axis = "y", labelsize = tksize)
            ax1[ie//2, ie%2].set_xlim([lim1, lim2])
            ax1[ie//2, ie%2].locator_params(axis = 'y', nbins = bins)
            ax1[ie//2, ie%2].locator_params(axis = 'x', nbins = bins)
            ie += 1
        fig1.suptitle(title, fontsize = stsize)
        fig1.tight_layout()
        plt.savefig(file_name+'.png')
    else:
        fig1, ax1 = plt.subplots(len(names), 1, figsize = [15, len(names)*25/5])
        ie = 0
        for el in names:
            ax1[ie].plot(L1, -1*np.log10(np.abs(x1[:, ie])), label = 't = {:4.0f} sec'.format(tt1), linewidth=lwsize)
            ax1[ie].plot(L2, -1*np.log10(np.abs(x2[:, ie])), label = 't = {:4.0f} sec'.format(tt2), linewidth=lwsize)
            ax1[ie].plot(L3, -1*np.log10(np.abs(x3[:, ie])), label = 't = {:4.0f} sec'.format(tt3), linewidth=lwsize)
            ax1[ie].legend(fontsize = lgsize)
            ax1[ie].set_title(el, fontsize = tsize)
            ax1[ie].set_xlabel(r'Distance ($\mu m$) ', fontsize = lsize)
            ax1[ie].set_ylabel(' ', fontsize = lsize)
            ax1[ie].tick_params(axis = "x", labelsize = tksize) 
            ax1[ie].tick_params(axis = "y", labelsize = tksize)
            ax1[ie].locator_params(axis = 'y', nbins = bins)
            ax1[ie].locator_params(axis = 'x', nbins = bins)       
            ax1[ie].set_xlim([lim1, lim2])
            ie += 1   
        fig1.suptitle(title, fontsize = stsize)
        fig1.tight_layout()
        plt.savefig(file_name+'.png')       
#********************************************************************************************      
def all_x_plotter_log10(names, L, tt, xs, lim1 = -10, lim2 = -10, title = '', file_name = 'all_values', ylab = "fractions"\
                  , fsize = 20, lsize = 20, tsize = 20, bins = 6, tksize = 20, stsize=20, lgsize=20, wlsize=3):
    
    if lim1 < 0: lim1 = L[0]
    if lim2 < 0: lim2 = L[-1]
    plt.figure(figsize = [15, 10])
    for i, x in enumerate(names):
        #ysmoothed = gaussian_filter1d(xs[:, i], sigma = 0.1)
        plt.plot(L, np.log10(xs[:, i]), linewidth=wlsize)
        plt.legend(names, fontsize = lgsize)
    plt.title(title, fontsize = tsize)
    plt.xlabel('Distance [um]', fontsize = lsize)
    plt.ylabel(ylab, fontsize = lsize)
    plt.xticks(fontsize = tksize)
    plt.yticks(fontsize = tksize)
    plt.xlim([lim1, lim2])
    plt.savefig(file_name+'.png')   
#********************************************************************************************
def plot_tc_results_xphf_phases(xph, allnpm, t, elnames, ph_names, L, lim1 = -10, lim2 = -10
                                , fsize = 20, lsize = 20, tsize = 20, bins = 6, tksize = 20, stsize=20, lgsize=20):
    '''grid point = data[n]
    phase names = data[n].keys() 
    phase fraction = data[n][keys][0]
    el fractions in phase = data[n][keys][1]
    ordered el names in phase = data[n][keys][2]
    ordered el fractions in phase = data[n][keys][3]'''
    params = {'text.usetex': False}
    plt.rcParams.update(params)
    if lim1 <= 0: lim1 = L[0]
    if lim2 <= 0: lim2 = L[-1]
    if len(ph_names)//2 > 1:
        if len(ph_names)%2 == 1:
            fig1, ax1 = plt.subplots(len(ph_names)//2+1, 2, figsize = [15, (len(ph_names)//2+1)*25/5])
        else:
            fig1, ax1 = plt.subplots(len(ph_names)//2, 2, figsize = [15, len(ph_names)//2*25/5])
        iph = 0
        for iph,ph in enumerate(ph_names):  
            for ie, el in enumerate(elnames):  
                ax1[iph//2, iph%2].plot(xph[ph][0], xph[ph][1][:, ie]*np.array(allnpm[ph]), 'o-', label = el)
            ax1[iph//2, iph%2].plot(xph[ph][0], allnpm[ph], 'silver', marker = 'd', label = 'npm {}'.format(ph))
            ax1[iph//2, iph%2].legend(fontsize = lgsize)
            ax1[iph//2, iph%2].set_title('X in {} t = {:2.1e}sec'.format(ph, t), fontsize = tsize)
            ax1[iph//2, iph%2].set_xlabel('Distance [um]', fontsize = lsize)
            ax1[iph//2, iph%2].set_ylabel('Fraction', fontsize = lsize)
            ax1[iph//2, iph%2].tick_params(axis = "x", labelsize = tksize) 
            ax1[iph//2, iph%2].tick_params(axis = "y", labelsize = tksize)
            ax1[iph//2, iph%2].locator_params(axis = 'y', nbins = bins)
            ax1[iph//2, iph%2].locator_params(axis = 'x', nbins = bins)
            ax1[iph//2, iph%2].set_xlim([lim1, lim2])
            iph += 1
        fig1.tight_layout()
        plt.savefig('XinPh-postprocessed-{:.0f}-sec.png'.format(t))  
    elif len(new_ph_names) == 1:
        ph = ph_names[0]
        plt.figure(figsize = [15, len(ph_names)*25/5])
        for ie, el in enumerate(elnames):
            plt.plot(xph[ph][0], xph[ph][1][:, ie]*np.array(allnpm[ph]), 'o-', label = el)
        plt.plot(xph[ph][0], allnpm[ph], 'silver', marker = 'd', label = 'npm {}'.format(ph))
        plt.legend(fontsize = lgsize)        
        plt.title('xs in {} t = {:2.1e}sec'.format(ph, TIM[t]), fontsize = tsize)
        plt.xlabel('Distance [um]', fontsize = lsize)
        plt.ylabel('Fractions', fontsize = lsize)
        plt.tick_params(axis = "x", labelsize = tksize) 
        plt.tick_params(axis = "y", labelsize = tksize)
        plt.locator_params(axis = 'y', nbins = bins)
        plt.locator_params(axis = 'x', nbins = bins)
        plt.xlim([lim1, lim2])
    else:
        iph = 0
        fig1, ax1 = plt.subplots(len(ph_names), 1, figsize = [15, len(ph_names)*25/5])
        for iph,ph in enumerate(ph_names):
            for ie, el in enumerate(elnames):
                ax1[iph].plot(xph[ph][0], xph[ph][1][:, ie]*np.array(allnpm[ph]), 'o-', label = el)
            ax1[iph].plot(xph[ph][0], allnpm[ph][1], 'silver', marker = 'd', label = 'npm {}'.format(ph))
            ax1[iph].legend(fontsize = lgsize)            
            ax1[iph].set_title('xs in {} t = {:2.1e}sec'.format(ph, t), fontsize = tsize)
            ax1[iph].set_xlabel('Distance [um]', fontsize = lsize)
            ax1[iph].set_ylabel('Fractions', fontsize = lsize)
            ax1[iph].tick_params(axis = "x", labelsize = tksize) 
            ax1[iph].tick_params(axis = "y", labelsize = tksize)
            ax1[iph].locator_params(axis = 'y', nbins = bins)
            ax1[iph].locator_params(axis = 'x', nbins = bins)
            ax1[iph].set_xlim([lim1, lim2])
            iph += 1
        fig1.tight_layout()
        plt.savefig('XinPh-postprocessed-{:4.0f}-sec.png'.format(t))        
#********************************************************************************************        
def plot_allx_2ts(x1, x2, L1, L2, tt1, tt2, names, lim1 = -10, lim2 = -10, title = '', file_name = '', fsize = 20, lsize = 20, tsize = 20, bins = 6, tksize = 20, stsize=20, lgsize=20, lwsize=20):
    NUM_COLORS = len(names)*2
    r = list(range(NUM_COLORS))
    random.shuffle(r)
    cm = plt.get_cmap('gist_rainbow')
    fig = plt.figure(figsize = [15, 10])
    ax = fig.add_subplot(111)
    ax.set_prop_cycle('color', [cm(1.*i/NUM_COLORS) for i in r])

    if lim1 < 0: lim1 = min(L1[0], L2[0])
    if lim2 < 0: lim2 = max(L1[-1], L2[-1])
    ie = 0
    for el in names:
        ax.plot(L1, x1[:, ie], linewidth=lwsize, label = '{}, {:.0f} sec'.format(el,tt2) )
        ax.plot(L2, x2[:, ie],'--', linewidth=lwsize, label = '{}, {:.0f} sec'.format(el,tt1) )
        ax.set_xlabel('Distance [um]', fontsize = lsize)
        ax.set_ylabel('Mole Fractions', fontsize = lsize)
        ax.tick_params(axis = "x", labelsize = tksize)
        ax.tick_params(axis = "y", labelsize = tksize)
        ax.set_xlim([lim1, lim2])
        ax.locator_params(axis = 'y', nbins = bins)
        ax.locator_params(axis = 'x', nbins = bins)
        ie += 1
    ax.legend(fontsize = lgsize)
    plt.savefig(file_name+'.png')
