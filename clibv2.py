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
import pdb
import random
import sys
from collections import defaultdict
#********************************************************************************************
def get_values_from_textfiles(*argv):
    '''Reads all values from a given path, 
    input:
        path(optional)
    output:dict{path, all_mfs, elnames, all_pts, times, n_pts, DICT_all_npms, 
                 DICT_phnames, DICT_all_mus, T, interstitials, substitutionals}
    '''
    data = {}
    if len(argv) == 1:
        data['path'] = argv[0]
        os.chdir(data['path'])
    else:
        data['path'] = get_path()
    data['all_mfs'] = np.loadtxt('MOLE_FRACTIONS.TXT', dtype = float)
    data['elnames'] = np.loadtxt('EL_NAMES.TXT', dtype = str)
    data['all_pts'] = np.loadtxt('VOLUME_MIDPOINTS.TXT', dtype = float)
    data['times'] = np.loadtxt('TIME.TXT', dtype = float)
    data['n_pts'] = np.loadtxt('VOLUMES_PER_REGION.TXT', dtype = float)
    data['DICT_all_npms'] = np.loadtxt('PHASE_FRACTIONS.TXT', dtype = float)
    phnames = []
    with open('PH_NAMES.TXT', 'r')as f:
        for line in f:
            phnames.append(line.strip())
    data['DICT_phnames'] = np.array(phnames)
    data['DICT_all_mus'] = np.loadtxt('CHEMICAL_POTENTIALS.TXT', dtype = float)
    data['T'] = np.loadtxt('T.DAT', dtype = int)+273
    int_idx = []
    if 'N' in data['elnames']:
        Nidx = np.where(data['elnames']  ==  'N')[0]
        int_idx.append(Nidx[0])
    if 'C' in data['elnames']:
        Cidx = np.where(data['elnames']  ==  'C')[0]
        int_idx.append(Cidx[0])
    if 'H' in data['elnames']:
        Hidx = np.where(data['elnames']  ==  'H')[0]
        int_idx.append(Hidx[0])
    if 'O' in data['elnames']:
        Oidx = np.where(data['elnames']  ==  'O')[0]
        int_idx.append(Oidx[0])
    if 'VA' in data['elnames']:
        VAidx = np.where(data['elnames']  ==  'VA')[0]
        int_idx.append(VAidx[0])
    sub_idx = list(np.arange(len(data['elnames'])))
    data['interstitials'] = [data['elnames'][int_idx], int_idx]
    for idx in data['interstitials'][1]: 
        sub_idx.pop(idx)
    data['substitutionals'] = [data['elnames'][sub_idx], sub_idx]
    return data
#********************************************************************************************
def get_path():
    '''change path to the designated folder, root is users home
    output:
        path
    '''
    flag = True
    while flag:
        root = Tk() # pting root to Tk() to use it as Tk() in program.
        root.withdraw() # Hides small tkinter window.
        root.attributes('-topmost', True) # Opened windows will be active. above all windows despite of selection.
        open_file = filedialog.askdirectory(initialdir = '~') # Returns opened path as str
        if open_file:
            flag = False
    os.chdir(open_file)
    return os.getcwd()
#********************************************************************************************
def get_timeStamp(times):
    time1 = float(input('Time to plot/{}/timesteps/{}/:'.format(times[-1], len(times))))
    if time1>= 0 and time1 <= times[-1]:
        tS, nearestTime = find_nearest(times, time1)
    return tS, nearestTime
#********************************************************************************************                      
def find_nearest(array, value):
    '''output: index, nearest value'''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]                      
#********************************************************************************************                      
def calculate_u_fractions(*argv):
    mf, sub_idx, elnames =  argv
    sub_sum = np.sum(mf[:, sub_idx], axis = 1)
    #uf_dict = {}
    uf = []
    for nel, el in enumerate(elnames):
        #uf_dict[el] = mf[:, nel]/sub_sum[:]
        uf.append(mf[:, nel]/sub_sum[:]) 
    return np.array(uf)
#********************************************************************************************
def get_tS_VLUs(dict_input, tS, nearestTime):
    '''get value of a given timestamp
    input: dict{all_values}, tS 
    output: dict{'tS_mfs', 'tS_DICT_npms', 'tS_DICT_mus', 'tS_pts', 'tS_DICT_phnames'}
    '''
    dict = copy.deepcopy(dict_input)
    dict['tS'] = tS
    dict['nearestTime'] = nearestTime
    mfs  = dict['all_mfs']
    elnames = dict['elnames']
    pts = dict['all_pts']
    n_pts  = dict['n_pts']
    npms = dict['DICT_all_npms']
    phnames = dict['DICT_phnames']
    mus = dict['DICT_all_mus']
    subs  = dict['substitutionals']   
    '''gets molefractions, chemicalpotentials, phasefractions, volumemidpts and time of the timeindex
    removes ZZDICT from phasenames, calculate u-fractions'''
    dict['tS_pts'] = pts[int(np.sum(n_pts[:tS])):int(np.sum(n_pts[:tS+1]))]*1e6    
    idx1 = int(np.sum(n_pts[:tS]))*int(len(phnames))
    idx2 = int(np.sum(n_pts[:tS+1]))*int(len(phnames))
    dict['tS_DICT_npms'] = npms[idx1:idx2].reshape((-1, int(len(phnames))))
    if 'ZZDICT_GHOST' in phnames[-1]:
        dict['tS_DICT_phnames'] = deepcopy(phnames[:-1])
    else:
        dict['tS_DICT_phnames'] = deepcopy(phnames)
    idx1 = int(np.sum(n_pts[:tS]))*len(elnames)
    idx2 = int(np.sum(n_pts[:tS+1]))*len(elnames)
    dict['tS_DICT_mfs'] = mfs[idx1:idx2].reshape((-1, len(elnames)))
    dict['tS_DICT_mus'] = mus[idx1:idx2].reshape((-1, len(elnames)))
    dict['tS_DICT_ufs'] = calculate_u_fractions(dict['tS_DICT_mfs'], subs[1], elnames).T 
    for key in ['all_mfs', 'all_pts', 'times', 'n_pts', 'DICT_all_npms' , 'DICT_phnames', 'DICT_all_mus']:
        dict.pop(key) 
    return dict
#********************************************************************************************
def tccalc(dict_input):
    ''' Calculate single equilibrium pt by pt with input condition
        input: dict{'path', 'elnames', 'T', 'interstitials', 'substitutionals', 
                    'tS', 'nearestTime', 'tS_pts', 'tS_DICT_npms', 
                    'tS_DICT_phnames', 'tS_mfs', 'tS_mus', 'tS_ufs', 
                    'xlim1', 'xlim2', 'phase_change', 'tc_Setting', 'name_paires'}
        output: dict{tS_TC_phnames, tS_TC_npms, tS_TC_vpvs, tS_TC_phXs, 
            tS_TC_acRef, tStmp_TC_acSER, tS_TC_mus, tS_TC_ws}
    '''
    try:
        dict = copy.deepcopy(dict_input)
        mfs = dict['tS_DICT_mfs']
        database = dict['tc_setting']['database']
        n_pts = len(dict['tS_pts'])
        phsToSus = dict['tc_setting']['phsToSus']
        acsRef = dict['tc_setting']['acRefs']
        T = dict['T']
        elnames = dict['elnames']
        pth = dict['path']
        p3flag = dict['tc_setting']['p3flag']
    except:
        print('error in TC input data')
        sys.exit()
    else:
        print('calculating thermodynaics')
    
    with TCPython() as start:
        if p3flag and os.path.isfile('{}/p.POLY3'.format(pth)):
            poly = start.select_database_and_elements(database, elnames).get_system().with_single_equilibrium_calculation()
            poly.run_poly_command('read {}/p.POLY3'.format(pth))
            #poly.run_poly_command('list-status , cps, ')
            poly.remove_all_conditions()
            poly.set_condition('N', 1)
            poly.set_condition('P', 1e5)
        else:
            print('read poly3 file, calculating')
            poly = start.set_cache_folder( "./_cache2").select_database_and_elements(database, elnames).get_system().with_single_equilibrium_calculation()
        if len(phsToSus)>0:
            for phase in phsToSus:
                poly.set_phase_to_suspended(phase)
        poly.set_condition("T", T)
        tc_phnames, tc_npms, tc_vpvs, tc_phXs = {}, {}, {}, {}
        tc_acRefs, tc_acSER, tc_mus, tc_ws = defaultdict(list), defaultdict(list), defaultdict(list), defaultdict(list)
        for pt in range(n_pts):
            for nel, el in enumerate(elnames[:-1]):
                poly.set_condition("X({})".format(el), mfs[pt, nel])
            pntEq = poly.calculate()
            stablePhs = pntEq.get_stable_phases()
            tc_phnames[pt] = stablePhs        
            if 'C' in elnames:
                for reference in acsRef:
                    tc_acRefs['ac(C, '+reference+')'].append(pntEq.get_value_of('ac(C, {})'.format(reference)))  
            for nel, el in enumerate(elnames[:]):
                tc_acSER[el].append(pntEq.get_value_of('ac({})'.format(el)))
                tc_mus[el].append(pntEq.get_value_of('mu({})'.format(el)))
                tc_ws[el].append(pntEq.get_value_of('w({})'.format(el)))
            for ph in stablePhs:      
                tc_npms['{}, {}'.format(pt, ph)] = pntEq.get_value_of('npm({})'.format(ph))
                tc_vpvs['{}, {}'.format(pt, ph)] = pntEq.get_value_of('vpv({})'.format(ph))
                temp1 = []
                for el2 in elnames[:]:
                    temp1.append(pntEq.get_value_of('X({}, {})'.format(ph, el2)))
                tc_phXs['{}, {}'.format(pt, ph)] = np.array(temp1)
    clear_output(wait = False)
    dict['tS_TC_phnames'] = tc_phnames
    dict['tS_TC_npms'] = tc_npms
    dict['tS_TC_vpvs'] = tc_vpvs
    dict['tS_TC_phXs'] = tc_phXs
    dict['tS_TC_acRef'] = tc_acRefs
    dict['tS_TC_acSER'] = tc_acSER
    dict['tS_TC_mus'] = tc_mus
    dict['tS_TC_ws'] = tc_ws
    ### trimming
    tc_NEAT_mfs, tc_NEAT_npms, tc_NEAT_vpvs, tc_NEAT_phXs = {}, {}, {}, {}
    tc_NEAT_phnames = []
    for nel, el in enumerate(elnames):
        tc_NEAT_mfs[el] = mfs[:, nel]
    for i in range(n_pts):
        for j in tc_phnames[i]:
               if j not in tc_NEAT_phnames:
                    tc_NEAT_phnames.append(j)
    for ph in tc_NEAT_phnames:
        tc_NEAT_npms[ph] = np.zeros([n_pts])
        tc_NEAT_vpvs[ph] = np.zeros([n_pts])
        tc_NEAT_phXs[ph] = np.zeros([n_pts, len(elnames)])
    for nph, ph in enumerate(tc_NEAT_phnames):
        for pt in range(n_pts):
            if '{}, {}'.format(str(pt), ph) in tc_npms:
                tc_NEAT_npms[ph][pt] = tc_npms['{}, {}'.format(str(pt), ph)]
                tc_NEAT_vpvs[ph][pt] = tc_vpvs['{}, {}'.format(str(pt), ph)]    
                for nel, el in enumerate(elnames):
                    tc_NEAT_phXs[ph][pt, nel] = tc_phXs['{}, {}'.format(str(pt), ph)][nel]
    dict['tS_TC_NEAT_phnames'] = tc_NEAT_phnames
    dict['tS_TC_NEAT_npms'] = tc_NEAT_npms
    dict['tS_TC_NEAT_vpvs'] = tc_NEAT_vpvs
    dict['tS_TC_NEAT_phXs'] = tc_NEAT_phXs
    dict['tS_TC_NEAT_mfs'] = tc_NEAT_mfs
    dict1 = correct_phase_indices(dict)
    dict2 = add_compSets(dict1)
    dict3 = phnameChange(dict2)
    dict_out = add_compSets_DICT(dict3)
    return dict_out
# #********************************************************************************************
def correct_phase_indices(dict_input):
    '''input:{'elnames', 'tS_TC_NEAT_phnames', 'tS_TC_NEAT_npms', 
    'tS_TC_NEAT_vpvs', 'tS_TC_NEAT_phXs} ''' 
    try: 
        dict = copy.deepcopy(dict_input)
        n_pts = len(dict['tS_pts']) 
        tS_TC_NEAT_npms = dict['tS_TC_NEAT_npms'] 
        tS_TC_NEAT_phXs = dict['tS_TC_NEAT_phXs'] 
        elnames = dict['elnames'] 
        nelements = len(elnames) 
        tS_TC_NEAT_phnames = dict['tS_TC_NEAT_phnames']
        phase_changes = dict['phase_changes'] 
    except:                     
        print('input error, correcting phase indices')
        sys.exit()    
    else:                       
        print('correcting indices')
                                
    CQT_tS_TC_NEAT_phXs = copy.deepcopy(tS_TC_NEAT_phXs)
    CQT_tS_TC_NEAT_npms = copy.deepcopy(tS_TC_NEAT_npms)
    for change in phase_changes:
        phase_to_change, search_element, cutoff = change
        for phase in tS_TC_NEAT_phnames:
            if phase_to_change in phase:
                for pt in range(n_pts):
                    phXs = tS_TC_NEAT_phXs[phase][pt, :]
                    sorted_indices = np.flip(sorted(range(len(tS_TC_NEAT_phXs[phase][pt, :])), key = lambda k: tS_TC_NEAT_phXs[phase][pt, k]))                               
                    sorted_phXs = np.flip(sorted(tS_TC_NEAT_phXs[phase][pt, :]))
                    sorted_elnames = elnames[sorted_indices]
                    sorted_searchElidx = np.where( sorted_elnames ==  search_element )[0][0]
                    if cutoff > 0 and cutoff<1:
                        if sorted_phXs[sorted_searchElidx] >  cutoff:
                                st = phase_to_change+'-'
                                for nel, el in enumerate(sorted_elnames[:2]):
                                    st +=  el
                                if st not in CQT_tS_TC_NEAT_npms:
                                    CQT_tS_TC_NEAT_phXs[st] = np.zeros((n_pts, nelements))
                                    CQT_tS_TC_NEAT_npms[st] = np.zeros(n_pts)
                                    CQT_tS_TC_NEAT_phXs[st][pt] = phXs
                                    CQT_tS_TC_NEAT_npms[st][pt] =  CQT_tS_TC_NEAT_npms[phase][pt]
                                    CQT_tS_TC_NEAT_phXs[phase][pt] = [0]*nelements 
                                    CQT_tS_TC_NEAT_npms[phase][pt] = 0
                    elif cutoff  ==  1:
                        if sorted_searchElidx == 0:
                            if sorted_phXs[sorted_searchElidx] >0:
                                st = phase_to_change+'-'
                                for nel, el in enumerate(sorted_elnames[:1]):
                                    st +=  el
                                if st not in CQT_tS_TC_NEAT_npms:
                                    CQT_tS_TC_NEAT_phXs[st] = np.zeros((n_pts, nelements))
                                    CQT_tS_TC_NEAT_npms[st] = np.zeros(n_pts)
                                    CQT_tS_TC_NEAT_phXs[st][pt, :] = phXs
                                    CQT_tS_TC_NEAT_npms[st][pt] =  CQT_tS_TC_NEAT_npms[phase][pt]
                                    CQT_tS_TC_NEAT_phXs[phase][pt] = [0]*nelements 
                                    CQT_tS_TC_NEAT_npms[phase][pt] = 0
    dict['CQT_tS_TC_NEAT_phXs'] = CQT_tS_TC_NEAT_phXs 
    dict['CQT_tS_TC_NEAT_npms'] = CQT_tS_TC_NEAT_npms 
    return dict
#********************************************************************************************
def add_compSets(dict_in):
    dict = copy.deepcopy(dict_in) 
    npms_dict = dict['CQT_tS_TC_NEAT_npms'] 
    keys = [key for key in npms_dict.keys()]
    phs_without_Csets  = []
    for key in keys:
        if '#1' in key:
            tmp = key.split('#')
            phs_without_Csets.append(tmp[0])
    for ph in phs_without_Csets:
        for i in np.arange(2, 10):
            for key in keys:
                 if key  ==  ph+'#{}'.format(i):
                    npms_dict[ph+'#1'] +=  npms_dict[ph+'#{}'.format(i)]
                    npms_dict.pop(ph+'#{}'.format(i))
    keys = [key for key in npms_dict.keys()]
    for i in np.arange(len(keys)):
        for j in np.arange(i, len(keys)):
             if sorted(keys[i])  ==  sorted(keys[j]):
                if keys[i] is not keys[j]:
                    print(keys[i], keys[j])
                    npms_dict[keys[i]] +=  npms_dict[keys[j]]
                    npms_dict.pop(keys[j])
    dict['sum_CQT_tS_TC_NEAT_npms'] = npms_dict
    return dict
#********************************************************************************************
def phnameChange(dict_in):
    dict = copy.deepcopy(dict_in)
    name_pairs = dict['name_pairs']
    npms_dict = dict['CQT_tS_TC_NEAT_npms'] 
    for name in name_pairs:
        if name[0] in npms_dict.keys():
            npms_dict[name[1]] = npms_dict[name[0]]
            del npms_dict[name[0]]
    dict['nameChanged_CQT_tS_TC_NEAT_npms'] = npms_dict
    return dict
#********************************************************************************************
def add_compSets_DICT(dict_in):
    dict = copy.deepcopy(dict_in)
    npms = dict['tS_DICT_npms'] 
    phs_names = dict['tS_DICT_phnames'] 
    npms_dict = {}
    for nph, key in enumerate(phs_names):
        npms_dict[key] = npms[:, nph]
    phs_without_Csets  = []
    keys = [key for key in npms_dict.keys()]
    for key in keys:
        if '#1' in key:
            tmp = key.split('#')
            phs_without_Csets.append(tmp[0])
    for ph in phs_without_Csets:
        for i in np.arange(2, 10):
            for key in keys:
                 if key  ==  ph+'#{}'.format(i):
                    npms_dict[ph+'#1'] +=  npms_dict[ph+'#{}'.format(i)]
                    npms_dict.pop(ph+'#{}'.format(i))
    keys = [key for key in npms_dict.keys()]
    for i in np.arange(len(keys)):
        for j in np.arange(i, len(keys)):
             if sorted(keys[i])  ==  sorted(keys[j]):
                if keys[i] is not keys[j]:
                    print(keys[i], keys[j])
                    npms_dict[keys[i]] +=  npms_dict[keys[j]]
                    npms_dict.pop(keys[j])
    dict['sum_tS_DICT_npms'] = npms_dict
    return dict  
#********************************************************************************************
def plot_dict(*argv):
    '''input: x, y, legend, title, xlims '''
    x = argv[0]
    y = argv[1]
    legend = argv[2]
    title = argv[3]
    xlims = argv[4]
    plt.figure(figsize = [8, 6])
    for key in legend:
        plt.plot(x, y[key])
    plt.legend(legend)
    plt.xlim([xlims[0], xlims[1]])
    plt.title(title)
#********************************************************************************************
def plot_list(*argv):
    '''input: x, y, legend, title, xlims '''
    x = argv[0]
    y = argv[1]
    legend = argv[2]
    title = argv[3]
    xlims = argv[4]
    plt.figure(figsize = [8, 6])
    for key in legend:
        plt.plot(x, y)
    plt.legend(legend)
    plt.xlim([xlims[0], xlims[1]])
    plt.title(title)
#********************************************************************************************
# def get_lims_gui(pts):
#     master = ttk.Tk()
# #     def get_entry_fields():
# #         print("xlim1: %s\nxlim2: %s" % (e1.get(), e2.get()))    
#     ttk.Label(master, text = "xlim1 min {:3.1f} or -1".format(pts[0])).grid(row = 0)
#     ttk.Label(master, text = "xlim2 max {:3.1f} or -1".format(pts[-1])).grid(row = 1)
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
# def set_xlim(pts1, pts3, lim1 = -2, lim2 = -2):
#     '''Setting plot xlimits, L < 0 : domain xlimits'''
#     if lim1 :w
#   -2 and lim2 !=  -2:
#         xlim1, xlim2 = lim1, lim2
#     else:
#         if   max(pts1[-1], pts3[-1]) > 700: xlim1, xlim2 = 350, 450
#         elif max(pts1[-1], pts3[-1]) <=  400 and max(pts1[-1], pts3[-1]) >=  200 : xlim1, xlim2 = 150, 250
#         elif max(pts1[-1], pts3[-1]) >=  250 and max(pts1[-1], pts3[-1]) <=  360 : xlim1, xlim2 = 280, max(pts1[-1], pts3[-1])
#         elif max(pts1[-1], pts3[-1]) <=  150 and max(pts1[-1], pts3[-1]) >=  100 : xlim1, xlim2 = 80, max(pts1[-1], pts3[-1])
#         elif max(pts1[-1], pts3[-1]) <=  100 : xlim1, xlim2 = 40, max(pts1[-1], pts3[-1])
#         else : xlim1, xlim2 = -1, -1
#     return xlim1, xlim2
# x1, x2 = get_lims_gui(tS_VLUs[3])
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
                  , fsize = 20, lsize = 20, tsize = 20, bins = 6, tksize = 20, stsize = 20, lgsize = 20, wlsize = 3):
    new_names = []
    for name in names:
        if '_' in name:
            name = name.replace('_', '\_')
        if '#' in name:
            name = name.replace('#', '\#')
        name = '$'+name+'$'
        new_names.append(name)
    if lim1 < 0: lim1 = L[0]
    if lim2 < 0: lim2 = L[-1]
    plt.figure(figsize = [15, 10])
    for i, x in enumerate(names):
        #ysmoothed = gaussian_filter1d(xs[:, i], sigma = 0.1)
        plt.plot(L, xs[:, i], linewidth = wlsize)
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
               fsize = 20, lsize = 20, tsize = 20, bins = 6, tksize = 20, stsize = 20, lgsize = 20):  
    if lim1 < 0: lim1 = min(L1[0], L2[0])
    if lim2 < 0: lim2 = max(L1[-1], L2[-1])
    if len(names)//2 > 1:
        if len(names)%2  ==  1:
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
            ie +=  1
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
            ie +=  1   
        fig1.suptitle(title, fontsize = stsize)
        fig1.tight_layout()
        plt.savefig(title+'.png')
#********************************************************************************************
def plot_tc_CQT_results_phf_all_phases(npm, xph, t, TIM, L, lim1 = -10, lim2 = -10
                                             , fsize = 20, lsize = 20, tsize = 20, bins = 6, tksize = 20, stsize = 20, lgsize = 20, lwsize = 3, filename = 'NPM-CQT-postprocessed_'):
    if lim1 <=  0: lim1 = L[0]
    if lim2 <=  0: lim2 = L[-1]
    plt.figure(figsize = [15, 15])
    for ph in npm.keys():
        if len(npm[ph]) !=  0:
            #ysmoothed = gaussian_filter1d(npm[ph], sigma = 0.1)
            plt.plot(xph[ph][0], npm[ph], linewidth = lwsize)
    plt.legend(npm.keys(), fontsize = lgsize)
    #plt.title(' Corrected NPM t = {:4.0f}sec, Post processed with all phases in the database'.format(TIM[t]), fontsize = 20)
    plt.xlabel(r'Distance ($\mu m$)', fontsize = lsize)
    plt.ylabel('Phase fraction', fontsize = lsize)
    plt.xticks(fontsize = tksize)
    plt.yticks(fontsize = tksize)
    plt.xlim([lim1, lim2])
    plt.savefig(filename+'-{:4.0f}-sec.png'.format(TIM[t]))
#********************************************************************************************
def plot_zero_padded_npm_CQT(Lx2, npm_pp_CQT, ttx2, lim1, lim2, L
                                   , fsize = 20, lsize = 20, tsize = 20, bins = 6, tksize = 20, stsize = 20, lgsize = 20, lwsize = 3, filename = 'NPM-zeroPadded-Corrected'):
    if lim1 <=  0: lim1 = L[0]
    if lim2 <=  0: lim2 = L[-1]
    r = np.zeros(len(Lx2)+1)
    for i in range(1, len(Lx2)+1):
        r[i] = (Lx2[i-1]-r[i-1])+Lx2[i-1]
    A = []
    for ph in npm_pp_CQT.keys():
        A.append(ph)
    plt.figure(figsize = [15, 10])
    for ph in A:
        vertices_npm = np.zeros(len(r))
        vertices_npm[0] = npm_pp_CQT[ph][0]
        for i in range(1, len(vertices_npm)-1):
                vertices_npm[i] = (npm_pp_CQT[ph][i]+npm_pp_CQT[ph][i-1])/2
        vertices_npm[-1] = npm_pp_CQT[ph][-1]
        ysmoothed = vertices_npm#gaussian_filter1d(vertices_npm, sigma = 0.1)
        plt.plot(r, ysmoothed, linewidth = lwsize)
    plt.legend(A, fontsize = lgsize)
    plt.xlabel(r'Distance ($\mu m$)', fontsize = lsize)
    plt.ylabel('Phase fraction', fontsize = lsize)
    plt.xticks(fontsize = tksize)
    plt.yticks(fontsize = tksize)
    plt.xlim([lim1, lim2])
    #plt.title('zero padded CQT Phase fraction t = {:2.1e} sec'.format(ttx2), fontsize = 20)
    plt.savefig(filename+'-{:4.0f}-sec.png'.format(ttx2))
#********************************************************************************************
def plot_x_3ts_log10(x1, x2, x3, L1, L2, L3, tt1, tt2, tt3, names, lim1 = -10, lim2 = -10, title = '', 
               file_name = '3ts', ylab = "fractions",  fsize = 20, lsize = 20, tsize = 20, bins = 6, 
               tksize = 20, stsize = 20, lgsize = 20, lwsize = 3):
    params = {'text.usetex': True}
    plt.rcParams.update(params)
    
    if lim1 < 0: lim1 = min(L1[0], L2[0])
    if lim2 < 0: lim2 = max(L1[-1], L2[-1])
    if len(names)//2 > 1:
        if len(names)%2  ==  1:
            fig1, ax1 = plt.subplots(len(names)//2+1, 2, figsize = [15, (len(names)//2+1)*25/5])
        else:
            fig1, ax1 = plt.subplots(len(names)//2, 2, figsize = [15, len(names)//2*25/5])
        ie = 0
        for el in names:
            ax1[ie//2, ie%2].plot(L1, -1*np.log10(np.abs(x1[:, ie])), label = 't = {:4.0f} sec'.format(tt1) , linewidth = lwsize)
            ax1[ie//2, ie%2].plot(L2, -1*np.log10(np.abs(x2[:, ie])), label = 't = {:4.0f} sec'.format(tt2) , linewidth = lwsize)
            ax1[ie//2, ie%2].plot(L3, -1*np.log10(np.abs(x3[:, ie])), label = 't = {:4.0f} sec'.format(tt3) , linewidth = lwsize)
            ax1[ie//2, ie%2].legend(fontsize = lgsize)
            ax1[ie//2, ie%2].set_title(el, fontsize = tsize)
            ax1[ie//2, ie%2].set_xlabel(r'Distance ($\mu m$)', fontsize = lsize)
            ax1[ie//2, ie%2].set_ylabel(ylab, fontsize = lsize)
            ax1[ie//2, ie%2].tick_params(axis = "x", labelsize = tksize)
            ax1[ie//2, ie%2].tick_params(axis = "y", labelsize = tksize)
            ax1[ie//2, ie%2].set_xlim([lim1, lim2])
            ax1[ie//2, ie%2].locator_params(axis = 'y', nbins = bins)
            ax1[ie//2, ie%2].locator_params(axis = 'x', nbins = bins)
            ie +=  1
        fig1.suptitle(title, fontsize = stsize)
        fig1.tight_layout()
        plt.savefig(file_name+'.png')
    else:
        fig1, ax1 = plt.subplots(len(names), 1, figsize = [15, len(names)*25/5])
        ie = 0
        for el in names:
            ax1[ie].plot(L1, -1*np.log10(np.abs(x1[:, ie])), label = 't = {:4.0f} sec'.format(tt1), linewidth = lwsize)
            ax1[ie].plot(L2, -1*np.log10(np.abs(x2[:, ie])), label = 't = {:4.0f} sec'.format(tt2), linewidth = lwsize)
            ax1[ie].plot(L3, -1*np.log10(np.abs(x3[:, ie])), label = 't = {:4.0f} sec'.format(tt3), linewidth = lwsize)
            ax1[ie].legend(fontsize = lgsize)
            ax1[ie].set_title(el, fontsize = tsize)
            ax1[ie].set_xlabel(r'Distance ($\mu m$) ', fontsize = lsize)
            ax1[ie].set_ylabel(' ', fontsize = lsize)
            ax1[ie].tick_params(axis = "x", labelsize = tksize) 
            ax1[ie].tick_params(axis = "y", labelsize = tksize)
            ax1[ie].locator_params(axis = 'y', nbins = bins)
            ax1[ie].locator_params(axis = 'x', nbins = bins)       
            ax1[ie].set_xlim([lim1, lim2])
            ie +=  1   
        fig1.suptitle(title, fontsize = stsize)
        fig1.tight_layout()
        plt.savefig(file_name+'.png')       
#********************************************************************************************      
def all_x_plotter_log10(names, L, tt, xs, lim1 = -10, lim2 = -10, title = '', file_name = 'all_values', ylab = "fractions"\
                  , fsize = 20, lsize = 20, tsize = 20, bins = 6, tksize = 20, stsize = 20, lgsize = 20, wlsize = 3):
    
    if lim1 < 0: lim1 = L[0]
    if lim2 < 0: lim2 = L[-1]
    plt.figure(figsize = [15, 10])
    for i, x in enumerate(names):
        #ysmoothed = gaussian_filter1d(xs[:, i], sigma = 0.1)
        plt.plot(L, np.log10(xs[:, i]), linewidth = wlsize)
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
                                , fsize = 20, lsize = 20, tsize = 20, bins = 6, tksize = 20, stsize = 20, lgsize = 20):
    '''grid pt = data[n]
    phase names = data[n].keys() 
    phase fraction = data[n][keys][0]
    el fractions in phase = data[n][keys][1]
    ordered el names in phase = data[n][keys][2]
    ordered el fractions in phase = data[n][keys][3]'''
    params = {'text.usetex': False}
    plt.rcParams.update(params)
    if lim1 <=  0: lim1 = L[0]
    if lim2 <=  0: lim2 = L[-1]
    if len(ph_names)//2 > 1:
        if len(ph_names)%2  ==  1:
            fig1, ax1 = plt.subplots(len(ph_names)//2+1, 2, figsize = [15, (len(ph_names)//2+1)*25/5])
        else:
            fig1, ax1 = plt.subplots(len(ph_names)//2, 2, figsize = [15, len(ph_names)//2*25/5])
        iph = 0
        for iph, ph in enumerate(ph_names):  
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
            iph +=  1
        fig1.tight_layout()
        plt.savefig('XinPh-postprocessed-{:.0f}-sec.png'.format(t))  
    elif len(new_ph_names)  ==  1:
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
        for iph, ph in enumerate(ph_names):
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
            iph +=  1
        fig1.tight_layout()
        plt.savefig('XinPh-postprocessed-{:4.0f}-sec.png'.format(t))        
#********************************************************************************************        
def plot_allx_2ts(x1, x2, L1, L2, tt1, tt2, names, lim1 = -10, lim2 = -10, title = '', file_name = '', fsize = 20, lsize = 20, tsize = 20, bins = 6, tksize = 20, stsize = 20, lgsize = 20, lwsize = 20):
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
        ax.plot(L1, x1[:, ie], linewidth = lwsize, label = '{}, {:.0f} sec'.format(el, tt2) )
        ax.plot(L2, x2[:, ie], '--', linewidth = lwsize, label = '{}, {:.0f} sec'.format(el, tt1) )
        ax.set_xlabel('Distance [um]', fontsize = lsize)
        ax.set_ylabel('Mole Fractions', fontsize = lsize)
        ax.tick_params(axis = "x", labelsize = tksize)
        ax.tick_params(axis = "y", labelsize = tksize)
        ax.set_xlim([lim1, lim2])
        ax.locator_params(axis = 'y', nbins = bins)
        ax.locator_params(axis = 'x', nbins = bins)
        ie +=  1
    ax.legend(fontsize = lgsize)
    plt.savefig(file_name+'.png')
