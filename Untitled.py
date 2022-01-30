#!/usr/bin/env python
# coding: utf-8
#%%
from clibv2 import * #%matplotlib qt#user_path = !eval echo ~$USER
import matplotlib.pyplot as plt

'''Pre-settings'''
phase_changes = [('BCC_A2', 'W', 1), ('FCC_A1', 'N', 0.03)] #phase_change = [('BCC_A2', 'W', 1), ('FCC_A1', 'C', 0.03)]#, ('FCC_A1', 'Ti', 0.05), ('FCC_A1', 'Al', 0.05)]
tc_setting = {'database':'tcfe9', 'acRefs':[], 'phsToSus':[], 'p3flag':True} 
name_pairs = [('BCC_A2#1', 'Ti64 BCC'), ('FCC_A1#1', 'Co FCC'), ('MC_SHP#1', 'WC'), ('M6C#1', 'Eta M6'), 
                ('BCC_A2-W', 'W BCC'), ('FCC_A1-CTI', 'TiC'), ('FCC_A1-CAL', 'AlC'), ('FCC_A1-ALC', 'AlC'), 
                ('LIQUID#1', 'LIQUID'), ('FCC_A1-TIC', 'TiC')]

'''get values from TXT files'''
try:
    #VLUs = get_values_from_textfiles() 
    VLUs = get_values_from_textfiles('/Users/salmasi/LUND/REALPCBN/test/')
    VLUs['phase_changes'] = phase_changes
    VLUs['tc_setting'] = tc_setting 
    VLUs['name_pairs'] = name_pairs 
except:
    print('Error reading text files')
else:
    print('Text files are loaded in memory') 

'''Select time steps'''
try:
    tS, nearestTime = get_timeStamp(VLUs['times'])
    print(tS, nearestTime)
except:
    print('Selected simulation time is not in range')
else:
    print('Last selected time is in range, timeStamp is {}'.format(tS))

'''get values of timesteps'''
tS_VLUs = get_tS_VLUs(VLUs, tS, nearestTime)

''' get plot limits and TC pre-settings'''
tS_VLUs['xlim1'] = float(input('xlims = /{}, {}/x1?/:\n'.format(tS_VLUs['tS_pts'][0], tS_VLUs['tS_pts'][-1])))
tS_VLUs['xlim2'] = float(input('xlims = /{}, {}/x2?/:\n'.format(tS_VLUs['tS_pts'][0], tS_VLUs['tS_pts'][-1])))
### todo : automatic alternating

''' TC calculations of tS_VLUs'''
tS_tc_VLUs = tccalc(tS_VLUs)

#%%
''' Plotting'''
plot_list(tS_tc_VLUs['tS_pts'], tS_tc_VLUs['tS_DICT_mfs'], tS_tc_VLUs['elnames'], 
            'MFs , t = {}'.format(tS_tc_VLUs['nearestTime']),
            (tS_tc_VLUs['xlim1'], tS_tc_VLUs['xlim2']))
plot_list(tS_tc_VLUs['tS_pts'], tS_tc_VLUs['tS_DICT_ufs'], tS_tc_VLUs['elnames'], 
            'TC calculated values, UFs , t = {}'.format(tS_tc_VLUs['nearestTime']),
            (tS_tc_VLUs['xlim1'], tS_tc_VLUs['xlim2']))
plot_dict(tS_tc_VLUs['tS_pts'], tS_tc_VLUs['tS_TC_ws'], tS_tc_VLUs['tS_TC_ws'].keys(), 
            'TC calculated values, WFs , t = {}'.format(tS_tc_VLUs['nearestTime']), 
            (tS_tc_VLUs['xlim1'], tS_tc_VLUs['xlim2']))
plot_list(tS_tc_VLUs['tS_pts'], tS_tc_VLUs['tS_DICT_npms'], 
            tS_tc_VLUs['tS_DICT_phnames'], 'DICT TXT FILES, t = {}'.format(tS_tc_VLUs['nearestTime']), 
            (tS_tc_VLUs['xlim1'], tS_tc_VLUs['xlim2']))
plot_dict(tS_tc_VLUs['tS_pts'], tS_tc_VLUs['sum_tS_DICT_npms'] ,tS_tc_VLUs['sum_tS_DICT_npms'].keys(),
            'DICT TXT FILES, Composition sets are summed, t = {}'.format(tS_tc_VLUs['nearestTime']), 
            (tS_tc_VLUs['xlim1'], tS_tc_VLUs['xlim2'])) 
plot_dict(tS_tc_VLUs['tS_pts'], tS_tc_VLUs['tS_TC_NEAT_npms'], tS_tc_VLUs['tS_TC_NEAT_npms'].keys(), 
            'TC calculated values, NPM, t = {}'.format(tS_tc_VLUs['nearestTime']), 
            (tS_tc_VLUs['xlim1'], tS_tc_VLUs['xlim2']))
plot_dict(tS_tc_VLUs['tS_pts'], tS_tc_VLUs['tS_TC_NEAT_vpvs'], tS_tc_VLUs['tS_TC_NEAT_vpvs'].keys(), 
            'TC calculated values, VPV, t = {}'.format(tS_tc_VLUs['nearestTime']), 
            (tS_tc_VLUs['xlim1'], tS_tc_VLUs['xlim2']))
plot_dict(tS_tc_VLUs['tS_pts'], tS_tc_VLUs['nameChanged_CQT_tS_TC_NEAT_npms'], tS_tc_VLUs['nameChanged_CQT_tS_TC_NEAT_npms'].keys(),
            'TC calculated values, summed NPMs changed names, t = {}'.format(tS_tc_VLUs['nearestTime']),
            (tS_tc_VLUs['xlim1'], tS_tc_VLUs['xlim2']))
plot_dict(tS_tc_VLUs['tS_pts'], tS_tc_VLUs['sum_CQT_tS_TC_NEAT_npms'], tS_tc_VLUs['sum_CQT_tS_TC_NEAT_npms'].keys(),
            'TC calculated values, summed corrected neat NPMs, t = {}'.format(tS_tc_VLUs['nearestTime']),
            (tS_tc_VLUs['xlim1'], tS_tc_VLUs['xlim2']))
plot_dict(tS_tc_VLUs['tS_pts'], tS_tc_VLUs['CQT_tS_TC_NEAT_npms'], tS_tc_VLUs['CQT_tS_TC_NEAT_npms'].keys(),
            'TC calculated values, corrected neat NPMs, t = {}'.format(tS_tc_VLUs['nearestTime']),
            (tS_tc_VLUs['xlim1'], tS_tc_VLUs['xlim2']))
plt.show()

# 'tS_pts'# 'xlim1', # 'xlim2', 
# 'tS_DICT_npms'
# 'tS_DICT_phnames'
# 'tS_DICT_mfs', 
# 'tS_DICT_mus'
# 'tS_DICT_ufs'
# 'sum_tS_DICT_npms'

# 'tS_TC_phnames', 
# 'tS_TC_npms', 
# 'tS_TC_vpvs', 
# 'tS_TC_phXs', 
# 'tS_TC_acRef', 
# 'tS_TC_acSER', 
# 'tS_TC_mus', 
# 'tS_TC_ws', 

# 'tS_TC_NEAT_phnames', 
# 'tS_TC_NEAT_npms', 
# 'tS_TC_NEAT_vpvs', 
# 'tS_TC_NEAT_phXs', 
# 'tS_TC_NEAT_mfs', 

# 'CQT_tS_TC_NEAT_phXs', 
# 'CQT_tS_TC_NEAT_npms', 
# 'sum_CQT_tS_TC_NEAT_npms', 

# 'nameChanged_CQT_tS_TC_NEAT_npms', 
# %%
