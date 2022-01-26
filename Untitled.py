#!/usr/bin/env python
# coding: utf-8
#todo: one dict data structure or object 
from clibv2 import *  #%matplotlib qt#user_path = !eval echo ~$USER
import matplotlib.pyplot as plt
'''reading values from TXT files'''
try:
    vals = get_values_from_textfiles() #vals = get_values_from_textfiles('/Users/salmasi/LUND/WC15Co-Ti64-200um200um/1000C2/')
except:
    print('Error reading text files')
else:
    print('Text files are loaded in memory')  
'''Selecting time steps'''
try:
    vals['tStamp'] = []
    vals['nearestTime'] = []
    tmp = get_timeStamp(vals['times'])
    vals['tStamp'].append(tmp[0])
    vals['nearestTime'].append(tmp[1])
except:
    print('Selected simulation time is not in range')
else:
    print('Last selected time is in range, timeStamp is {}'.format(vals['tStamp'][-1]))
'''get values of selected time steps from Data'''

tStamp_vals = get_tStamp_vals(vals)
# todo: add uf_calc to get_values_of_timeStamp
'''plot limits'''
tStamp_vals['xlim1'] = float(input('xlims =/{},{}/x1?/:'.format(tStamp_vals['tStamp_points'][0],tStamp_vals['tStamp_points'][-1])))
tStamp_vals['xlim2'] = float(input('xlims =/{},{}/x2?/:'.format(tStamp_vals['tStamp_points'][0],tStamp_vals['tStamp_points'][-1])))
'''Pre settings'''
tStamp['phase_change'] = [('BCC_A2','W',1),('FCC_A1','N',0.03)] #phase_change = [('BCC_A2','W',1),('FCC_A1','C',0.03)]#,('FCC_A1','Ti',0.05),('FCC_A1','Al',0.05)]
tStamp['tc_Setting']=['tcfe9',[],[]] #database,activityRefs,phsToSuspend# phsToSuspend =['gas']#phsToSuspend =['graphite']
tStamp.['name_paires'] = [('BCC_A2#1','Ti64 BCC'),('FCC_A1#1','Co FCC'), ('MC_SHP#1','WC'),('M6C#1','Eta M6'),('BCC_A2-W','W BCC'),('FCC_A1-CTI','TiC'),('FCC_A1-CAL','AlC'),('FCC_A1-ALC','AlC'),('LIQUID#1','LIQUID'),('FCC_A1-TIC','TiC')]
# todo : automatic alternating

tStamp_vals_summed = add_phs_compositionDets_DICTRA(tStamp_vals)

'''Thermodynaic calculations of the timestamp'''
T,subidx,elnames,pth,poly3_file = vals[9],vals[10],vals[2],vals[0],False
tStamp_tcProc_vals = tccalc(tStamp_vals, database, phsToSuspend, activityRefs, T, elnames, pth, poly3_file)
# todo : calculate u-fractions
trimed_tcCalculated_vals = trim_tcCalculated_values(tStamp_vals,tStamp_tcProc_vals,elnames,subidx)

'''Correct phase indices to account for carbides and nitrides or any change in composition sets'''
_,corrected_phase_indices_vals=correct_phase_indices(trimed_tcCalculated_vals,phase_change)
summed_npms_dict = add_composition_sets(corrected_phase_indices_vals)
name_changed_npms_dict = ph_namechange(summed_npms_dict,name_paires)

plot_list(tStamp_vals[3],tStamp_vals[1],tStamp_vals[4],'DICTRA TXT FILES, t={}'.format(nearestTime),(x1,x2))
plot_dict(tStamp_vals[3], tStamp_vals_summed, tStamp_vals_summed.keys(), 'DICTRA TXT FILES, Composition sets are summed, t={}'.format(nearestTime),(x1,x2))    
plot_dict(trimed_tcCalculated_vals[0], trimed_tcCalculated_vals[2], trimed_tcCalculated_vals[2].keys(), 'TC calculated values, NPM, t={}'.format(nearestTime), (x1,x2))
plot_dict(trimed_tcCalculated_vals[0], trimed_tcCalculated_vals[3], trimed_tcCalculated_vals[3].keys(), 'TC calculated values, VPV, t={}'.format(nearestTime), (x1,x2))
plot_dict(trimed_tcCalculated_vals[0], summed_npms_dict, summed_npms_dict.keys(), 'TC calculated values,summed NPMs, t={}'.format(nearestTime), (x1,x2))
plot_dict(trimed_tcCalculated_vals[0], name_changed_npms_dict, name_changed_npms_dict.keys(), 'TC calculated values,summed NPMs changed names, t={}'.format(nearestTime), (x1,x2))
plot_dict(trimed_tcCalculated_vals[0], trimed_tcCalculated_vals[1], trimed_tcCalculated_vals[11], 'TC calculated values, MFs , t={}'.format(nearestTime), (x1,x2))
plot_dict(trimed_tcCalculated_vals[0], uf_dict, uf_dict.keys(), 'TC calculated values, UFs , t={}'.format(nearestTime), (x1,x2))

plt.show()