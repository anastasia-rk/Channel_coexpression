# In this file, we generating conductance multipliers independently from log-normal distributions with sample statistics
# from the mass action kinetics data
import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
# import pandas as pd
import myokit as mk
import os
import csv
import pickle
import gc
from tqdm import tqdm
matplotlib.use('Agg')
from mass_action_correct import *

# definitions
def run_simulation_w_points(expressions_to_vary,multiplier_pairs,simModel,DataFolderName,FigureFolderName,folderName,nTries=1000):
    print('Expressions to vary: ', expressions_to_vary)
    print('Using pre-generated conductance multipliers')
    # expression gains reset to match the baseline model operation
    expression_gains = [1, 1, 1, 1, 1]
    # simgas of lognormal distributions from which the gains will be sampled
    all_sigmas = ['NA', 'NA', 'NA', 'NA', 'NA']  # this is only for storing
    # get the indices of the currents to vary
    indices = [currentNames.index(expression_to_vary) for expression_to_vary in expressions_to_vary]
    # create a dictionary to store the simulation results
    allKeys = gainNames + ['Time'] + currentNames + ['Voltage']
    simulationResults = dict.fromkeys(allKeys)
    # this dictionary will contain the lists of values for each simulation index
    for key in simulationResults.keys():
        simulationResults[key] = []
    # counters
    countAlternans = 0
    # initialise biomarker dictionary
    biomarkers = dict.fromkeys(
        ['APD90', 'APD60', 'APD40', 'APA', 'TRI60', 'TRI40', 'EAD', '90PERCENT', '90to90', 'Alternan'])
    # initialise the biomarker dictionary with empty lists
    for key in biomarkers.keys():
        biomarkers[key] = []
    # figure for plotting currents
    nThingsToPlot = 1 + len(currentNamesToPlot)
    nRows = int(np.ceil(nThingsToPlot / 3))
    height = 5
    fig, axs = plt.subplots(nRows, 3, figsize=(height*3, height*nRows))
    axs = axs.flatten()
    # loop over the number of simulations
    for iTry in tqdm(range(nTries)):
        # read the gains on conductance from the list
        gain_values = multiplier_pairs[iTry]
        for iGain, gain_value in enumerate(gain_values):
            expression_gains[indices[iGain]] = gain_value
        # get the list of conductances with the new expression gains
        conductances_sampled = [baselineConductances[i] * expression_gains[i] for i in range(len(baselineConductances))]
        # set constants in the mmt file to new values
        for iConductance, conductance in enumerate(conductances_sampled):
            simModel.set_constant(conductanceNamesInMMT[iConductance], conductance)
        # run the simulation
        d = simModel.run(t_end, log=vars_to_log)
        # get the times as array
        times = np.array(d.time())
        # find the index of d.time() from which the last 4 paces start
        last_paces_start = np.where(times > t_start_storing)[0][
                               0] - 1  # include one time point before because we dont have regular time increments
        # print('time at which storing starts: ', times[last_paces_start])
        # keep only the last 4 paces
        times = times[last_paces_start:]
        V_last_four_paces = d['membrane.v'][last_paces_start:]
        # plot the simulation outputs
        axs[0].plot(times / 1000, V_last_four_paces)  # plot the voltage separately
        for iCurrent, currentName in enumerate(currentNamesInMMT):
            if currentName in currentNamesToPlot:
                # plot the current
                axs[currentNamesToPlot.index(currentName) + 1].plot(times / 1000, d[currentName][last_paces_start:])
            # store the current values in the dictionary
            simulationResults[currentNames[iCurrent]].append(d[currentName][last_paces_start:])
        del d
        gc.collect()
        # store the rest of the simulation results
        for iGain, gainName in enumerate(gainNames):
            simulationResults[gainName].append(expression_gains[iGain])
        simulationResults['Time'].append(times)
        simulationResults['Voltage'].append(V_last_four_paces)
        ################################################################################################################
        # analysing the AP traces
        # split the last four paces into two groups by finding indexes of elements in times
        APs = []
        AP_times = []
        AP_interpolators = []
        for i_beat in range(len(times_of_beats) - 1):
            indexes = \
                np.where((times >= t_end - times_of_beats[i_beat]) & (times < t_end - times_of_beats[i_beat + 1]))[0]
            np.where((times >= t_end - times_of_beats[i_beat]) & (times < t_end - times_of_beats[i_beat + 1]))[0]
            APs.append(V_last_four_paces[indexes[0]:indexes[-1]])
            AP_times.append(times[indexes[0]:indexes[-1]])
            bs_function = sp.interpolate.CubicSpline(AP_times[i_beat] - (t_end - times_of_beats[i_beat]), APs[i_beat])
            AP_interpolators.append(bs_function(np.arange(10, 900, 0.1)))  # interpolate at 1 ms intervals
        # check for alternans
        isAnAlternan = False
        for iTrace, Trace in enumerate(AP_interpolators):
            for jTrace in range(iTrace + 1, len(AP_interpolators)):
                # сheck if the difference between the traces is greater than 1 mV
                if np.any(np.abs(np.diff(Trace - AP_interpolators[jTrace])) >= 1):
                    isAnAlternan = True
                    countAlternans += 1
                    biomarkers['Alternan'].append(isAnAlternan)
                    # add a nan to all biomarkers
                    for key in biomarkers.keys():
                        if key != 'Alternan':
                            biomarkers[key].append(np.nan)
                    print('Alternan detected. Total number of alternans: ', countAlternans)
                    break
        if ~isAnAlternan:
            # use the last beat to calculate biomarkers
            iTrace = -1
            Trace = APs[iTrace]
            # check if early afterdepolarisation was triggered
            peaks, peak_properties = sp.signal.find_peaks(Trace, height=0)
            if len(peaks) > 1:
                time_between_peaks = AP_times[iTrace][peaks[-1]] - AP_times[iTrace][peaks[0]]
                if time_between_peaks >= 50:
                    biomarkers['EAD'].append(time_between_peaks)
            #  find maximum voltage and its index
            maxVoltage = np.max(Trace)
            maxVoltageIndex = np.where(Trace == maxVoltage)[0][0]
            # find the minimum voltage and its index after the peak
            minVoltageAfterPeak = np.min(Trace[maxVoltageIndex:])
            minVoltageIndex = np.where(Trace == minVoltageAfterPeak)[0][0]
            biomarkers['APA'].append(maxVoltage - minVoltageAfterPeak)
            # find the 90% repolarisation voltage
            repolarisationVoltage90 = 0.1 * (maxVoltage - minVoltageAfterPeak) + minVoltageAfterPeak
            repolarisationVoltage60 = 0.4 * (maxVoltage - minVoltageAfterPeak) + minVoltageAfterPeak
            repolarisationVoltage40 = 0.6 * (maxVoltage - minVoltageAfterPeak) + minVoltageAfterPeak
            # find the index at which the voltage crosses the 90% repolarisation voltage
            repolarisationIndex90 = np.where(Trace[maxVoltageIndex:] <= repolarisationVoltage90)[0][0] + maxVoltageIndex
            repolarisationIndex60 = np.where(Trace[maxVoltageIndex:] <= repolarisationVoltage60)[0][0] + maxVoltageIndex
            repolarisationIndex40 = np.where(Trace[maxVoltageIndex:] <= repolarisationVoltage40)[0][0] + maxVoltageIndex
            # find the index of the voltage below 90% repolarisation voltage before the max value
            beforepeakIndex90 = np.where(Trace[:maxVoltageIndex] >= repolarisationVoltage90)[0][0]
            # find APD90, APD60, APD40
            APD90 = AP_times[iTrace][repolarisationIndex90] - AP_times[iTrace][maxVoltageIndex]
            APD60 = AP_times[iTrace][repolarisationIndex60] - AP_times[iTrace][maxVoltageIndex]
            APD40 = AP_times[iTrace][repolarisationIndex40] - AP_times[iTrace][maxVoltageIndex]
            APD90to90 = AP_times[iTrace][repolarisationIndex90] - AP_times[iTrace][beforepeakIndex90]
            # compute trianglation slopes as average slope of the trace from APD40 and APD60 to APD90
            triangulationSlope40 = np.mean(np.diff(Trace[repolarisationIndex40:repolarisationIndex90]) / np.diff(
                AP_times[iTrace][repolarisationIndex40:repolarisationIndex90]))
            triangulationSlope60 = np.mean(np.diff(Trace[repolarisationIndex60:repolarisationIndex90]) / np.diff(
                AP_times[iTrace][repolarisationIndex60:repolarisationIndex90]))
            # store the biomarkers
            biomarkers['APD90'].append(APD90)
            biomarkers['APD60'].append(APD60)
            biomarkers['APD40'].append(APD40)
            biomarkers['90PERCENT'].append(repolarisationVoltage90)
            biomarkers['90to90'].append(APD90to90)
            biomarkers['TRI60'].append(triangulationSlope60)
            biomarkers['TRI40'].append(triangulationSlope40)
            biomarkers['Alternan'].append(isAnAlternan)
            # find peaks in the AP trace
        ################################################################################################################
        # reset the simulation
        simModel.reset()
    #  end for over simulations
    for ax in axs:
        ax.set_xlabel('time,s')
        ax.set_ylabel(ylabels[axs.tolist().index(ax)])
    # fig.suptitle(r'Covarying $g_{Kr}$ and $g_{CaL}$', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
    figName = FigureFolderName + '/' + folderName + '/'
    for ax in axs:
        ax.set_xlim((t_start_storing) / 1000, t_end / 1000)
    plt.savefig(figName + 'last_four_paces.png', dpi=300)
    for ax in axs:
        ax.set_xlim((t_end - 1000) / 1000, t_end / 1000)
    plt.savefig(figName + 'last_pace.png', dpi=300)
    # plt.show()
    del fig
    # create a figure and make subplots for each biomarker as scutter plots
    fig1, axs1 = plt.subplots(2, 4, figsize=(15, 10))
    axs1 = axs1.ravel()
    # plot the biomarkers
    iAxs = 0
    for iBiomarker, biomarker in enumerate(biomarkers.keys()):
        if biomarker == 'Alternan' or biomarker == 'EAD':
            continue
        axs1[iAxs].scatter(simulationResults['gain_kr'], biomarkers[biomarker], s=5, label=r'Multiplier on $g_{Kr}$')
        axs1[iAxs].scatter(simulationResults['gain_ca'], biomarkers[biomarker], s=5, label=r'Multiplier on $g_{CaL}$')
        axs1[iAxs].set_xlabel(r'Expression gain values')
        axs1[iAxs].set_ylabel(biomarker)
        axs1[iAxs].set_xscale('symlog', base=2, linthresh=0.125)
        axs1[iAxs].set_xticks([0.25, 0.5, 1, 2, 4, 8])
        if iAxs == 2:
            axs1[iAxs].legend()
        iAxs += 1
    fig1.suptitle(r'Biomarkers of AP traces obtained from gain covariation', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
    plt.savefig(figName + 'biomarkers.png', dpi=300)
    # plot histograms of the biomarkers
    fig2, axs2 = plt.subplots(2, 4, figsize=(15, 10))
    axs2 = axs2.ravel()
    ####################################################################################################################
    # plot the biomarkers
    iAxs = 0
    for iBiomarker, biomarker in enumerate(biomarkers.keys()):
        if biomarker == 'Alternan' or biomarker == 'EAD':
            continue
        axs2[iAxs].hist(biomarkers[biomarker], bins=20, density=True, histtype='step')
        # add mean, median and std to the plot
        axs2[iAxs].axvline(np.mean(biomarkers[biomarker]), color='k', linestyle='solid', linewidth=1)
        axs2[iAxs].axvline(np.median(biomarkers[biomarker]), color='k', linestyle='dashed', linewidth=1)
        axs2[iAxs].axvline(np.median(biomarkers[biomarker]) + np.std(biomarkers[biomarker]), color='grey',
                           linestyle='dotted', linewidth=1)
        axs2[iAxs].axvline(np.mean(biomarkers[biomarker]) - np.std(biomarkers[biomarker]), color='grey',
                           linestyle='dotted', linewidth=1)
        ymin, ymax = axs2[iAxs].get_ylim()
        axs2[iAxs].text(np.mean(biomarkers[biomarker]), ymax * 0.9,
                        'mean: {:.2f}'.format(np.mean(biomarkers[biomarker])))
        axs2[iAxs].text(np.median(biomarkers[biomarker]), ymax * 0.8,
                        'md: {:.2f}'.format(np.median(biomarkers[biomarker])))
        axs2[iAxs].text(np.median(biomarkers[biomarker]) + np.std(biomarkers[biomarker]), ymax * 0.7,
                        'std: {:.2f}'.format(np.std(biomarkers[biomarker])))
        axs2[iAxs].set_xlabel(biomarker)
        iAxs += 1
    fig2.suptitle(r'Histogramms of biomarkers obtained from gain covariation', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
    plt.savefig(figName + 'biomarker_hists.png', dpi=300)
    ########################################################################################################################
    # save the simulation results and biomarkers into pickle files
    if not os.path.exists(DataFolderName + '/' + folderName):
        os.makedirs(DataFolderName + '/' + folderName)
    pickle.dump(simulationResults, open(DataFolderName + '/' + folderName + '/simulationResults.pkl', 'wb'))
    pickle.dump(biomarkers, open(DataFolderName + '/' + folderName + '/biomarkers.pkl', 'wb'))
    del simulationResults, fig1, fig2, axs1, axs2, biomarkers
    gc.collect()
    return True

# main
if __name__ == '__main__':
    # if there is no folder for figures, create one
    FigureFolderName = 'Figures_Tomek_correct'
    DataFolderName = 'Simulated_data_Tomek_correct'
    if not os.path.exists(FigureFolderName):
        os.makedirs(FigureFolderName)
    # if there is no folder for output, create one
    if not os.path.exists(DataFolderName):
        os.makedirs(DataFolderName)
    # load the model and protocol
    model, protocol, embedded_script = mk.load('Tomek_fixed_cellml/Tomek_dynCl_endo.mmt')
    currentNamesInMMT = ['IKr.IKr', 'IKs.IKs', 'Ito.Ito', 'INa.INa', 'ICaL.ICaL']
    s = mk.Simulation(model, protocol)
    vars_to_log = ['environment.time', 'membrane.v'] + currentNamesInMMT
    d = s.run(1000,log=vars_to_log)
    filename = 'Tomek_cellml/test_output.csv'
    d.save_csv(filename, precision=64, order=None, delimiter=',', header=True)
    s.reset() # resets the state of the simulation

    # plot the results for the original and modified model
    currentNamesInMMT = ['IKr.IKr', 'IKs.IKs', 'Ito.Ito', 'INa.INa', 'ICaL.ICaL']
    ylabels = [r'$V$, mV', r'$I_{Kr}$, ', r'$I_{Ks}$, A/F', r'$I_{to}$, A/F', r'$I_{Na}$, A/F', r'$I_{CaL}$, A/F']
    fig, axs = plt.subplots(2, 3, figsize=(20, 10))
    axs = axs.flatten()
    axs[0].plot(d.time(), d['membrane.v'])  # plot the voltage separately
    for iCurrent, currentName in enumerate(currentNamesInMMT):
        # plot the current
        if iCurrent == 1:
            axs[iCurrent + 1].plot(d.time(), d[currentName])
        else:
            axs[iCurrent + 1].plot(d.time(), d[currentName])
    for ax in axs:
        ax.set_xlabel('time, s')
        ax.set_ylabel(ylabels[axs.tolist().index(ax)])
    fig.suptitle(r'Corrected Tomek model', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(FigureFolderName + '/Tomek_baseline.png', dpi=300)
    # plt.show()

    ## try altering conductances in the model
    model, protocol, embedded_script = mk.load('Tomek_fixed_cellml/Tomek_dynCl_endo_changing_conductances.mmt')
    ########################################################################################################################
    # set up the simulation - this is the same for all configurations
    s = mk.Simulation(model, protocol)
    ## these variables are shared by all other configurations and should not be reset
    nTries = 1000  # number of simulations
    # OHR model takes about 450 paces to equilibrate - from Mann et al. 2016
    t_end = 454000  # end time of the simulation - there could be a smarter way to define this in terms of paces
    times_of_beats = [4000, 3000, 2000, 1000, 0]  # these are in descending order because they will be substraced from t_end
    # for storage only: record the last 4 paces
    t_start_storing = t_end - 4000
    # baseline conductances of the currents of interest
    currentNames = ['IKr', 'IKs', 'Ito', 'INa', 'ICaL']
    currentNamesInMMT = ['IKr.IKr', 'IKs.IKs', 'Ito.Ito', 'INa.INa', 'ICaL.ICaL']
    currentNamesToPlot = ['IKr.IKr', 'ICaL.ICaL']
    conductanceNamesInMMT = ['IKr.GKr_b', 'IKs.GKs_b', 'Ito.Gto_b', 'INa.GNa', 'ICaL.PCa_b']
    gainNames = ['gain_kr', 'gain_ks', 'gain_to', 'gain_na', 'gain_ca']
    baselineConductances = [0.0321, 0.0011, 0.16, 11.7802, 8.37570000000000046e-5]
    print('baseline conductances before modification: ', baselineConductances)
    ylabels = [r'$V$, mV', r'$I_{Kr}$, ', r'$I_{Ks}$, A/F', r'$I_{to}$, A/F', r'$I_{Na}$, A/F', r'$I_{CaL}$, A/F']
########################################################################################################################
        # independently varying expression of I_Kr and I_CaL
########################################################################################################################
    print('Independently varied conductances for I_Kr and I_CaL')
    # create folder to store selected simulation outputs
    # expression gains reset to match the baseline model operation
    expression_gains = [1, 1, 1, 1, 1]
    # simgas of lognormal distributions from which the gains will be sampled
    all_sigmas = ['NA', 'NA', 'NA', 'NA', 'NA']  # this is only for storing
    # baseline conductances of the currents of interest
    expressions_to_vary = ['IKr', 'ICaL']
    folderName = expressions_to_vary[0] + '_' + expressions_to_vary[1] + '_independent_test'
    if not os.path.exists(DataFolderName + '/' + folderName):
        os.makedirs(DataFolderName + '/' + folderName)
    if not os.path.exists(FigureFolderName + '/' + folderName):
        os.makedirs(FigureFolderName + '/' + folderName)
    covars = pickle.load(open('Pickles/gtex_covariances.pkl', 'rb'))
    Sigma_bivariate = covars['CACNA1C']
    # we only want the diagonal element from the
    sigmas = [np.sqrt(Sigma_bivariate[0, 0]), np.sqrt(Sigma_bivariate[1, 1])]
    print('Expressions to vary: ', expressions_to_vary)
    print('Sigmas: ', sigmas)
    # using the sigmas generate a list of pairs of lognormally distributed values
    gain_values = np.random.lognormal(mean=0.0, sigma=sigmas, size=(nTries, len(sigmas)))
    SimulationSuccess = run_simulation_w_points(expressions_to_vary, gain_values, s, DataFolderName, FigureFolderName,
                                                folderName, nTries)
########################################################################################################################
    # covaried coexpression of I_Kr and I_CaL
########################################################################################################################
    print('Covaried conductances for I_Kr and I_CaL')
    # create folder to store selected simulation outputs
    # baseline conductances of the currents of interest
    expressions_to_vary = ['IKr', 'ICaL']
    folderName = expressions_to_vary[0] + '_' + expressions_to_vary[1] + '_cotranscripted'
    if not os.path.exists(DataFolderName + '/' + folderName):
        os.makedirs(DataFolderName + '/' + folderName)
    if not os.path.exists(FigureFolderName + '/' + folderName):
        os.makedirs(FigureFolderName + '/' + folderName)
    # load covariance from the pickle file created by the script 'read_gtex.py'
    covars = pickle.load(open('Pickles/gtex_covariances.pkl', 'rb'))
    Sigma_bivariate = covars['CACNA1C']  # np.array([[0.25, 0.25*(0.5/0.7)], [0.25*(0.5/0.7), 0.25]])
    sampled_from_normal = np.random.multivariate_normal(mean=[0, 0], cov=Sigma_bivariate, size=nTries)
    transcript_levels = np.exp(sampled_from_normal)
    SimulationSuccess = run_simulation_w_points(expressions_to_vary,transcript_levels, s, DataFolderName, FigureFolderName,
                                                folderName, nTries)
    ########################################################################################################################
    # covaried coexpression of I_Kr and I_CaL with cotransaltion covariance
    ########################################################################################################################
    print('Covaried conductances for I_Kr and I_CaL (with cotransaltion covariance)')
    # create folder to store selected simulation outputs
    # baseline conductances of the currents of interest
    expressions_to_vary = ['IKr', 'ICaL']
    folderName = expressions_to_vary[0] + '_' + expressions_to_vary[1] + '_cotranslated'
    if not os.path.exists(DataFolderName + '/' + folderName):
        os.makedirs(DataFolderName + '/' + folderName)
    if not os.path.exists(FigureFolderName + '/' + folderName):
        os.makedirs(FigureFolderName + '/' + folderName)
    # Sigma_bivariate = np.array([[0.25, 0.25*(0.5/0.7)], [0.25*(0.5/0.7), 0.25]])
    # load covariance from the pickle file created by the script 'read_gtex.py'
    Sigma_bivariate = pickle.load(open('Pickles/hERG_CACNA1C_cotranslational_covariance.pkl', 'rb'))
    # sample gains from bivarite lognormal distribution
    sampled_from_normal = np.random.multivariate_normal(mean=[0, 0], cov=Sigma_bivariate, size=nTries)
    protein_levels = np.exp(sampled_from_normal)
    SimulationSuccess = run_simulation_w_points(expressions_to_vary,protein_levels,s,DataFolderName,FigureFolderName,folderName,nTries)
########################################################################################################################
    # Dependent variation of ion channel conductances
########################################################################################################################
    print('Dependent variation of conductances for I_Kr and I_CaL (with cotransaltion covariance)')
    # create folder to store selected simulation outputs
    # baseline conductances of the currents of interest
    expressions_to_vary = ['IKr', 'ICaL']
    folderName = expressions_to_vary[0] + '_' + expressions_to_vary[1] + '_dependent'
    if not os.path.exists(DataFolderName + '/' + folderName):
        os.makedirs(DataFolderName + '/' + folderName)
    if not os.path.exists(FigureFolderName + '/' + folderName):
        os.makedirs(FigureFolderName + '/' + folderName)
    # sample gains from bivarite lognormal distribution
    gain_values = np.random.lognormal(mean=0.0, sigma=0.25, size=(nTries, 1))
    gain_values = np.concatenate((gain_values, gain_values), axis=1)
    SimulationSuccess = run_simulation_w_points(expressions_to_vary,gain_values,s,DataFolderName,FigureFolderName,folderName,nTries)