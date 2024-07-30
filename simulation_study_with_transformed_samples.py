# In this file, instead of generating points independently between the transcription and translation levels
# we generate the transcription levels and the use the mass action transform on them to obtain multipliers for the
# final scenario with co-translational_complexes
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
def run_simulation_w_points(expressions_to_vary,multiplier_values,simModel,DataFolderName,FigureFolderName,folderName,nTries=1000):
    print('Expressions to vary: ', expressions_to_vary)
    print('Using pre-generated conductance multipliers')
    # expression gains reset to match the baseline model operation
    expression_gains = [1, 1, 1, 1, 1]
    # simgas of lognormal distributions from which the gains will be sampled
    all_sigmas = ['NA', 'NA', 'NA', 'NA', 'NA']  # this is only for storing
    # get the indices of the currents to vary
    indices = [currentNames.index(expression_to_vary) for expression_to_vary in expressions_to_vary]
    # create a dictionary to store the simulation results
    allKeys = ['Index','Time'] + currentNames + ['Voltage']
    simulationResults = dict.fromkeys(allKeys)
    simRes_pro_arrhythmic = simulationResults.copy()
    # this dictionary will contain the lists of values for each simulation index
    for key in simulationResults.keys():
        simulationResults[key] = []
        simRes_pro_arrhythmic[key] = []
    # counters
    countAlternans = 0
    # initialise biomarker dictionary
    biomarkerNames = ['APD90', 'APD60', 'APD40', 'APA', 'TRI60', 'TRI40', 'EAD', '90PERCENT', '90to90', 'Alternan','Pro-Arrhythmic']
    columnNames = ['Index'] + gainNames + biomarkerNames
    biomarkers = pd.DataFrame(columns=columnNames)
    # figure for plotting currents
    nThingsToPlot = 1 + len(currentNamesToPlot)
    nRows = int(np.ceil(nThingsToPlot / 3))
    height = 5
    fig, axs = plt.subplots(nRows, 3, figsize=(height*3, height*nRows))
    axs = axs.flatten()
    # loop over the number of simulations
    for iTry in tqdm(range(nTries)):
        # read the gains on conductance from the list
        gain_values = multiplier_values[iTry]
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
        last_four_start = np.where(times > t_last_four)[0][0] - 1  # include one time point before because we dont have regular time increments
        last_pace_start = np.where(times > t_start_storing)[0][0] - 1
        # keep only the last 4 paces
        times_last_four = times[last_four_start:]
        times_last_pace = times[last_pace_start:]
        V_last_four_paces = d['membrane.v'][last_four_start:]
        # plot the simulation outputs
        axs[0].plot(times_last_four / 1000, V_last_four_paces , label=str(iTry+1))  # plot the voltage separately
        for iCurrent, currentName in enumerate(currentNamesInMMT):
            if currentName in currentNamesToPlot:
                # plot the current
                axs[currentNamesToPlot.index(currentName) + 1].plot(times_last_four / 1000, d[currentName][last_four_start:])
            # store the current values in the dictionary
            simulationResults[currentNames[iCurrent]].append(d[currentName][last_pace_start:])
        # store the rest of the simulation results
        simulationResults['Index'].append(iTry)
        simulationResults['Time'].append(times_last_pace)
        simulationResults['Voltage'].append(d['membrane.v'][last_pace_start:])
        del d
        gc.collect()
        ################################################################################################################
        # analysing the AP traces
        # split the last four paces into two groups by finding indexes of elements in times
        APs = []
        AP_times = []
        AP_interpolators = []
        for i_beat in range(len(times_of_beats) - 1):
            indexes = np.where((times_last_four >= t_end - times_of_beats[i_beat]) & (times_last_four < t_end - times_of_beats[i_beat + 1]))[0]
            APs.append(V_last_four_paces[indexes[0]:indexes[-1]])
            AP_times.append(times_last_four[indexes[0]:indexes[-1]])
            times_of_interp = AP_times[i_beat] - (t_end - times_of_beats[i_beat])
            bs_function = sp.interpolate.CubicSpline(times_of_interp, APs[i_beat],extrapolate=False)
            AP_interpolators.append(bs_function(np.linspace(0.1, 900, 1000)))  # interpolate at 1 ms intervals
        # check for alternans
        isAnAlternan = False
        isProArrhythmic = False
        for iTrace, Trace in enumerate(AP_interpolators):
            for jTrace in range(iTrace + 1, len(AP_interpolators)):
                # сheck if the difference between the traces is greater than 1 mV
                if np.any(np.abs(np.diff(Trace - AP_interpolators[jTrace])) >= 1):
                    isAnAlternan = True
                    isProArrhythmic = True
                    countAlternans += 1
                    newRow = [iTry] + expression_gains + [np.nan] * len(biomarkerNames[:-2]) + [isAnAlternan, isProArrhythmic]
                    biomarkers.loc[iTry] = newRow
                    print('Alternan detected. Total number of alternans: ', countAlternans)
                    break
            if isAnAlternan:
                break # once we found a difference between beats we can stop checking
        if not isAnAlternan:
            # print(not isAnAlternan)
            # use the last beat to calculate biomarkers
            iTrace = -1
            Trace = APs[iTrace]
            # check if early afterdepolarisation was triggered
            peaks, peak_properties = sp.signal.find_peaks(Trace, height=0)
            EAD = np.NaN
            if len(peaks) > 1:
                time_between_peaks = AP_times[iTrace][peaks[-1]] - AP_times[iTrace][peaks[0]]
                # if the time between peaks is greater than 50 ms we have an EAD
                if time_between_peaks >= 50:
                    EAD = time_between_peaks
                    isProArrhythmic = True
            #  find maximum voltage and its index
            maxVoltage = np.max(Trace)
            maxVoltageIndex = np.where(Trace == maxVoltage)[0][0]
            # find the minimum voltage and its index after the peak
            minVoltageAfterPeak = np.min(Trace[maxVoltageIndex:])
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
            newRow = [iTry] + expression_gains + [APD90, APD60, APD40, maxVoltage - minVoltageAfterPeak, triangulationSlope60, triangulationSlope40, EAD, repolarisationVoltage90, APD90to90, isAnAlternan, isProArrhythmic]
            biomarkers.loc[iTry] = newRow
        ################################################################################################################
        # stupid way to store stuff in the dictionary but we have timeseries in each entry
        if isProArrhythmic:
            # populate the proarrythmic dictionaries
            for key in simulationResults.keys():
                simRes_pro_arrhythmic[key].append(simulationResults[key][-1])
        # reset the simulation
        simModel.reset()
    #  end for over simulations
    ####################################################################################################################
    # plot the last four paces
    for ax in axs:
        ax.set_xlabel('time,s')
        ax.set_ylabel(ylabels_of_varied[axs.tolist().index(ax)])
    plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
    figName = FigureFolderName + '/' + folderName + '/'
    for ax in axs:
        ax.set_xlim((t_last_four) / 1000, t_end / 1000)
    plt.savefig(figName + 'last_four_paces.png', dpi=300)
    for ax in axs:
        ax.set_xlim((t_end - 1000) / 1000, t_end / 1000)
    if isATest:
        axs[0].legend() # show which tries have alternans
    plt.savefig(figName + 'last_pace.png', dpi=300)
    del fig
    # create a figure and make subplots for each biomarker as scutter plots
    fig1, axs1 = plt.subplots(2, 4, figsize=(15, 10))
    axs1 = axs1.ravel()
    # plot the biomarkers
    iAxs = 0
    # store only dictionary values for which pro-arrhythmicity was not detected
    well_behaved = biomarkers[biomarkers['Pro-Arrhythmic']==False]
    pro_arrhythmic = biomarkers[biomarkers['Pro-Arrhythmic']==True]
    for iBiomarker, biomarker in enumerate(biomarkerNames):
        if biomarker == 'Alternan' or biomarker == 'EAD' or biomarker == 'Pro-Arrhythmic':
            continue
        for iGain, gainName in enumerate(gains_to_vary):
            axs1[iAxs].scatter(well_behaved[gainName], well_behaved[biomarker], s=5, label='Multiplier on '+ gain_labels[iGain])
            axs1[iAxs].scatter(pro_arrhythmic[gainName], pro_arrhythmic[biomarker], s=10, color='magenta')
            if iGain==len(gains_to_vary):
                axs1[iAxs].scatter(pro_arrhythmic[gainName], pro_arrhythmic[biomarker], s=5, color='k',label='Pro-Arrhythmic')
        axs1[iAxs].set_xlabel(r'Multiplier  values')
        axs1[iAxs].set_ylabel(biomarker)
        axs1[iAxs].set_xscale('symlog', base=2, linthresh=0.125)
        axs1[iAxs].set_xticks([0.25, 0.5, 1, 2, 4, 8])
        if iAxs == 2:
            axs1[iAxs].legend()
        iAxs += 1
    fig1.suptitle(r'Biomarkers of AP traces obtained from gain covariation', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
    plt.savefig(figName + 'biomarkers.png', dpi=300)
    ####################################################################################################################
    # plot histograms of the biomarkers
    fig2, axs2 = plt.subplots(2, 4, figsize=(15, 10))
    axs2 = axs2.ravel()
    iAxs = 0
    for iBiomarker, biomarker in enumerate(biomarkerNames):
        if biomarker == 'Alternan' or biomarker == 'EAD' or biomarker == 'Pro-Arrhythmic':
            continue
        axs2[iAxs].hist(biomarkers[biomarker], bins=20, density=True, histtype='step')
        # add mean, median and std to the plot
        axs2[iAxs].axvline(np.nanmean(biomarkers[biomarker]), color='k', linestyle='solid', linewidth=1)
        axs2[iAxs].axvline(np.nanmedian(biomarkers[biomarker]), color='k', linestyle='dashed', linewidth=1)
        axs2[iAxs].axvline(np.nanmean(biomarkers[biomarker]) + np.nanstd(biomarkers[biomarker]), color='grey',
                           linestyle='dotted', linewidth=1)
        axs2[iAxs].axvline(np.nanmean(biomarkers[biomarker]) - np.nanstd(biomarkers[biomarker]), color='grey',
                           linestyle='dotted', linewidth=1)
        ymin, ymax = axs2[iAxs].get_ylim()
        axs2[iAxs].text(np.nanmean(biomarkers[biomarker]), ymax * 0.9,
                        'mean: {:.2f}'.format(np.nanmean(biomarkers[biomarker])))
        axs2[iAxs].text(np.nanmedian(biomarkers[biomarker]), ymax * 0.8,
                        'md: {:.2f}'.format(np.nanmedian(biomarkers[biomarker])))
        axs2[iAxs].text(np.mean(biomarkers[biomarker]) + np.nanstd(biomarkers[biomarker]), ymax * 0.7,
                        'std: {:.2f}'.format(np.nanstd(biomarkers[biomarker])))
        axs2[iAxs].set_xlabel(biomarker)
        iAxs += 1
    fig2.suptitle(r'Histogramms of biomarkers obtained from gain covariation', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
    plt.savefig(figName + 'biomarker_hists.png', dpi=300)
    ####################################################################################################################
    # plot the cloud of conductance multiplier values and the voltage traces - this will only be plotted for the first two conductances in expressions to vary
    fig3, axs3 = plt.subplots(1, 2, figsize=(10, 5))
    axs3 = axs3.ravel()
    axs3[0].scatter(biomarkers[gains_to_vary[0]], biomarkers[gains_to_vary[1]], s=5, color='k', alpha=0.5,label='Stable')
    axs3[0].scatter(pro_arrhythmic[gains_to_vary[0]], pro_arrhythmic[gains_to_vary[1]], s=15, marker='x', color='magenta',label='Pro-Arrhythmic')
    axs3[0].set_xlabel(gain_labels[0])
    axs3[0].set_ylabel(gain_labels[1])
    axs3[0].set_xscale('symlog', base=2, linthresh=0.125)
    axs3[0].set_xticks([0.25, 0.5, 1, 2, 4, 8])
    axs3[0].set_yscale('symlog', base=2, linthresh=0.125)
    axs3[0].set_yticks([0.25, 0.5, 1, 2, 4, 8])
    axs3[0].legend(loc='best')
    for iTry in range(nTries):
        axs3[1].plot(simulationResults['Time'][iTry] / 1000, simulationResults['Voltage'][iTry], color='k')
    for iArrythmic in range(len(simRes_pro_arrhythmic['Index'])):
        axs3[1].plot(simRes_pro_arrhythmic['Time'][iArrythmic] / 1000, simRes_pro_arrhythmic['Voltage'][iArrythmic], color='magenta')
    axs3[1].set_xlabel('time, s')
    axs3[1].set_ylabel('Voltage, mV')
    plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
    plt.savefig(figName + 'conductance_multipliers.png', dpi=300)
    ########################################################################################################################
    # save the simulation results and biomarkers into pickle files
    if not os.path.exists(DataFolderName + '/' + folderName):
        os.makedirs(DataFolderName + '/' + folderName)
    pickle.dump(simulationResults, open(DataFolderName + '/' + folderName + '/simulationResults.pkl', 'wb'))
    pickle.dump(simRes_pro_arrhythmic, open(DataFolderName + '/' + folderName + '/pro_arrhythmic_results.pkl', 'wb'))
    fileName = DataFolderName + '/' + folderName + '/biomarkers.csv'
    biomarkers.to_csv(fileName, index=False)
    # save the biomarkers that are pro-arrhythmic
    fileName = DataFolderName + '/' + folderName + '/pro_arrhythmic_biomarkers.csv'
    pro_arrhythmic.to_csv(fileName, index=False)
    # pickle.dump(pro_arrhythmic, open(DataFolderName + '/' + folderName + '/pro_arrhythmic_biomarkers.pkl', 'wb'))
    del simulationResults, fig1, fig2, axs1, axs2, biomarkers, pro_arrhythmic, simRes_pro_arrhythmic
    gc.collect()
    return True

# main
if __name__ == '__main__':
    # are we in the test mode?
    isATest = False #turning this on will save results into a test folder and run only 10 simulations
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
    if isATest:
        nTries = 10  # number of simulations
    else:
        nTries = 1000
    # OHR model takes about 450 paces to equilibrate - from Mann et al. 2016
    t_end = 454000  # end time of the simulation - there could be a smarter way to define this in terms of paces
    times_of_beats = [4000, 3000, 2000, 1000, 0]  # these are in descending order because they will be substraced from t_end
    # for storage only: record the last 4 paces
    t_start_storing = t_end - 1000
    t_last_four = t_end - 4000
    # baseline conductances of the currents of interest
    currentNames = ['IKr', 'IKs', 'Ito', 'INa', 'ICaL']
    currentNamesInMMT = ['IKr.IKr', 'IKs.IKs', 'Ito.Ito', 'INa.INa', 'ICaL.ICaL']
    currentNamesToPlot = ['IKr.IKr', 'ICaL.ICaL']
    conductanceNamesInMMT = ['IKr.GKr_b', 'IKs.GKs_b', 'Ito.Gto_b', 'INa.GNa', 'ICaL.PCa_b']
    gainNames = ['gain_kr', 'gain_ks', 'gain_to', 'gain_na', 'gain_ca']
    baselineConductances = [0.0321, 0.0011, 0.16, 11.7802, 8.37570000000000046e-5]
    print('baseline conductances before modification: ', baselineConductances)
    ylabels = [r'$V$, mV', r'$I_{Kr}, A/F$, ', r'$I_{Ks}$, A/F', r'$I_{to}$, A/F', r'$I_{Na}$, A/F', r'$I_{CaL}$, A/F']
    gain_labels  = [r'$g_{Kr}$', r'$g_{Ks}$', r'$g_{to}$', r'$g_{Na}$', r'$g_{CaL}$']
    gene_names = ['KCNH2', 'KCNQ1', 'KCND3', 'SCN5A', 'CACNA1C']
########################################################################################################################
# select which expression levels we are going to vary
########################################################################################################################
    # always include IKr as we know all ratios for it
    expressions_to_vary = ['IKr', 'ICaL']
    iExpressions = [currentNames.index(expression) for expression in expressions_to_vary]
    gain_labels = [gain_labels[i] for i in iExpressions]
    gains_to_vary = [gainNames[i] for i in iExpressions]
    ylabels_of_varied = [ylabels[0]] + [ylabels[i+1] for i in iExpressions]
    gene_names = [gene_names[i] for i in iExpressions[1:]]
########################################################################################################################
        # independently varying expression of I_Kr and I_CaL
########################################################################################################################
    print('Independently varied conductances for I_Kr and I_CaL')
    # select which conductances we wish to vary
    # create the folder for storing the results
    folderName = expressions_to_vary[0] + '_' + expressions_to_vary[1] + '_independent'
    if isATest:
        folderName = folderName + '_test'
    if not os.path.exists(DataFolderName + '/' + folderName):
        os.makedirs(DataFolderName + '/' + folderName)
    if not os.path.exists(FigureFolderName + '/' + folderName):
        os.makedirs(FigureFolderName + '/' + folderName)
    covars = pickle.load(open('Pickles/gtex_covariances.pkl', 'rb'))
    Sigma_bivariate = covars[gene_names[0]]
    # we only want the diagonal element from the
    if isATest:
        sigmas = [3, 3] # for testing pro-arrythimic highlighting - make sure to have large sigmas to generate some alternans
    else:
        sigmas = [2*np.sqrt(Sigma_bivariate[0, 0]), 2*np.sqrt(Sigma_bivariate[1, 1])]
    print('Expressions to vary: ', expressions_to_vary)
    print('Sigmas: ', sigmas)
    # using the sigmas generate a list of pairs of lognormally distributed values
    gain_values = np.random.lognormal(mean=0.0, sigma=sigmas, size=(nTries, len(sigmas)))
    SimulationSuccess = run_simulation_w_points(expressions_to_vary, gain_values, s, DataFolderName, FigureFolderName,
                                                folderName, nTries)
########################################################################################################################
    # covaried coexpression of I_Kr and I_CaL at transcription level
########################################################################################################################
    print('Covaried conductances for I_Kr and I_CaL')
    # create folder to store selected simulation outputs
    # baseline conductances of the currents of interest
    folderName = expressions_to_vary[0] + '_' + expressions_to_vary[1] + '_cotranscripted'
    if isATest:
        folderName = folderName + '_test'
    if not os.path.exists(DataFolderName + '/' + folderName):
        os.makedirs(DataFolderName + '/' + folderName)
    if not os.path.exists(FigureFolderName + '/' + folderName):
        os.makedirs(FigureFolderName + '/' + folderName)
    # load covariance from the pickle file created by the script 'read_gtex.py'
    covars = pickle.load(open('Pickles/gtex_covariances.pkl', 'rb'))
    Sigma_bivariate = covars[gene_names[0]]
    print(Sigma_bivariate)
    means = [0]*len(expressions_to_vary)
    sampled_from_normal = np.random.multivariate_normal(mean=means, cov=Sigma_bivariate, size=nTries)
    transcript_levels = np.exp(sampled_from_normal)
    SimulationSuccess = run_simulation_w_points(expressions_to_vary, transcript_levels, s, DataFolderName, FigureFolderName,
                                                folderName, nTries)
    ########################################################################################################################
    # covaried coexpression of I_Kr and I_CaL with cotransaltion covariance
    ########################################################################################################################
    print('Covaried conductances for I_Kr and I_CaL (converted from the co-transcription samples)')
    # create folder to store selected simulation outputs
    folderName = expressions_to_vary[0] + '_' + expressions_to_vary[1] + '_cotranslated'
    if isATest:
        folderName = folderName + '_test'
    if not os.path.exists(DataFolderName + '/' + folderName):
        os.makedirs(DataFolderName + '/' + folderName)
    if not os.path.exists(FigureFolderName + '/' + folderName):
        os.makedirs(FigureFolderName + '/' + folderName)
    # transform the cloud of transctipts
    l1_sample = l1_as_fn_of_k1k2_hERG_val(transcript_levels[:, 0], transcript_levels[:, 1])
    l2_sample = l2_as_fn_of_k1k2_SCN5A_val(transcript_levels[:, 0], transcript_levels[:, 1])
    protein_levels = np.concatenate((l1_sample.reshape(-1,1), l2_sample.reshape(-1,1)), axis=1)
    SimulationSuccess = run_simulation_w_points(expressions_to_vary,protein_levels,s,DataFolderName,FigureFolderName,folderName,nTries)
########################################################################################################################
    # Dependent variation of ion channel conductances
########################################################################################################################
    print('Dependent variation of conductances for I_Kr and I_CaL')
    # create folder to store selected simulation outputs
    folderName = expressions_to_vary[0] + '_' + expressions_to_vary[1] + '_dependent'
    if isATest:
        folderName = folderName + '_test'
    if not os.path.exists(DataFolderName + '/' + folderName):
        os.makedirs(DataFolderName + '/' + folderName)
    if not os.path.exists(FigureFolderName + '/' + folderName):
        os.makedirs(FigureFolderName + '/' + folderName)
    # sample gains from univariate lognormal distribution
    gain_values = np.random.lognormal(mean=0.0, sigma=sigmas[0], size=(nTries, 1))
    gain_values = np.concatenate((gain_values, gain_values), axis=1)
    SimulationSuccess = run_simulation_w_points(expressions_to_vary,gain_values,s,DataFolderName,FigureFolderName,folderName,nTries)

