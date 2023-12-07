import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
# import pandas as pd
import myokit as mk
import os
# import csv
import gc
import pickle
from tqdm import tqdm
matplotlib.use('Agg') # to turn off interractive mode

# main
if __name__ == '__main__':
    # if there is no folder for figures, create one
    FigureFolderName = 'Figures_Ballouz'
    DataFolderName = 'Simulated_data_Ballouz'
    if not os.path.exists(FigureFolderName):
        os.makedirs(FigureFolderName)
    # if there is no folder for output, create one
    if not os.path.exists(DataFolderName):
        os.makedirs(DataFolderName)
    # load the model and protocol
    model, protocol, embedded_script = mk.load('OHara_cellml/OHR_chaning_conductances.mmt')
    currentNamesInMMT = ['IKr.IKr', 'IKs.IKs', 'Ito.Ito', 'INa.INa', 'ICaL.ICaL']
    s = mk.Simulation(model, protocol)
    # d = s.run(1000)
    # s.reset()      # resets the state of the simulation
    vars_to_log = ['environment.time', 'membrane.v'] + currentNamesInMMT
    d_OHR_basic = s.run(1000,log=vars_to_log)
    # change the constants to mean values suggested in Tbale 1 in Mann et al. 2016
    ohr_baseline_conductances = []
    multipliers_from_mann = [1.00, 5.75, 2.01, 1.00, 2.95, 9.12] #the ones corresponding to conductances in the paper
    multipliers_from_krogh = [1.17, 8.09, 3.57, 3.05, 1.7, 1.91]
    additional_multipliers_from_mann = [4.48, 1.55, 4.69] # the ones denoted by S varibales in the paper
    constant_names_in_OHR = ['IKr.GKr_b', 'IKs.GKs_b','ICaL.PCa_b', 'INaL.GNaL_b', 'INaCa_i.Gncx_b', 'INaK.Pnak_b']
    # baseline_conductances_OHR = [4.65854545454545618e-2, 6.35800000000000080e-3, 0.0001007, 1.99574999999999753e-2, 0.0008, 30] # taken from the mmt file
    baseline_conductances_OHR = [0.046, 0.0034, 0.0001, 0.0075, 0.0008,30] # taken directly from Rudy Lab
    # simulate baseline model with params from Rudy lab
    # get new conductances by multiplying the baseline conductances by the multipliers
    modified_conductances_OHR = []
    for i in range(len(baseline_conductances_OHR)):
        modified_conductances_OHR.append(baseline_conductances_OHR[i] * multipliers_from_mann[i])
    # set the new conductances in the model
    for iConductance, conductance in enumerate(modified_conductances_OHR):
        s.set_constant(constant_names_in_OHR[iConductance], conductance)
    # run the simulation
    s.reset()  # resets the state of the simulation
    d_OHR_modified = s.run(1000,log=vars_to_log)
    # try with the krogh values for multipliers
    modified_conductances_OHR = []
    for i in range(len(baseline_conductances_OHR)):
        modified_conductances_OHR.append(baseline_conductances_OHR[i] * multipliers_from_krogh[i])
    # set the new conductances in the model
    for iConductance, conductance in enumerate(modified_conductances_OHR):
        s.set_constant(constant_names_in_OHR[iConductance], conductance)
    # run the simulation
    s.reset()      # resets the state of the simulation
    d_OHR_krogh = s.run(1000,log=vars_to_log)
    # further modify the conductances of excangers - placeholder fo this here!!


    # plot the results for the original and modified model
    currentNamesInMMT = ['IKr.IKr', 'IKs.IKs', 'Ito.Ito', 'INa.INa', 'ICaL.ICaL']
    ylabels = [r'$V$, mV', r'$I_{Kr}$, ', r'$I_{Ks}$, A/F', r'$I_{to}$, A/F', r'$I_{Na}$, A/F', r'$I_{CaL}$, A/F']
    fig, axs = plt.subplots(2, 3, figsize=(20, 10))
    axs = axs.flatten()
    axs[0].plot(d_OHR_basic.time(), d_OHR_basic['membrane.v'])  # plot the voltage separately
    axs[0].plot(d_OHR_modified.time(), d_OHR_modified['membrane.v'])
    axs[0].plot(d_OHR_krogh.time(), d_OHR_krogh['membrane.v'])
    for iCurrent, currentName in enumerate(currentNamesInMMT):
        # plot the current
        if iCurrent == 1:
            axs[iCurrent + 1].plot(d_OHR_basic.time(), d_OHR_basic[currentName], label="OHR original")
            axs[iCurrent + 1].plot(d_OHR_modified.time(), d_OHR_modified[currentName], label="Mann et.al. multipliers")
            axs[iCurrent + 1].plot(d_OHR_krogh.time(), d_OHR_krogh[currentName], label="Krogh-Madsen et.al.  multipliers")
        else:
            axs[iCurrent + 1].plot(d_OHR_basic.time(), d_OHR_basic[currentName])
            axs[iCurrent + 1].plot(d_OHR_modified.time(), d_OHR_modified[currentName])
            axs[iCurrent + 1].plot(d_OHR_krogh.time(), d_OHR_krogh[currentName])
    for ax in axs:
        ax.set_xlabel('time, s')
        ax.set_ylabel(ylabels[axs.tolist().index(ax)])
    # legend for one of the axes only
    axs[2].legend()
    fig.suptitle(r'Original OHara-Rudy model vs Modifications', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(FigureFolderName+'/OHR_models_one_pace.png', dpi=300)
    # plt.show()


    ## try altering conductances in the model
    model, protocol, embedded_script = mk.load('OHara_cellml/OHR_chaning_conductances.mmt')
    ########################################################################################################################
    # set up the simulation - this is the same for all configurations
    s = mk.Simulation(model, protocol)
    ## these variables are shared by all other configurations and should not be reset
    nTries = 500 # number of simulations
    # OHR model takes about 450 paces to equilibrate - from Mann et al. 2016
    t_end = 454000 # end time of the simulation - there could be a smarter way to define this in terms of paces
    times_of_beats = [4000, 3000, 2000, 1000, 0]  # these are in descending order because they will be substraced from t_end
    # for storage only: record the last 4 paces
    t_start_storing = t_end - 4000
    # baseline conductances of the currents of interest
    currentNames = ['IKr', 'IKs', 'Ito', 'INa', 'ICaL']
    currentNamesInMMT = ['IKr.IKr', 'IKs.IKs', 'Ito.Ito', 'INa.INa', 'ICaL.ICaL']
    conductanceNamesInMMT = ['IKr.GKr_b', 'IKs.GKs_b', 'Ito.Gto_b', 'INa.GNa', 'ICaL.PCa_b']
    gainNames = ['gain_kr', 'gain_ks', 'gain_to', 'gain_na', 'gain_ca']
    baselineConductances = [4.65854545454545618e-2, 6.35800000000000080e-3, 0.02, 75, 0.0001007]
    print('baseline conductances before modification: ', baselineConductances)
    # in the baselineConductances list, replace the values by the value in modified_conductances_OHR
    for iConductance, conductanceName in enumerate(conductanceNamesInMMT):
        if conductanceName in constant_names_in_OHR:
            baselineConductances[iConductance] = modified_conductances_OHR[constant_names_in_OHR.index(conductanceName)]
    print('baseline conductances after modification: ', baselineConductances)
    ylabels = [r'$V$, mV', r'$I_{Kr}$, ', r'$I_{Ks}$, A/F', r'$I_{to}$, A/F', r'$I_{Na}$, A/F', r'$I_{CaL}$, A/F']
########################################################################################################################
    # # try changing conductances and pass them into the file
    # fig, axs = plt.subplots(2, 3, figsize=(20, 10))
    # axs = axs.flatten()
    # for gain_value in np.arange(0.5, 1.5, 0.1):
    #     # get the list of conductances with the new expression gains
    #     conductances_sampled = [baselineConductances[i] * gain_value for i in range(len(baselineConductances))]
    #     # set constants in the mmt file to new values
    #     for iConductance, conductance in enumerate(conductances_sampled):
    #         s.set_constant(conductanceNamesInMMT[iConductance], conductance)
    #     # run the simulation
    #     d = s.run(t_end,log=vars_to_log)
    #     # get the times as array
    #     times = np.array(d.time())
    #     # plot the simulation outputs
    #     axs[0].plot(times/1000, d['membrane.v'])  # plot the voltage separately
    #     for iCurrent, currentName in enumerate(currentNamesInMMT):
    #         # plot the current
    #         if iCurrent==1:
    #             axs[iCurrent + 1].plot(times/1000, d[currentName], label='gain = %.3f' % gain_value)
    #         else:
    #             axs[iCurrent + 1].plot(times/1000, d[currentName])
    #     # reset the simulation
    #     s.reset()
    # #  end for over simulations
    # for ax in axs:
    #     ax.set_xlabel('time, s')
    #     ax.set_ylabel(ylabels[axs.tolist().index(ax)])
    # # legend for one of the axes only
    # axs[2].legend(bbox_to_anchor=(0.9, 0.2, 0.2, 1), loc="lower left",
    #               borderaxespad=0, ncol=1)
    # fig.suptitle(r'Changing all conductances simulteneously', fontsize=14)
    # plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
    # figName = FigureFolderName+'/all_conductances_changing'
    # for ax in axs:
    #     ax.set_xlim(t_start_storing/1000, t_end/1000)
    # plt.savefig(figName + '_last_four_paces.png', dpi=300)
    # for ax in axs:
    #     ax.set_xlim((t_end - 1000)/1000, t_end/1000)
    # plt.savefig(figName + '_last_pace.png', dpi=300)
    # plt.show()
    # print('pause here')
########################################################################################################################
    ## from here onward things can be changed
    # independently varying expression of I_Kr and I_CaL
    print('Independently varied conductances for I_Kr and I_CaL')
    # create folder to store selected simulation outputs
    folderName = 'Kr_CaL_independent'
    # expression gains reset to match the baseline model operation
    expression_gains = [1, 1, 1, 1, 1]
    # simgas of lognormal distributions from which the gains will be sampled
    all_sigmas = ['NA', 'NA', 'NA', 'NA', 'NA']  # this is only for storing
    # baseline conductances of the currents of interest
    expressions_to_vary = ['IKr', 'ICaL']
    sigmas = [np.sqrt(0.25), np.sqrt(0.25)]
    print('Expressions to vary: ', expressions_to_vary)
    print('Sigmas: ', sigmas)
    # get the indices of the currents to vary
    indices = [currentNames.index(expression_to_vary) for expression_to_vary in expressions_to_vary]
    # assign the sigmas to the corresponding currents
    for iSigma, sigma in enumerate(sigmas):
        all_sigmas[indices[iSigma]] = sigma
    # create a dictionary to store the simulation results
    allKeys = gainNames + ['Time'] + currentNames + ['Voltage']
    simulationResults = dict.fromkeys(allKeys)
    # this dictionary will contain the lists of values for each simulation index
    for key in simulationResults.keys():
        simulationResults[key] = []
    # counters
    countAlternans = 0
    # initialise biomarker dictionary
    biomarkers = dict.fromkeys(['APD90', 'APD60', 'APD40', 'APA', 'TRI60', 'TRI40', 'EAD', 'Alternan'])
    # initialise the biomarker dictionary with empty lists
    for key in biomarkers.keys():
        biomarkers[key] = []
    # figure for plotting currents
    fig, axs = plt.subplots(2, 3, figsize=(20, 10))
    axs = axs.flatten()
    # loop over the number of simulations
    for iTry in  tqdm(range(nTries)):
        # sample gains from lognormal distributions for the currents of interest
        for iSigma, sigma in enumerate(sigmas):
            if sigma != 'NA':
                gain_value = np.random.lognormal(mean=0.0, sigma=sigma, size=1)[0]
                expression_gains[indices[iSigma]] = gain_value
        # get the list of conductances with the new expression gains
        conductances_sampled = [baselineConductances[i] * expression_gains[i] for i in range(len(baselineConductances))]
        # set constants in the mmt file to new values
        for iConductance, conductance in enumerate(conductances_sampled):
            s.set_constant(conductanceNamesInMMT[iConductance], conductance)
        # run the simulation
        d = s.run(t_end,log=vars_to_log)
        # get the times as array
        times = np.array(d.time())
        # find the index of d.time() from which the last 4 paces start
        last_paces_start = np.where(times > t_start_storing)[0][0] - 1 # include one time point before because we dont have regular time increments
        # print('time at which storing starts: ', times[last_paces_start])
        # keep only the last 4 paces
        times = times[last_paces_start:]
        V_last_four_paces =  d['membrane.v'][last_paces_start:]
        # plot the simulation outputs
        axs[0].plot(times/1000, V_last_four_paces) # plot the voltage separately
        for iCurrent, currentName in enumerate(currentNamesInMMT):
            # plot the current
            axs[iCurrent+1].plot(times/1000, d[currentName][last_paces_start:])
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
        for i_beat in range(len(times_of_beats)-1):
            indexes = np.where((times >= t_end - times_of_beats[i_beat]) & (times < t_end - times_of_beats[i_beat+1]))[0]
            APs.append(V_last_four_paces[indexes[0]-1:indexes[-1]])
            AP_times.append(times[indexes[0]-1:indexes[-1]])
            bs_function = sp.interpolate.CubicSpline(AP_times[i_beat]-(t_end-times_of_beats[i_beat]), APs[i_beat])
            AP_interpolators.append(bs_function(np.arange(10, 900, 0.1))) # interpolate at 1 ms intervals
        # check for alternans
        isAnAlternan = False
        for iTrace,Trace in enumerate(AP_interpolators):
            for jTrace in range(iTrace+1,len(AP_interpolators)):
                # сheck if the difference between the traces is greater than 1 mV
                if np.any(np.abs(np.diff(Trace-AP_interpolators[jTrace]))>=1):
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
            secondDerivative = np.diff(np.diff(Trace))
            d2Vdt2 = np.diff(np.diff(Trace)) /  (np.diff(AP_times[iTrace])[1:]**2)
            # placeholder for finding EADs
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
            # find APD90, APD60, APD40
            APD90 = AP_times[iTrace][repolarisationIndex90] - AP_times[iTrace][maxVoltageIndex]
            APD60 = AP_times[iTrace][repolarisationIndex60] - AP_times[iTrace][maxVoltageIndex]
            APD40 = AP_times[iTrace][repolarisationIndex40] - AP_times[iTrace][maxVoltageIndex]
            # compute trianglation slopes as average slope of the trace from APD40 and APD60 to APD90
            triangulationSlope40 = np.mean(np.diff(Trace[repolarisationIndex40:repolarisationIndex90])/np.diff(AP_times[iTrace][repolarisationIndex40:repolarisationIndex90]))
            triangulationSlope60 = np.mean(np.diff(Trace[repolarisationIndex60:repolarisationIndex90])/np.diff(AP_times[iTrace][repolarisationIndex60:repolarisationIndex90]))
            # store the biomarkers
            biomarkers['APD90'].append(APD90)
            biomarkers['APD60'].append(APD60)
            biomarkers['APD40'].append(APD40)
            biomarkers['TRI60'].append(triangulationSlope60)
            biomarkers['TRI40'].append(triangulationSlope40)
            biomarkers['Alternan'].append(isAnAlternan)
        ################################################################################################################
        # reset the simulation
        s.reset()
    #  end for over simulations
    for ax in axs:
        ax.set_xlabel('time,s')
        ax.set_ylabel(ylabels[axs.tolist().index(ax)])
    fig.suptitle(r'Changing $I_{Kr}$ and $I_{CaL}$ conductances independently',fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
    figName = FigureFolderName+'/Kr_and_CaL_gains_ind_variation_sigma2_025'
    for ax in axs:
        ax.set_xlim((t_start_storing)/1000, t_end/1000)
    plt.savefig(figName + '_last_four_paces.png', dpi=300)
    for ax in axs:
        ax.set_xlim((t_end - 1000)/1000, t_end/1000)
    plt.savefig(figName + '_last_pace.png', dpi=300)
    # plt.show()
    del fig
    # create a figure and make subplots for each biomarker as scutter plots
    fig1, axs1 = plt.subplots(2, 3, figsize=(15, 10))
    axs1 = axs1.ravel()
    # plot the biomarkers
    iAxs = 0
    for iBiomarker, biomarker in enumerate(biomarkers.keys()):
        if biomarker == 'Alternan' or biomarker == 'EAD':
            continue
        axs1[iAxs].scatter(simulationResults['gain_kr'], biomarkers[biomarker], s=5, label=r'Gain on $I_{Kr}$')
        axs1[iAxs].scatter(simulationResults['gain_ca'], biomarkers[biomarker], s=5, label=r'Gain on $I_{CaL}$')
        axs1[iAxs].set_xlabel(r'Expression gain values')
        axs1[iAxs].set_ylabel(biomarker)
        if iAxs == 2:
            axs1[iAxs].legend()
        iAxs += 1
    fig1.suptitle(r'Biomarkers of AP traces obtained from independent gain variation', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
    plt.savefig(figName + '_biomarkers.png', dpi=300)
    # plot histograms of the biomarkers
    fig2, axs2 = plt.subplots(2, 3, figsize=(15, 10))
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
        axs2[iAxs].axvline(np.median(biomarkers[biomarker])+np.std(biomarkers[biomarker]), color='grey', linestyle='dotted', linewidth=1)
        axs2[iAxs].axvline(np.mean(biomarkers[biomarker])- np.std(biomarkers[biomarker]), color='grey',
                           linestyle='dotted', linewidth=1)
        ymin, ymax = axs2[iAxs].get_ylim()
        axs2[iAxs].text(np.mean(biomarkers[biomarker]), ymax*0.9, 'mean: {:.2f}'.format(np.mean(biomarkers[biomarker])))
        axs2[iAxs].text(np.median(biomarkers[biomarker]), ymax*0.8, 'md: {:.2f}'.format(np.median(biomarkers[biomarker])))
        axs2[iAxs].text(np.median(biomarkers[biomarker])+np.std(biomarkers[biomarker]), ymax*0.7, 'std: {:.2f}'.format(np.std(biomarkers[biomarker])))
        axs2[iAxs].set_xlabel(biomarker)
        iAxs += 1
    fig2.suptitle(r'Histogramms of biomarkers obtained from independent gain variation', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
    plt.savefig(figName + '_biomarker_hists.png', dpi=300)
########################################################################################################################
    # save the simulation results and biomarkers into pickle files
    if not os.path.exists(DataFolderName + '/' + folderName):
        os.makedirs(DataFolderName + '/' + folderName)
    pickle.dump(simulationResults, open(DataFolderName+ '/' + folderName + '/simulationResults.pkl', 'wb'))
    pickle.dump(biomarkers, open(DataFolderName+ '/' + folderName + '/biomarkers.pkl', 'wb'))
    del simulationResults, fig1, fig2, axs1, axs2, biomarkers
    gc.collect()

########################################################################################################################
    # dependent co-expression of I_Kr and I_CaL
########################################################################################################################
    # gain factors
    print('Jointly varied conductances for I_Kr and I_CaL')
    folderName = 'Kr_CaL_joint'
    # reset expression gains to baseline values
    expression_gains = [1, 1, 1, 1, 1]
    # currents of interest
    expressions_to_vary = ['IKr','ICaL']
    all_sigmas = ['NA', 'NA', 'NA', 'NA', 'NA']  # this is only for storing
    sigmas = [np.sqrt(0.25)] * 2
    print('Expressions to vary: ', expressions_to_vary)
    print('Sigmas: ', sigmas)
    # get the indices of the currents to vary
    indices = [currentNames.index(expression_to_vary) for expression_to_vary in expressions_to_vary]
    # assign the sigmas to the corresponding currents
    for iSigma, sigma in enumerate(sigmas):
        all_sigmas[indices[iSigma]] = sigma
    # create a dictionary to store the simulation results
    allKeys = gainNames + ['Time'] + currentNames + ['Voltage']
    simulationResults = dict.fromkeys(allKeys)
    # this dictionary will contain the lists of values for each simulation index
    for key in simulationResults.keys():
        simulationResults[key] = []
    # counters
    countAlternans = 0
    # initialise biomarker dictionary
    biomarkers = dict.fromkeys(['APD90', 'APD60', 'APD40', 'APA', 'TRI60', 'TRI40', 'EAD', 'Alternan'])
    # initialise the biomarker dictionary with empty lists
    for key in biomarkers.keys():
        biomarkers[key] = []
    # figure for plotting currents
    fig, axs = plt.subplots(2, 3, figsize=(20, 10))
    axs = axs.flatten()
    # loop over the number of simulations
    for iTry in tqdm(range(nTries)):
        # sample a single value from a lognormal distribution
        gain_value = np.random.lognormal(mean=0.0, sigma=sigmas[0], size=1)[0]
        # assign the gain value to the currents of interest
        for iSigma, sigma in enumerate(sigmas):
                expression_gains[indices[iSigma]] = gain_value
                # get the list of conductances with the new expression gains
        conductances_sampled = [baselineConductances[i] * expression_gains[i] for i in range(len(baselineConductances))]
        # set constants in the mmt file to new values
        for iConductance, conductance in enumerate(conductances_sampled):
            s.set_constant(conductanceNamesInMMT[iConductance], conductance)
        # run the simulation
        d = s.run(t_end, log=vars_to_log)
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
            # plot the current
            axs[iCurrent + 1].plot(times / 1000, d[currentName][last_paces_start:])
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
            APs.append(V_last_four_paces[indexes[0] - 1:indexes[-1]])
            AP_times.append(times[indexes[0] - 1:indexes[-1]])
            bs_function = sp.interpolate.CubicSpline(AP_times[i_beat] - (t_end - times_of_beats[i_beat]),
                                                     APs[i_beat])
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
            secondDerivative = np.diff(np.diff(Trace))
            d2Vdt2 = np.diff(np.diff(Trace)) / (np.diff(AP_times[iTrace])[1:] ** 2)
            # placeholder for finding EADs
            #  find maximum voltage and its index
            maxVoltage = np.max(Trace)
            maxVoltageIndex = np.where(Trace == maxVoltage)[0][0]
            biomarkers['APA'].append(maxVoltage)
            # find the minimum voltage and its index after the peak
            minVoltageAfterPeak = np.min(Trace[maxVoltageIndex:])
            minVoltageIndex = np.where(Trace == minVoltageAfterPeak)[0][0]
            # find the 90% repolarisation voltage
            repolarisationVoltage90 = 0.1 * (maxVoltage - minVoltageAfterPeak) + minVoltageAfterPeak
            repolarisationVoltage60 = 0.4 * (maxVoltage - minVoltageAfterPeak) + minVoltageAfterPeak
            repolarisationVoltage40 = 0.6 * (maxVoltage - minVoltageAfterPeak) + minVoltageAfterPeak
            # find the index at which the voltage crosses the 90% repolarisation voltage
            repolarisationIndex90 = np.where(Trace[maxVoltageIndex:] <= repolarisationVoltage90)[0][
                                        0] + maxVoltageIndex
            repolarisationIndex60 = np.where(Trace[maxVoltageIndex:] <= repolarisationVoltage60)[0][
                                        0] + maxVoltageIndex
            repolarisationIndex40 = np.where(Trace[maxVoltageIndex:] <= repolarisationVoltage40)[0][
                                        0] + maxVoltageIndex
            # find APD90, APD60, APD40
            APD90 = AP_times[iTrace][repolarisationIndex90] - AP_times[iTrace][maxVoltageIndex]
            APD60 = AP_times[iTrace][repolarisationIndex60] - AP_times[iTrace][maxVoltageIndex]
            APD40 = AP_times[iTrace][repolarisationIndex40] - AP_times[iTrace][maxVoltageIndex]
            # compute trianglation slopes as average slope of the trace from APD40 and APD60 to APD90
            triangulationSlope40 = np.mean(np.diff(Trace[repolarisationIndex40:repolarisationIndex90]) / np.diff(
                AP_times[iTrace][repolarisationIndex40:repolarisationIndex90]))
            triangulationSlope60 = np.mean(np.diff(Trace[repolarisationIndex60:repolarisationIndex90]) / np.diff(
                AP_times[iTrace][repolarisationIndex60:repolarisationIndex90]))
            # store the biomarkers
            biomarkers['APD90'].append(APD90)
            biomarkers['APD60'].append(APD60)
            biomarkers['APD40'].append(APD40)
            biomarkers['TRI60'].append(triangulationSlope60)
            biomarkers['TRI40'].append(triangulationSlope40)
            biomarkers['Alternan'].append(isAnAlternan)
        ################################################################################################################
        # reset the simulation
        s.reset()
    #  end for over simulations
    for ax in axs:
        ax.set_xlabel('time,s')
        ax.set_ylabel(ylabels[axs.tolist().index(ax)])
    fig.suptitle(r'Changing $I_{Kr}$ and $I_{CaL}$ conductances jointly', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
    figName = FigureFolderName+'/Kr_and_CaL_gains_joint_variation_sigma2_025'
    for ax in axs:
        ax.set_xlim(t_start_storing/1000, t_end/1000)
    plt.savefig(figName + '_last_four_paces.png', dpi=300)
    for ax in axs:
        ax.set_xlim((t_end - 1000)/1000, t_end/1000)
    plt.savefig(figName + '_last_pace.png', dpi=300)
    # plt.show()
    del fig
    # create a figure and make subplots for each biomarker as scutter plots
    fig1, axs1 = plt.subplots(2, 3, figsize=(15, 10))
    axs1 = axs1.ravel()
    # plot the biomarkers
    iAxs = 0
    for iBiomarker, biomarker in enumerate(biomarkers.keys()):
        if biomarker == 'Alternan' or biomarker == 'EAD':
            continue
        axs1[iAxs].scatter(simulationResults['gain_kr'], biomarkers[biomarker], s=5, label=r'Gain on $I_{Kr}$')
        axs1[iAxs].scatter(simulationResults['gain_ca'], biomarkers[biomarker], s=5, label=r'Gain on $I_{CaL}$')
        axs1[iAxs].set_xlabel(r'Expression gain values')
        axs1[iAxs].set_ylabel(biomarker)
        if iAxs == 2:
            axs1[iAxs].legend()
        iAxs += 1
    fig1.suptitle(r'Biomarkers of AP traces obtained from independent gain variation', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
    plt.savefig(figName + '_biomarkers.png', dpi=300)
    # plot histograms of the biomarkers
    fig2, axs2 = plt.subplots(2, 3, figsize=(15, 10))
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
        axs2[iAxs].axvline(np.median(biomarkers[biomarker])+np.std(biomarkers[biomarker]), color='grey', linestyle='dotted', linewidth=1)
        axs2[iAxs].axvline(np.mean(biomarkers[biomarker])- np.std(biomarkers[biomarker]), color='grey',
                           linestyle='dotted', linewidth=1)
        ymin, ymax = axs2[iAxs].get_ylim()
        axs2[iAxs].text(np.mean(biomarkers[biomarker]), ymax*0.9, 'mean: {:.2f}'.format(np.mean(biomarkers[biomarker])))
        axs2[iAxs].text(np.median(biomarkers[biomarker]), ymax*0.8, 'md: {:.2f}'.format(np.median(biomarkers[biomarker])))
        axs2[iAxs].text(np.median(biomarkers[biomarker])+np.std(biomarkers[biomarker]), ymax*0.7, 'std: {:.2f}'.format(np.std(biomarkers[biomarker])))
        axs2[iAxs].set_xlabel(biomarker)
        iAxs += 1
    fig2.suptitle(r'Histogramms of biomarkers obtained from independent gain variation', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
    plt.savefig(figName + '_biomarker_hists.png', dpi=300)
########################################################################################################################
    # save the simulation results and biomarkers into pickle files
    if not os.path.exists(DataFolderName + '/' + folderName):
        os.makedirs(DataFolderName + '/' + folderName)
    pickle.dump(simulationResults, open(DataFolderName+ '/' + folderName + '/simulationResults.pkl', 'wb'))
    pickle.dump(biomarkers, open(DataFolderName+ '/' + folderName + '/biomarkers.pkl', 'wb'))
    del simulationResults, fig1, fig2, axs1, axs2, biomarkers
    gc.collect()
