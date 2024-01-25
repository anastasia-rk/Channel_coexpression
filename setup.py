import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
# import pandas as pd
import myokit as mk
import os
import pandas as pd
import csv
import pickle
import gc
from tqdm import tqdm
matplotlib.use('Agg')

# main
if __name__ == '__main__':
    # if there is no folder for figures, create one
    FigureFolderName = 'Figures_test'
    DataFolderName = 'Simulated_data_test'
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
    filename = 'Tomek_fixed_cellml/test_output.csv'
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
    t_end = 1500000  # end time of the simulation - there could be a smarter way to define this in terms of paces
    nExtraPaces = 200  # extra time to run the simulation for after the last beat to check if the state has equilibrated
    t_extra = 1000  # extra time to run the simulation for after the last beat to check if the state has equilibrated
    times_of_beats = [0, 1000, 2000, 3000, 4000]  # these are in descending order because they will be substraced from t_end
    # for storage only: record the last 4 paces
    t_start_storing = t_end - 4000
    # baseline conductances of the currents of interest
    currentNames = ['IKr', 'IKs', 'Ito', 'INa', 'ICaL']
    currentNamesInMMT = ['IKr.IKr', 'IKs.IKs', 'Ito.Ito', 'INa.INa', 'ICaL.ICaL']
    conductanceNamesInMMT = ['IKr.GKr_b', 'IKs.GKs_b', 'Ito.Gto_b', 'INa.GNa', 'ICaL.PCa_b']
    gainNames = ['gain_kr', 'gain_ks', 'gain_to', 'gain_na', 'gain_ca']
    baselineConductances = [0.0321, 0.0011, 0.16, 11.7802, 8.37570000000000046e-5]
    print('baseline conductances before modification: ', baselineConductances)
    ylabels = [r'$V$, mV', r'$I_{Kr}$, ', r'$I_{Ks}$, A/F', r'$I_{to}$, A/F', r'$I_{Na}$, A/F', r'$I_{CaL}$, A/F']
########################################################################################################################
    # plot the pdf of the lognormal distribution with given sigma and zero  mean
    sigma = 0.5
    x = np.linspace(0, 5, 100)
    pdf = np.exp(-(np.log(x) - 0) ** 2 / (2 * sigma ** 2)) / (x * sigma * np.sqrt(2 * np.pi))
    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    ax.plot(x, pdf)
    # draw vertical line at x=0.5 and x=2
    ax.axvline(x=0.5, color='k', linestyle='dashed', linewidth=1)
    ax.axvline(x=2, color='k', linestyle='dashed', linewidth=1)
    # add values of pdf at x=0.5 and x=2 as labels on x axis
    plt.xlabel('x')
    plt.ylabel('pdf')
    plt.title('Lognormal distribution with sigma = %.2f' % sigma)
    plt.savefig(FigureFolderName + '/lognormal_sigma_%.2f.png' % sigma, dpi=300)
    #######################################################################################################################
    # # try changing conductances and pass them into the file
    # fig, axs = plt.subplots(2, 3, figsize=(20, 10))
    # axs = axs.flatten()
    # for gain_value in np.arange(0.5, 2, 0.15):
    #     # get the list of conductances with the new expression gains
    #     conductances_sampled = [baselineConductances[i] * gain_value for i in range(len(baselineConductances))]
    #     # set constants in the mmt file to new values
    #     for iConductance, conductance in enumerate(conductances_sampled):
    #         s.set_constant(conductanceNamesInMMT[iConductance], conductance)
    #     # run the simulation
    #     d = s.run(t_end,log=vars_to_log)
    #     t_end_new = t_end
    #     old_state = s.state()
    #     # run for extra 20 paces and check if the state has equilibrated
    #     for iPaces in range(nExtraPaces):
    #         t_end_new = t_end_new + t_extra
    #         d = s.run(t_extra, log=vars_to_log)
    #         d_states = np.array(s.state()) - np.array(old_state)
    #         if np.max(np.abs(d_states) < 1e-4):
    #             print('State equilibrated after ', iPaces+1, ' extra paces. Max(x_end - x_0): ', np.max(np.abs(d_states)))
    #             break
    #         else:
    #             old_state = s.state()
    #     # run for extra 4 paces to analyse them
    #     d = s.run(4*t_extra, log=vars_to_log)
    #     # get the times as array
    #     times = np.array(d.time()) - t_end_new
    #     # plot the simulation outputs
    #     axs[0].plot(times, d['membrane.v'])  # plot the voltage separately
    #     for iCurrent, currentName in enumerate(currentNamesInMMT):
    #         # plot the current
    #         if iCurrent==1:
    #             axs[iCurrent + 1].plot(times, d[currentName], label='gain = %.3f' % gain_value)
    #         else:
    #             axs[iCurrent + 1].plot(times, d[currentName])
    #     # reset the simulation
    #     del d
    #     s.reset()
    # #  end for over simulations
    # for ax in axs:
    #     ax.set_xlabel('time, ms')
    #     ax.set_ylabel(ylabels[axs.tolist().index(ax)])
    # # legend for one of the axes only
    # axs[2].legend(bbox_to_anchor=(0.9, 0.2, 0.2, 1), loc="lower left",
    #               borderaxespad=0, ncol=1)
    # fig.suptitle(r'Changing all conductances simulteneously', fontsize=14)
    # plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
    # figName = FigureFolderName+'/all_conductances_changing'
    # for ax in axs:
    #     ax.set_xlim(times[0],times[-1])
    # plt.savefig(figName + '_last_four_paces.png', dpi=300)
    # for ax in axs:
    #     ax.set_xlim(times[np.where(times > 3000)[0][0]],times[-1])
    # plt.savefig(figName + '_last_pace.png', dpi=300)
    # plt.show()
    print('pause here')
    #######################################################################################################################
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
    sigmas = [0.5, 0.5]
    print('Expressions to vary: ', expressions_to_vary)
    print('Sigmas: ', sigmas)
    # get the indices of the currents to vary
    indices = [currentNames.index(expression_to_vary) for expression_to_vary in expressions_to_vary]
    # assign the sigmas to the corresponding currents
    for iSigma, sigma in enumerate(sigmas):
        all_sigmas[indices[iSigma]] = sigma
    # create a dataframe to store the simulation results into csv file
    allKeys = ['Simulation'] +  gainNames + ['Time'] + currentNames + ['Voltage']
    all_simulations_df = pd.DataFrame(columns=allKeys)
    # store gains for plotting
    store_gains = dict.fromkeys(gainNames)
    # this dictionary will contain the lists of values for each simulation index
    for key in store_gains.keys():
        store_gains[key] = []
    # counters
    countAlternans = 0
    # initialise biomarker dictionary
    biomarkers = dict.fromkeys(['APD90', 'APD60', 'APD40', 'APA', 'TRI60', 'TRI40', 'EAD','90PERCENT','90to90','Alternan'])
    # initialise the biomarker dictionary with empty lists
    for key in biomarkers.keys():
        biomarkers[key] = []
    # figure for plotting currents
    fig, axs = plt.subplots(2, 3, figsize=(20, 10))
    axs = axs.flatten()
    # loop over the number of simulations
    for iTry in tqdm(range(nTries)):
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
        d = s.run(t_end, log=vars_to_log)
        t_end_new = t_end
        old_state = s.state()
        # run for extra 20 paces and check if the state has equilibrated
        for iPaces in range(nExtraPaces):
            t_end_new = t_end_new + t_extra
            d = s.run(t_extra, log=vars_to_log)
            d_states = np.array(s.state()) - np.array(old_state)
            if np.max(np.abs(d_states) < 1e-4):
                break
            else:
                old_state = s.state()
        # run for extra 4 paces to analyse them
        d = s.run(4 * t_extra, log=vars_to_log)
        # get the times as array
        times = np.array(d.time()) - t_end_new
        # find the index of d.time() from which the last 4 paces start
        V_last_four_paces = d['membrane.v']
        # plot the simulation outputs
        axs[0].plot(times / 1000, V_last_four_paces)  # plot the voltage separately
        # store the current values in the dataframe
        simulationResults_df = pd.DataFrame(columns=allKeys)
        for iCurrent, currentName in enumerate(currentNamesInMMT):
            # plot the current
            axs[iCurrent + 1].plot(times / 1000, d[currentName])
            # store the current name in simulationResults_df
            simulationResults_df[currentNames[iCurrent]] = d[currentName]
        del d
        gc.collect()
        # store the rest of the simulation results
        for iGain, gainName in enumerate(gainNames):
            store_gains[gainName].append(expression_gains[iGain])
            simulationResults_df[gainName] = expression_gains[iGain]*np.ones(len(times))
        simulationResults_df['Time'] = times
        simulationResults_df['Voltage'] = V_last_four_paces
        simulationResults_df['Simulation'] = iTry*np.ones(len(times))
        # add the simulation results to the dataframe all_simulations_df
        all_simulations_df = pd.concat([all_simulations_df,simulationResults_df])
        ################################################################################################################
        # analysing the AP traces
        # split the last four paces into two groups by finding indexes of elements in times
        APs = []
        AP_times = []
        AP_interpolators = []
        for i_beat in range(len(times_of_beats) - 1):
            indexes = np.where((times >= times_of_beats[i_beat]) & (times < times_of_beats[i_beat + 1]))[0]
            APs.append(V_last_four_paces[indexes[0]:indexes[-1]])
            AP_times.append(times[indexes[0]:indexes[-1]])
            bs_function = sp.interpolate.CubicSpline(AP_times[i_beat] - (times_of_beats[i_beat]), APs[i_beat])
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
                if time_between_peaks  >= 50:
                    biomarkers['EAD'].append(time_between_peaks)
            #  find maximum voltage and its index
            maxVoltage = np.max(Trace)
            maxVoltageIndex = np.where(Trace == maxVoltage)[0][0]
            # find the minimum voltage and its index after the peak
            minVoltageAfterPeak = np.min(Trace[maxVoltageIndex:])
            minVoltageIndex = np.where(Trace == minVoltageAfterPeak)[0][0]
            biomarkers['APA'].append(maxVoltage-minVoltageAfterPeak)
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
        s.reset()
    #  end for over simulations
    for ax in axs:
        ax.set_xlabel('time,s')
        ax.set_ylabel(ylabels[axs.tolist().index(ax)])
    fig.suptitle(r'Changing $I_{Kr}$ and $I_{CaL}$ conductances independently', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
    figName = FigureFolderName + '/Kr_and_CaL_gains_ind_variation_sigma2_025'
    for ax in axs:
        ax.set_xlim((t_start_storing) / 1000, t_end / 1000)
    plt.savefig(figName + '_last_four_paces.png', dpi=300)
    for ax in axs:
        ax.set_xlim((t_end - 1000) / 1000, t_end / 1000)
    plt.savefig(figName + '_last_pace.png', dpi=300)
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
        #     get unique values of gain_kr and gain_ca from dataframe

        axs1[iAxs].scatter(store_gains['gain_kr'], biomarkers[biomarker], s=5, label=r'Gain on $I_{Kr}$')
        axs1[iAxs].scatter(store_gains['gain_ca'], biomarkers[biomarker], s=5, label=r'Gain on $I_{CaL}$')
        axs1[iAxs].set_xlabel(r'Expression gain values')
        axs1[iAxs].set_ylabel(biomarker)
        if iAxs == 2:
            axs1[iAxs].legend()
        iAxs += 1
    fig1.suptitle(r'Biomarkers of AP traces obtained from independent gain variation', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
    plt.savefig(figName + '_biomarkers.png', dpi=300)
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
    fig2.suptitle(r'Histogramms of biomarkers obtained from independent gain variation', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
    plt.savefig(figName + '_biomarker_hists.png', dpi=300)
    ########################################################################################################################
    # save the simulation results and biomarkers into pickle files
    if not os.path.exists(DataFolderName + '/' + folderName):
        os.makedirs(DataFolderName + '/' + folderName)
    #  store all simulations into a csv file
    all_simulations_df.to_csv(DataFolderName + '/' + folderName + '/all_simulated_data.csv')
    pickle.dump(biomarkers, open(DataFolderName + '/' + folderName + '/biomarkers.pkl', 'wb'))
    del all_simulations_df, fig1, fig2, axs1, axs2, biomarkers
    gc.collect()