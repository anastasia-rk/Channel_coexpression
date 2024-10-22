import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import myokit as mk
import os
# import csv
import gc
import pickle
from tqdm import tqdm
import probscale
matplotlib.use('Agg') # to turn off interractive mode


if __name__ == '__main__':
    # if there is no folder for figures, create one
    FigureFolderName = 'Figures_Tomek_correct'
    DataFolderName = 'Simulated_data_Tomek_correct'
    expressions_to_vary = ['IKr', 'ICaL']
    gains_to_vary = ['gain_kr', 'gain_ca']
    folderNamePrefix = '_'.join(expressions_to_vary)

    # possible currents that we can vary
    biomarkerNames = ['APD90', 'APD60', 'APD40', 'APA', 'TRI60', 'TRI40', 'EAD', '90PERCENT', '90to90', 'Alternan',
                      'Pro-Arrhythmic']
    currentNames = ['IKr', 'IKs', 'Ito', 'INa', 'ICaL']
    gainNames = ['gain_kr', 'gain_ks', 'gain_to', 'gain_na', 'gain_ca']
    ylabels = [r'$I_{Kr}$, ', r'$I_{Ks}$, A/F', r'$I_{to}$, A/F', r'$I_{Na}$, A/F', r'$I_{CaL}$, A/F']
    gain_labels  = [r'$g_{Kr}$', r'$g_{Ks}$', r'$g_{to}$', r'$g_{Na}$', r'$g_{CaL}$']

    # find_indeces of expressions to vary
    iExpressions = [currentNames.index(expression) for expression in expressions_to_vary]
    gain_labels = [gain_labels[i] for i in iExpressions]
    ylabels = [ylabels[i] for i in iExpressions]

    folderNames = []
    if folderNamePrefix == 'all':
        # in case we vary all - we don't yet have the translation implemented
        names = ['independent', 'cotranscripted', 'dependent']
    else:
        names = ['independent', 'cotranscripted', 'cotranslated', 'dependent']
    for name in names:
        forlderName = folderNamePrefix + '_' + name
        folderNames.append(forlderName)
    #################################################################################################################
    # create a figure for plotting APD90 only for each of the folders
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    axs = axs.ravel()
    for iScenario, folderName in enumerate(folderNames):
        figName = FigureFolderName + '/' + folderName + '/'
        with open(os.path.join(DataFolderName + '/' + folderName +  '/pro_arrhythmic_results.pkl'), 'rb') as handle:
            simRes_pro_arrhythmic = pickle.load(handle)
        with open(os.path.join(DataFolderName + '/' + folderName +  '/simulationResults.pkl'), 'rb') as handle:
            simulationResults = pickle.load(handle)
        biomarkers = pd.read_csv(DataFolderName + '/' + folderName +  '/biomarkers.csv')
        pro_arrythmic = biomarkers.loc[biomarkers['Pro-Arrhythmic'] == True]
        well_behaved = biomarkers.loc[biomarkers['Pro-Arrhythmic'] == False]

        # plot the scatter matrix of all columns that contain the word 'gain' split by pro-arrythmic and well-behaved on the same figure
        # extact the df that only has gains we varied
        well_behaved_varied_gains = well_behaved[gains_to_vary]
        pro_arrythmic_varied_gains = pro_arrythmic[gains_to_vary]
        scatter_matrix = pd.plotting.scatter_matrix(well_behaved_varied_gains,marker='o', color='k', alpha=0.2, figsize=(10, 10), diagonal='kde')
        # if the pro-arrythmic is not empty, plot it as well
        if not pro_arrythmic.empty:
            scatter_matrix = pd.plotting.scatter_matrix(pro_arrythmic_varied_gains, marker='x', color='orange', alpha=0.2, figsize=(10, 10), diagonal='kde')
        # plt.tight_layout()
        plt.savefig(figName + 'scatter_matrix.png', dpi=300)
        print(folderName + ' scatter matrix saved')
        #################################################################################################################

        # # plot the biomarkers
        # # create a figure and make subplots for each biomarker as scutter plots
        # fig1, axs1 = plt.subplots(2, 4, figsize=(15, 10))
        # axs1 = axs1.ravel()
        # # plot the biomarkers
        # iAxs = 0
        # for iBiomarker, biomarker in enumerate(biomarkerNames):
        #     if biomarker == 'Alternan' or biomarker == 'EAD' or biomarker == 'Pro-Arrhythmic':
        #         continue
        #     for iGain, gainName in enumerate(gains_to_vary):
        #         axs1[iAxs].scatter(well_behaved[gainName], well_behaved[biomarker], s=5,
        #                            label='Multiplier on ' + gain_labels[iGain])
        #         axs1[iAxs].scatter(pro_arrythmic[gainName], pro_arrythmic[biomarker], s=5, color='k')
        #         if iGain == len(gains_to_vary):
        #             axs1[iAxs].scatter(pro_arrythmic[gainName], pro_arrythmic[biomarker], s=5, color='k',
        #                                label='Pro-Arrhythmic')
        #     axs1[iAxs].set_xlabel(r'Multiplier  values')
        #     axs1[iAxs].set_ylabel(biomarker)
        #     axs1[iAxs].set_xscale('symlog', base=2, linthresh=0.125)
        #     axs1[iAxs].set_xticks([0.25, 0.5, 1, 2, 4, 8])
        #     if iAxs == 2:
        #         axs1[iAxs].legend()
        #     iAxs += 1
        # fig1.suptitle(r'Biomarkers of AP traces obtained from independent gain variation', fontsize=14)
        # plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
        # plt.savefig(figName + 'biomarkers.png', dpi=300)
        #################################################################################################################
        # # plot histograms of the biomarkers
        # fig2, axs2 = plt.subplots(2, 4, figsize=(15, 10))
        # axs2 = axs2.ravel()
        # # plot the biomarkers
        # iAxs = 0
        # for iBiomarker, biomarker in enumerate(biomarkerNames):
        #     if biomarker == 'Alternan' or biomarker == 'EAD' or biomarker == 'Pro-Arrhythmic':
        #         continue
        #     axs2[iAxs].hist(biomarkers[biomarker], bins=20, density=True, histtype='step')
        #     # add mean, median and std to the plot
        #     axs2[iAxs].axvline(np.mean(biomarkers[biomarker]), color='k', linestyle='solid', linewidth=1)
        #     axs2[iAxs].axvline(np.median(biomarkers[biomarker]), color='k', linestyle='dashed', linewidth=1)
        #     axs2[iAxs].axvline(np.median(biomarkers[biomarker]) + np.std(biomarkers[biomarker]), color='grey',
        #                        linestyle='dotted', linewidth=1)
        #     axs2[iAxs].axvline(np.mean(biomarkers[biomarker]) - np.std(biomarkers[biomarker]), color='grey',
        #                        linestyle='dotted', linewidth=1)
        #     ymin, ymax = axs2[iAxs].get_ylim()
        #     axs2[iAxs].text(np.mean(biomarkers[biomarker]), ymax * 0.9,
        #                     'mean: {:.2f}'.format(np.mean(biomarkers[biomarker])))
        #     axs2[iAxs].text(np.median(biomarkers[biomarker]), ymax * 0.8,
        #                     'md: {:.2f}'.format(np.median(biomarkers[biomarker])))
        #     axs2[iAxs].text(np.median(biomarkers[biomarker]) + np.std(biomarkers[biomarker]), ymax * 0.7,
        #                     'std: {:.2f}'.format(np.std(biomarkers[biomarker])))
        #     axs2[iAxs].set_xlabel(biomarker)
        #     iAxs += 1
        # fig2.suptitle(r'Histogramms of biomarkers obtained from independent gain variation', fontsize=14)
        # plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
        # plt.savefig(figName + 'biomarker_hists.png', dpi=300)
        ################################################################################################################
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # this plot currently only works for pair-wise assoctiation - we have a 2d plot of conductance multipliers
        # fig3, axs3 = plt.subplots(1, 2, figsize=(10, 4))
        # axs3 = axs3.ravel()
        # axs3[0].plot([0.125, 8],[0.125, 8], color='k', linestyle='dashed', linewidth=0.5,label='$\lambda_1 = \lambda_2$')
        # axs3[0].scatter(biomarkers[gains_to_vary[0]], biomarkers[gains_to_vary[1]], s=5,label='($\lambda_1,\lambda_2$) used in simulatiion')
        # axs3[0].scatter(biomarkers[gains_to_vary[0]], biomarkers[gains_to_vary[1]], s=5, color='k',label='Pro-Arrhythmic')
        # # axs3[0].scatter(simulationResults['gain_kr'], simulationResults['gain_ca'], s=5,label='($\lambda_1,\lambda_2$) used in simulatiion')
        # axs3[0].set_xlabel(r'Multiplier on '+ gains_to_vary[0])
        # axs3[0].set_ylabel(r'Multiplier on '+ gains_to_vary[1])
        # axs3[0].set_xscale('symlog', base=2, linthresh=0.125)
        # axs3[0].set_yscale('symlog', base=2, linthresh=0.125)
        # axs3[0].set_xticks([0.125, 0.25, 0.5, 1, 2, 4, 8])
        # axs3[0].set_yticks([0.125, 0.25, 0.5, 1, 2, 4, 8])
        # axs3[0].set_xlim([0.125, 8])
        # axs3[0].set_ylim([0.125, 8])
        # axs3[0].legend(loc='upper left')
        # # plot the ap traces
        # for iTrace in range(1000):
        #     times = simulationResults['Time'][iTrace]/1000
        #     axs3[1].plot(times, simulationResults['Voltage'][iTrace],linewidth=0.5)
        #     # axs3[2].plot(times, simulationResults[expressions_to_vary[0]][iTrace],linewidth=0.5)
        #     # axs3[3].plot(times, simulationResults[expressions_to_vary[1]][iTrace],linewidth=0.5)
        # axs3[1].set_xlabel('Time [s]')
        # axs3[1].set_ylabel('Voltage [mV]')
        # axs3[1].set_xlim([453.0, 454.0])
        # # axs3[2].set_xlabel('Time [s]')
        # # axs3[2].set_ylabel(ylabels[0])
        # # axs3[2].set_xlim([453.0, 454.0])
        # # axs3[3].set_xlabel('Time [s]')
        # # axs3[3].set_ylabel(ylabels[1])
        # # axs3[3].set_xlim([453.0, 454.0])
        # plt.tight_layout()
        # plt.savefig(figName + 'conductance_levels.png', dpi=300)
        #################################################################################################################
        # plot the APD90 hists in the axs
        axs[iScenario].hist(biomarkers['APD90'], bins=20, density=True, histtype='step',color='k')
        # add mean, median and std to the plot
        axs[iScenario].axvline(np.mean(biomarkers['APD90']), color='orange', linestyle='solid', linewidth=1)
        axs[iScenario].axvline(np.median(biomarkers['APD90']), color='orange', linestyle='dashed', linewidth=1)
        ymin, ymax = axs[iScenario].get_ylim()
        axs[iScenario].axvline(np.median(biomarkers['APD90']) + np.std(biomarkers['APD90']), color='orange',
                              linestyle='dotted', linewidth=1)
        axs[iScenario].axvline(np.mean(biomarkers['APD90']) - np.std(biomarkers['APD90']), color='orange',
                                linestyle='dotted', linewidth=1)
        axs[iScenario].text(np.mean(biomarkers['APD90']), ymax * 0.9,
                        'mean: {:.2f}'.format(np.mean(biomarkers['APD90'])))
        # axs[iScenario].text(np.median(biomarkers['APD90']), ymax * 0.8,
        #                 'md: {:.2f}'.format(np.median(biomarkers['APD90'])))
        axs[iScenario].text(np.median(biomarkers['APD90']) + np.std(biomarkers['APD90']), ymax * 0.7,
                        'std: {:.2f}'.format(np.std(biomarkers['APD90'])))
        axs[iScenario].set_xlabel('APD90')
        axs[iScenario].set_title(names[iScenario])
        axs[iScenario].set_xlim([0, 600])
    fig.tight_layout()
    fig.savefig(FigureFolderName + '/' + folderNamePrefix + '_APD90_hists.png', dpi=400)
