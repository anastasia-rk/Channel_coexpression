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
import probscale
matplotlib.use('Agg') # to turn off interractive mode


if __name__ == '__main__':
    # if there is no folder for figures, create one
    FigureFolderName = 'Figures_Tomek_correct'
    DataFolderName = 'Simulated_data_Tomek_correct'
    expressions_to_vary = ['IKr', 'ICaL']
    gainNames = ['gain_kr', 'gain_ca']
    folderNames = []
    folderNames.append(expressions_to_vary[0] + '_' + expressions_to_vary[1] + '_independent')
    folderNames.append(expressions_to_vary[0] + '_' + expressions_to_vary[1] + '_cotranscripted_test')
    folderNames.append(expressions_to_vary[0] + '_' + expressions_to_vary[1] + '_cotranslated_test')
    folderNames.append(expressions_to_vary[0] + '_' + expressions_to_vary[1] + '_dependent_test')
    names = ['independent', 'cotranscripted', 'cotranslated', 'dependent']
    # only plot the
    # create a figure for plotting APD90 only for each of the folders
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    axs = axs.ravel()
    for iScenario, folderName in enumerate(folderNames):
        with open(os.path.join(DataFolderName + '/' + folderName +  '/biomarkers.pkl'), 'rb') as handle:
            biomarkers = pickle.load(handle)
        with open(os.path.join(DataFolderName + '/' + folderName +  '/simulationResults.pkl'), 'rb') as handle:
            simulationResults = pickle.load(handle)
        #################################################################################################################
        figName = FigureFolderName + '/' + folderName + '/'
        # plot the biomarkers
        # create a figure and make subplots for each biomarker as scutter plots
        fig1, axs1 = plt.subplots(2, 4, figsize=(15, 10))
        axs1 = axs1.ravel()
        # plot the biomarkers
        iAxs = 0
        for iBiomarker, biomarker in enumerate(biomarkers.keys()):
            if biomarker == 'Alternan' or biomarker == 'EAD':
                continue
            axs1[iAxs].scatter(simulationResults['gain_kr'], biomarkers[biomarker], s=5,
                               label=r'Multiplier on $g_{Kr}$')
            axs1[iAxs].scatter(simulationResults['gain_ca'], biomarkers[biomarker], s=5,
                               label=r'Multiplier on $g_{CaL}$')
            axs1[iAxs].set_xlabel(r'Expression gain values')
            axs1[iAxs].set_ylabel(biomarker)
            axs1[iAxs].set_xscale('symlog', base=2, linthresh=0.125)
            axs1[iAxs].set_xticks([0.25, 0.5, 1, 2, 4, 8])
            if iAxs == 2:
                axs1[iAxs].legend()
            iAxs += 1
        fig1.suptitle(r'Biomarkers of AP traces obtained from independent gain variation', fontsize=14)
        plt.tight_layout(rect=[0, 0.03, 0.98, 0.95])
        plt.savefig(figName + 'biomarkers.png', dpi=300)
        #################################################################################################################
        # plot histograms of the biomarkers
        fig2, axs2 = plt.subplots(2, 4, figsize=(15, 10))
        axs2 = axs2.ravel()
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
        plt.savefig(figName + 'biomarker_hists.png', dpi=300)

        fig3, axs3 = plt.subplots(2, 2, figsize=(10, 8))
        axs3 = axs3.ravel()
        axs3[0].plot([0.125, 8],[0.125, 8], color='k', linestyle='dashed', linewidth=0.5,label='$\lambda_1 = \lambda_2$')
        axs3[0].scatter(simulationResults['gain_kr'], simulationResults['gain_ca'], s=5,label='($\lambda_1,\lambda_2$) used in simulatiion')
        axs3[0].set_xlabel(r'Multiplier on $g_{Kr}$')
        axs3[0].set_ylabel(r'Multiplier on $g_{CaL}$')
        axs3[0].set_xscale('symlog', base=2, linthresh=0.125)
        axs3[0].set_yscale('symlog', base=2, linthresh=0.125)
        axs3[0].set_xticks([0.125, 0.25, 0.5, 1, 2, 4, 8])
        axs3[0].set_yticks([0.125, 0.25, 0.5, 1, 2, 4, 8])
        # fig3[0].suptitle(r'Conductance levels for $g_{Kr}$ and $g_{CaL}$', fontsize=14)
        axs3[0].set_xlim([0.125, 8])
        axs3[0].set_ylim([0.125, 8])
        axs3[0].legend(loc='upper left')
        # plot the ap traces
        for iTrace in range(1000):
            times = simulationResults['Time'][iTrace]/1000
            axs3[1].plot(times, simulationResults['Voltage'][iTrace],linewidth=0.5)
            axs3[2].plot(times, simulationResults['IKr'][iTrace],linewidth=0.5)
            axs3[3].plot(times, simulationResults['ICaL'][iTrace],linewidth=0.5)
        axs3[1].set_xlabel('Time [s]')
        axs3[1].set_ylabel('Voltage [mV]')
        axs3[2].set_xlabel('Time [s]')
        axs3[2].set_ylabel('$I_Kr$[nA]')
        axs3[3].set_xlabel('Time [s]')
        axs3[3].set_ylabel('$I_{CaL}$[nA]')
        axs3[1].set_xlim([453.0, 454.0])
        axs3[2].set_xlim([453.0, 454.0])
        axs3[3].set_xlim([453.0, 454.0])
        plt.tight_layout()
        plt.savefig(figName + 'conductance_levels.png', dpi=300)
        #################################################################################################################
        # plot the APD90 hists in the axs
        axs[iScenario].hist(biomarkers['APD90'], bins=20, density=True, histtype='step')
        # add mean, median and std to the plot
        axs[iScenario].axvline(np.mean(biomarkers['APD90']), color='k', linestyle='solid', linewidth=1)
        axs[iScenario].axvline(np.median(biomarkers['APD90']), color='k', linestyle='dashed', linewidth=1)
        axs[iScenario].axvline(np.median(biomarkers['APD90']) + np.std(biomarkers['APD90']), color='grey',
                              linestyle='dotted', linewidth=1)
        axs[iScenario].axvline(np.mean(biomarkers['APD90']) - np.std(biomarkers['APD90']), color='grey',
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
    fig.savefig(FigureFolderName + '/APD90_hists.png', dpi=400)
