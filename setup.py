import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
# import pandas as pd
import myokit  as mk

# main
if __name__ == '__main__':
    model, protocol, embedded_script = mk.load('Tomek_cellml/Tomek_mid.mmt')
    s = mk.Simulation(model, protocol)
    d = s.run(1000)
    # def_state = s.default_state()  # Returns the original state
    # end_state = s.state()  # Returns the state at t=1000
    # t_end = s.time()
    filename = 'Tomek_cellml/test_output.csv'
    d.save_csv(filename, precision=64, order=None, delimiter=',', header=True)
    s.reset() # resets the state of the simulation

    # Plot simluation output
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    axs = axs.flatten()
    axs[0].plot(d.time(), d['membrane.v'])
    axs[1].plot(d.time(), d['ICaL.ICaL'])
    axs[2].plot(d.time(), d['IKr.IKr'])
    axs[3].plot(d.time(), d['IKs.IKs'])
    axs[4].plot(d.time(), d['INa.INa'])
    axs[5].plot(d.time(), d['Ito.Ito'])
    axs[0].set_ylabel(r'$V$')
    axs[1].set_ylabel(r'$I_{CaL}$')
    axs[2].set_ylabel(r'$I_{Kr}$')
    axs[3].set_ylabel(r'$I_{Ks}$')
    axs[4].set_ylabel(r'$I_{Na}$')
    axs[5].set_ylabel(r'$I_{to}$')
    for ax in axs:
        ax.set_xlabel('time,ms')
    fig.suptitle('Tomek mid model - baseline conductances',fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig('Figures/baseline_model.png',dpi=300)
    plt.show()

    ## try altering conductances in the model
    model, protocol, embedded_script = mk.load('Tomek_cellml/Tomek_mid_changing_conductances.mmt')
    # baseline conductances of the currents of interest
    conductance_kr_baseline = 0.0321
    conductance_ks_baseline = 0.0011
    conductance_to_baseline = 0.16
    conductance_na_baseline = 11.7802
    t_end = 10000
    # try changing conductances and pass them into the file
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    axs = axs.flatten()
    for gain_value in np.arange(0.5, 1.5, 0.1):
        # new conductances
        conductance_kr = conductance_kr_baseline * gain_value
        conductance_ks = conductance_ks_baseline * gain_value
        conductance_to = conductance_to_baseline * gain_value
        conductance_na = conductance_na_baseline * gain_value
        # placeholder for ICal conductance
        # run the simulation
        s = mk.Simulation(model, protocol)
        # set constants to new values
        s.set_constant('IKr.GKr_b', conductance_kr)
        s.set_constant('IKs.GKs_b', conductance_ks)
        s.set_constant('Ito.Gto_b', conductance_to)
        s.set_constant('INa.GNa', conductance_na)
        d = s.run(t_end)
        axs[0].plot(d.time(), d['membrane.v'])
        axs[1].plot(d.time(), d['ICaL.ICaL'])
        axs[2].plot(d.time(), d['IKr.IKr'],label='gain = %.3f' % gain_value)
        axs[3].plot(d.time(), d['IKs.IKs'])
        axs[4].plot(d.time(), d['INa.INa'])
        axs[5].plot(d.time(), d['Ito.Ito'])
        s.reset()
    axs[0].set_ylabel(r'$V$')
    axs[1].set_ylabel(r'$I_{CaL}$')
    axs[2].set_ylabel(r'$I_{Kr}$')
    axs[3].set_ylabel(r'$I_{Ks}$')
    axs[4].set_ylabel(r'$I_{Na}$')
    axs[5].set_ylabel(r'$I_{to}$')
    for ax in axs:
        ax.set_xlabel('time,ms')
    axs[2].legend(bbox_to_anchor=(0.9, 0.2, 0.2, 1), loc="lower left",
                borderaxespad=0,ncol=1)
    fig.suptitle('Simultaneous change in conductances',fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    figName = 'Figures/Simultaneous_change_in_conductances'
    plt.savefig(figName + '.png', dpi=300)
    for ax in axs:
        ax.set_xlim(t_end - 1000, t_end)
    plt.savefig(figName + '_last_pace.png', dpi=300)
    plt.show()
    print('pause here')


    # replicating Ballouz et al. 2021 - kinda, first scenario only
    nTries = 50
    # gain factors
    gain_kr, gain_ks, gain_to, gain_ks, gain_na = 1, 1, 1, 1, 1
    sigma_kr = np.sqrt(0.25)
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    axs = axs.flatten()
    for iTry in range(nTries):
        # sample a number from lognormal distribution
        gain_kr = np.random.lognormal(mean=0.0, sigma=sigma_kr, size=1)
        #  and then multiply the baseline conductances by the gain
        conductance_kr = conductance_kr_baseline * gain_kr
        conductance_ks = conductance_ks_baseline * gain_ks
        conductance_to = conductance_to_baseline * gain_to
        conductance_na = conductance_na_baseline * gain_na
        # plcaeholder for ICal conductance
        # run the simulation
        s = mk.Simulation(model, protocol)
        # set constants to new values
        s.set_constant('IKr.GKr_b', conductance_kr)
        s.set_constant('IKs.GKs_b', conductance_ks)
        s.set_constant('Ito.Gto_b', conductance_to)
        s.set_constant('INa.GNa', conductance_na)
        d = s.run(t_end)
        axs[0].plot(d.time(), d['membrane.v'])
        axs[1].plot(d.time(), d['ICaL.ICaL'])
        axs[2].plot(d.time(), d['IKr.IKr'])
        axs[3].plot(d.time(), d['IKs.IKs'])
        axs[4].plot(d.time(), d['INa.INa'])
        axs[5].plot(d.time(), d['Ito.Ito'])
        s.reset()
    axs[0].set_ylabel(r'$V$')
    axs[1].set_ylabel(r'$I_{CaL}$')
    axs[2].set_ylabel(r'$I_{Kr}$')
    axs[3].set_ylabel(r'$I_{Ks}$')
    axs[4].set_ylabel(r'$I_{Na}$')
    axs[5].set_ylabel(r'$I_{to}$')
    for ax in axs:
        ax.set_xlabel('time,ms')
    fig.suptitle(r'Changing $I_{Kr}$ conductances',fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    figName = 'Figures/G_Kr_gain__lognormal_sigma2_025'
    # for all ax, set xlim from 9000 to 10000
    plt.savefig(figName + '.png', dpi=300)
    for ax in axs:
        ax.set_xlim(t_end - 1000, t_end)
    plt.savefig(figName + '_last_pace.png', dpi=300)
    plt.show()

    # independently varying expression of I_Kr and I_Na
    # gain factors
    gain_kr, gain_ks, gain_to, gain_ks, gain_na = 1, 1, 1, 1, 1
    sigma_kr, sigma_na = np.sqrt(0.25), np.sqrt(0.25)
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    axs = axs.flatten()
    for iTry in range(nTries):
        # sample gains from lognormal distribution
        gain_kr = np.random.lognormal(mean=0.0, sigma=sigma_kr, size=1)
        gain_na = np.random.lognormal(mean=0.0, sigma=sigma_na, size=1)
        #  and then multiply the baseline conductances by the gain
        conductance_kr = conductance_kr_baseline * gain_kr
        conductance_ks = conductance_ks_baseline * gain_ks
        conductance_to = conductance_to_baseline * gain_to
        conductance_na = conductance_na_baseline * gain_na
        # plcaeholder for ICal conductance
        # run the simulation
        s = mk.Simulation(model, protocol)
        # set constants to new values
        s.set_constant('IKr.GKr_b', conductance_kr)
        s.set_constant('IKs.GKs_b', conductance_ks)
        s.set_constant('Ito.Gto_b', conductance_to)
        s.set_constant('INa.GNa', conductance_na)
        d = s.run(t_end)
        axs[0].plot(d.time(), d['membrane.v'])
        axs[1].plot(d.time(), d['ICaL.ICaL'])
        axs[2].plot(d.time(), d['IKr.IKr'])
        axs[3].plot(d.time(), d['IKs.IKs'])
        axs[4].plot(d.time(), d['INa.INa'])
        axs[5].plot(d.time(), d['Ito.Ito'])
        s.reset()
    axs[0].set_ylabel(r'$V$')
    axs[1].set_ylabel(r'$I_{CaL}$')
    axs[2].set_ylabel(r'$I_{Kr}$')
    axs[3].set_ylabel(r'$I_{Ks}$')
    axs[4].set_ylabel(r'$I_{Na}$')
    axs[5].set_ylabel(r'$I_{to}$')
    for ax in axs:
        ax.set_xlabel('time,ms')
    fig.suptitle(r'Changing $I_{Kr}$ and $I_{Na}$ conductances independently',fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    figName = 'Figures/Kr_and_Na_gains_indy_variation_sigma2_025'
    plt.savefig(figName+'.png',dpi=300)
    for ax in axs:
        ax.set_xlim(t_end-1000,t_end)
    plt.savefig(figName+'_last_pace.png', dpi=300)
    plt.show()

    # dependent co-expression of I_Kr and I_Na
    # gain factors
    gain_kr, gain_ks, gain_to, gain_ks, gain_na = 1, 1, 1, 1, 1
    sigma_kr, sigma_na = np.sqrt(0.25), np.sqrt(0.25)
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    axs = axs.flatten()
    for iTry in range(nTries):
        # sample a number from lognormal distribution and assing to all affected gains
        gain_kr = np.random.lognormal(mean=0.0, sigma=sigma_kr, size=1)
        gain_na = gain_kr.copy()
        #  and then multiply the baseline conductances by the gain
        conductance_kr = conductance_kr_baseline * gain_kr
        conductance_ks = conductance_ks_baseline * gain_ks
        conductance_to = conductance_to_baseline * gain_to
        conductance_na = conductance_na_baseline * gain_na
        # plcaeholder for ICal conductance
        # run the simulation
        s = mk.Simulation(model, protocol)
        # set constants to new values
        s.set_constant('IKr.GKr_b', conductance_kr)
        s.set_constant('IKs.GKs_b', conductance_ks)
        s.set_constant('Ito.Gto_b', conductance_to)
        s.set_constant('INa.GNa', conductance_na)
        d = s.run(t_end)
        axs[0].plot(d.time(), d['membrane.v'])
        axs[1].plot(d.time(), d['ICaL.ICaL'])
        axs[2].plot(d.time(), d['IKr.IKr'])
        axs[3].plot(d.time(), d['IKs.IKs'])
        axs[4].plot(d.time(), d['INa.INa'])
        axs[5].plot(d.time(), d['Ito.Ito'])
        s.reset()
    axs[0].set_ylabel(r'$V$')
    axs[1].set_ylabel(r'$I_{CaL}$')
    axs[2].set_ylabel(r'$I_{Kr}$')
    axs[3].set_ylabel(r'$I_{Ks}$')
    axs[4].set_ylabel(r'$I_{Na}$')
    axs[5].set_ylabel(r'$I_{to}$')
    for ax in axs:
        ax.set_xlabel('time,ms')
    fig.suptitle(r'Changing $I_{Kr}$ and $I_{Na}$ conductances independently',fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    figName = 'Figures/Kr_and_Na_gains_dependent_variation_sigma2_025'
    plt.savefig(figName+'.png',dpi=300)
    for ax in axs:
        ax.set_xlim(t_end-1000,t_end)
    plt.savefig(figName+'_last_pace.png', dpi=300)
    plt.show()

    # # independently varying expresssion of I_Kr and I_Ks
    # # gain factors
    # gain_kr, gain_ks, gain_to, gain_ks, gain_na = 1, 1, 1, 1, 1
    # sigma_kr, sigma_ks = np.sqrt(0.25), np.sqrt(0.25)
    # fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    # axs = axs.flatten()
    # for iTry in range(nTries):
    #     # sample a number from lognormal distribution
    #     gain_kr = np.random.lognormal(mean=0.0, sigma=sigma_kr, size=1)
    #     gain_ks = np.random.lognormal(mean=0.0, sigma=sigma_ks, size=1)
    #     #  and then multiply the baseline conductances by the gain
    #     conductance_kr = conductance_kr_baseline * gain_kr
    #     conductance_ks = conductance_ks_baseline * gain_ks
    #     conductance_to = conductance_to_baseline * gain_to
    #     conductance_na = conductance_na_baseline * gain_na
    #     # plcaeholder for ICal conductance
    #     # run the simulation
    #     s = mk.Simulation(model, protocol)
    #     # set constants to new values
    #     s.set_constant('IKr.GKr_b', conductance_kr)
    #     s.set_constant('IKs.GKs_b', conductance_ks)
    #     s.set_constant('Ito.Gto_b', conductance_to)
    #     s.set_constant('INa.GNa', conductance_na)
    #     d = s.run(t_end)
    #     axs[0].plot(d.time(), d['membrane.v'])
    #     axs[1].plot(d.time(), d['ICaL.ICaL'])
    #     axs[2].plot(d.time(), d['IKr.IKr'])
    #     axs[3].plot(d.time(), d['IKs.IKs'])
    #     axs[4].plot(d.time(), d['INa.INa'])
    #     axs[5].plot(d.time(), d['Ito.Ito'])
    #     s.reset()
    # axs[0].set_ylabel(r'$V$')
    # axs[1].set_ylabel(r'$I_{CaL}$')
    # axs[2].set_ylabel(r'$I_{Kr}$')
    # axs[3].set_ylabel(r'$I_{Ks}$')
    # axs[4].set_ylabel(r'$I_{Na}$')
    # axs[5].set_ylabel(r'$I_{to}$')
    # for ax in axs:
    #     ax.set_xlabel('time,ms')
    # fig.suptitle(r'Changing $I_{Kr}$ and $I_{Ks}$ conductances independently',fontsize=14)
    # plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    # figName = 'Figures/Kr_and_Ks_gains_indy_variation_sigma2_025'
    # plt.savefig(figName+'.png',dpi=300)
    # for ax in axs:
    #     ax.set_xlim(t_end-1000,t_end)
    # plt.savefig(figName+'_last_pace.png', dpi=300)
    # plt.show()
    print('pause here')





    # we want to sample amplifying gains from a lognormal distribution
    # placeholder for sampling
