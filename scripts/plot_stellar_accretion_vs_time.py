#Routine to plot accretion rate onto binary and compare it
#to the accretion rate at a different point in the gas disk
#or to a different simulation set

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import sys
import os
from string import split


#rarely changed parameters
outer_rad = 70.0
outer_density = 1.0238548701E-4 # *fixed* density at outer radius
l = 1.0  #  power-law index of temperature profile
outer_accretion_coeff=1.07

ADD_PERICENTER_SEP=0
PLOT_TOTAL_ACCRETION=0

if __name__ == '__main__':


    filename=sys.argv[1]
    if (not "binary-accretion-rate" in filename):
        print "ERROR: data file of the wrong type"
        exit()
    if ("_q" in filename): qb= float(split(split(filename,"_q")[1],"_")[0])
    if ("_e" in filename): eb= float(split(split(filename,"_e")[1],"_")[0])
    if ("_h" in filename):h= float(split(split(filename,"_h")[1],"_")[0])
    if ("_alpha" in filename):alpha= float(split(split(filename,"_alpha")[1],"_")[0])
    if ("_eta" in filename):eta= float(split(split(filename,"_eta")[1],".txt")[0])
    if ("inbound" in filename): boundary_label = 'inbound_'
    else: boundary_label = ''
    if ("sink" in filename): sink_label='-sink'
    else: sink_label=''

        

    Mdot0 = -3 * np.pi * outer_density * alpha * h**2 * outer_rad**(1.5-l)
    #outer_accretion_coeff = 1.0
    Mdot0 *= outer_accretion_coeff
    print Mdot0, outer_accretion_coeff

    times =  np.loadtxt(filename)[:,0]/(2*np.pi)
    rates = np.loadtxt(filename)[:,1:]/np.abs(Mdot0)

    ind  = times >= 3250
    print rates.shape
    print rates[ind,0].mean(),rates[ind,1].mean(),(rates[ind,0]+rates[ind,1]).mean()


    Mdot1_mean = rates[times >= 2000,0].mean()
    Mdot2_mean = rates[times >= 2000,1].mean()
    
    

    Npanels=1
    if (len(sys.argv) > 2):
        xlim0 = float(sys.argv[2])
        xlim1 = float(sys.argv[3])
        
        
        if (len(sys.argv) > 4):
            Npanels = int(sys.argv[4])
        else:
            Npanels = 1
    else:
        xlim0 = times[0]
        xlim1 = times[-1]

        
    #prepare figure
    if (Npanels > 1):
        fig=plt.figure(1,figsize=(12,3.8*Npanels))
        fig.subplots_adjust(right=0.97,top=0.99,left=0.11,bottom=0.05*5/Npanels,hspace=0.1)
    else:
        fig=plt.figure(1,figsize=(12,3.4))
        fig.subplots_adjust(right=0.97,top=0.97,left=0.075,bottom=0.17,hspace=0.05)
        ax1 = fig.add_subplot(111)
    cmap=cm.jet


    if (ADD_PERICENTER_SEP):
        fig.subplots_adjust(right=0.93)
        import orbital

        

    
    for kk in range(Npanels):
        ax1 = fig.add_subplot(Npanels,1,kk+1)
        
        ind = ((times < (xlim1-xlim0)*(kk+1)*1.0/Npanels + xlim0) & 
               (times >= (xlim1-xlim0)*kk*1.0/Npanels + xlim0))
        binary_mdot1 = rates[ind,0]
        binary_mdot2 = rates[ind,1]
        
        if not (PLOT_TOTAL_ACCRETION):
            ax1.plot(times[ind],binary_mdot1,color='b',lw=1.4)
            ax1.plot(times[ind],binary_mdot2,color='crimson',lw=1.8,alpha=0.7)
        else:
            ax1.plot(times[ind],binary_mdot1+binary_mdot2,color='k',lw=2.0)


        #compute local maxima
        indlocalmax = (np.diff(np.sign(np.diff(binary_mdot1))) < 0).nonzero()[0] +1
        localmax = binary_mdot1[indlocalmax]
        tlocalmax = times[ind][indlocalmax]
        tlocalmax = tlocalmax[localmax > 2.0]
        localmax = localmax[localmax > 2.0]

        #compute time of pericenter
        peritime = np.arange(times[ind][0]+0.5,times[ind][-1]+0.5,1.0)
        for jj,peak in enumerate(localmax):
            peritime = peritime[np.abs(peritime-tlocalmax[jj]) == np.abs(peritime-tlocalmax[jj]).min()]
            #print tlocalmax[jj] - peritime
            

        #ax1.plot(tlocalmax,localmax,'ko',ms=7.0)
        
        if (ADD_PERICENTER_SEP):
            #ecc_anom = np.vectorize(orbital.KeplerEquation)(np.linspace(times[ind][0],times[ind][-1],300) + np.pi,eb)
            ecc_anom = np.vectorize(orbital.KeplerEquation)(times[ind]*2*np.pi + np.pi,eb)
            sep = (1.0 - eb*np.cos(ecc_anom))
            #find minima in separation
            indlocalmin = (np.diff(np.sign(np.diff(sep))) > 0).nonzero()[0] +1
            seplocalmin = sep[indlocalmin]
            #print seplocalmin
            tlocalmin = times[ind][indlocalmin]
            print tlocalmin

            ax2 = ax1.twinx()
            ax2.plot(times[ind],sep,color='purple',ls='--')
            ax2.set_xlim(xlim0,xlim1)
            ax2.set_ylim(1-eb,1+eb)
            ax2.set_ylabel(r'$|\mathbf{r}_2-\mathrm{r}_1|/a_b$',size=20)
            #ax2.plot(tlocalmin,seplocalmin,'ro',ms=7.0)
        
        [tick.label.set_fontsize(14) for tick in ax1.xaxis.get_major_ticks()]
        [tick.label.set_fontsize(16) for tick in ax1.yaxis.get_major_ticks()]
        
        ax1.set_xlim((xlim1-xlim0)*kk*1.0/Npanels + xlim0,(xlim1-xlim0)*(kk+1)*1.0/Npanels + xlim0)
        ax1.set_ylim(0.2,1.4)
        #ax1.set_ylim(0.0,5.0)
        yticks = ax1.get_yticks()
        ax1.set_yticks(yticks[1:-1])
        
        if (kk==20):
            ax1.text(0.01,0.85,r'$\langle\dot{M}_1\rangle/|\dot{M}_\mathrm{out}|=%.2f$   ' % Mdot1_mean,
                     color='blue',size=22,transform=ax1.transAxes)
            ax1.text(0.3,0.85,r'$\langle\dot{M}_2\rangle/|\dot{M}_\mathrm{out}|=%.2f$   ' % Mdot2_mean,
                     color='crimson',size=22,transform=ax1.transAxes)
            ax1.text(0.8,.87,r'$e_b=%.1f$' % eb,size=34,transform=ax1.transAxes)

    #ax1.legend(loc='upper right',prop={"size":12},ncol=3,columnspacing=0.8,handletextpad=0.1)

    #ax1.set_ylim(0.5,2.3)
    #ax1.set_xlim(2203,2230)
    

    ax1.set_xlabel(r"$t/P_b$",size=24,labelpad=-2)

    if (Npanels > 1): labelpad = 42
    else: labelpad = 9

    if (Npanels == 1):
        ax1.set_ylabel(r"$\dot{M}/|\dot{M}_\mathrm{out}|$",size=20,labelpad=labelpad)
    else:
        ax0 = fig.add_subplot(111,frameon=False)
        ax0.set_xticks([])
        ax0.set_yticks([])
        ax0.set_ylabel(r"$\dot{M}/|\dot{M}_\mathrm{out}|$",size=34,labelpad=labelpad)


    figlabels=''
    if (PLOT_TOTAL_ACCRETION): figlabels+='-total'
    if (ADD_PERICENTER_SEP): figlabels+='-separation'

    figname="figure-accretion-rate-evolution_binary"+figlabels+"_norbits%i-%i_q%3.1f_e%3.1f_h%3.1f_alpha%3.1f_eta%.2f.pdf"\
            % (xlim0,xlim1,qb,eb,h,alpha,eta)
    
    print "Saving figure:", figname
    fig.savefig(figname)
