import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib as mpl

sys.path.append('../../../py')
from STREAM import STREAMOutput

FONTSIZE = 20
mpl.rcParams.update({'text.usetex': True, 'font.family': 'sans', 'font.size': FONTSIZE})

def drawPeakAndLoss():
    """
    Draw a plot with a output from the given STREAMOutput object.
    """
    factor = np.genfromtxt('DataUncertaintySweep/uncertaintyfactor.csv', delimiter=',').T

    peak = np.genfromtxt('DataUncertaintySweep/peakvalueNre.csv', delimiter=',').T

    loss = np.genfromtxt('DataUncertaintySweep/RunawayLossWithAvalanche.csv', delimiter=',').T


    fig, axs = plt.subplots(1, 2, figsize=(9,4))

    plotInternal(axs[0], factor, peak, ylabel=r'Peak $n_{\rm RE}$ value', color='tab:blue', showlabel=False, label=r' ')
    axs[0].set_ylim([0, 2e13])
    axs[0].set_xticks([0.001, 0.01, 0.1, 1])
    axs[0].set_yticks([0e13, 1e13, 2e13])
    plotInternal(axs[1], factor, loss, ylabel=r'Runaway loss fraction', color='tab:blue', showlabel=False, label=r' ')
    axs[1].set_xticks([0.001, 0.01, 0.1, 1])
    axs[1].set_ylim([0, 1.05])
    for i in range(axs.shape[0]):
        axs[i].grid(True, color='gainsboro', linewidth=0.5)
        axs[i].tick_params(axis='both', labelsize=FONTSIZE)

    fig.tight_layout()
    #fig.savefig('../runawaycase_REparam.pdf')

    plt.show()

def plotInternal(ax, x, y, ylabel, ylbl=True, xlbl=True, ylim=None, log=False, showlabel=False, label=None, *args, **kwargs):
    if label is not None and showlabel == False:
        label = None

    if log:
        ax.semilogy(x, y, label=label, *args, **kwargs)
    else:
        ax.semilogx(x, y, label=label, *args, **kwargs)

    if xlbl:
        ax.set_xlabel(r'Runaway loss factor $f_{\rm u}$')
    if ylbl:
        ax.set_ylabel(ylabel)

    if ylim is not None:
        ax.set_ylim(ylim)

def drawEoverED():
    legend = ['0.001', '0.01', '0.1', '1.0']
    for ext in [0.001, 0.01, 0.1, 1.0]:
        so1 = STREAMOutput(f'SweepSettingsEED/output1_{ext}.h5')
        EoverED1 = so1.eqsys.E_field.norm('ED')[:, 0] * 100
        t1 = so1.grid.t[:]
        dt1 = (t1[-1] - t1[0]) / (len(t1) - 1)
        dreicer1 = so1.other.fluid.gammaDreicer[:, 0]
        dreicer_int1 = np.sum(dreicer1*dt1)
        so2 = STREAMOutput(f'SweepSettingsEED/output2_{ext}.h5')
        EoverED2 = so2.eqsys.E_field.norm('ED')[:, 0] * 100
        t2 = so2.grid.t[:] + 1e-4
        dt2 = (t2[-1] - t2[0]) / (len(t2) - 1)
        dreicer2 = so2.other.fluid.gammaDreicer[:, 0]
        dreicer_int2 = np.sum(dreicer2*dt2)

        EoverED = np.append(EoverED1, EoverED2)
        t = np.append(t1, t2)

        print('factor=' + str(ext) + ', int dreicer=' + str(np.round_((dreicer_int1+dreicer_int2)/1e13,6)) + ', n_re=' + str(np.round_(so2.eqsys.n_re[-1][0]/1e13,6)) + ', loss fraction=' + str(np.round_(1-so2.eqsys.n_re[-1][0]/(dreicer_int1+dreicer_int2), 6)))

        plt.loglog(t, EoverED)
    plt.legend(legend)
    plt.show()

drawPeakAndLoss()