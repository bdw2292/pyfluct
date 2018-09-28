import matplotlib.pyplot as plt
import os
import scipy as sp
import numpy as np

script_dir = os.path.dirname(__file__)
results_dir = os.path.join(script_dir, 'Plots/')
title_font = {'fontname':'Arial', 'size':'14', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'} # Bottom vertical alignment for more space
axis_font = {'fontname':'Arial', 'size':'14'}

def AutoDetectSteadyState(array):
    tol=.01
    mean=np.mean(array)
    for i in range(len(array)):
        element=array[i]
        if np.abs(mean-element)<=tol:
            return i
    
def SingleColorNoise(model): #Plot single color noise simulation from the deffinition of intrinsic,extrinsic and total noise
    for i in range(len(model.intnoise)):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        line1,=ax1.plot(model.t,model.intnoise[i],color='red',label='$\sigma^{2}_{int}$')
        line2,=ax1.plot(model.t,model.extnoise[i],color='blue',label='$\sigma^{2}_{ext}$')
        line3,=ax1.plot(model.t,model.totnoise[i],color='green',label='$\sigma^{2}_{tot}$')
        ax1.legend(handles=[line1,line2,line3],loc='best')
        if model.static==True and i==model.drawnrateindex:
            model.K[model.drawnrateindex]=model.mean
            plt.title('IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font)
        else:
            plt.title('IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font)
        if model.birthdist==None:
            plt.xlabel(r'Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font)
        else:
            plt.xlabel(r'Dis=%s, Obs=%s, Var=%s, ObsVar=%s, Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.birthdist,model.observations,model.var,model.obsvar,model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font)
        plt.ylabel('Single Color $\sigma^{2}$, Spec=%s'%(model.ICspec[i]),**axis_font)
        fig = plt.gcf()
        plt.show()
        fig.savefig(results_dir + model.reac + 'SingleColorNoiseSimulation' + model.ICspec[i] + '.png', format='png')
        
def TwoColorNoise(model): #Plot Two color noise simulation from derived formulas
    model.trials=1
    for i in range(len(model.intnoise)):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        line1,=ax1.plot(model.t,model.intnoise[i],color='red',label='$\sigma^{2}_{int}$')
        line2,=ax1.plot(model.t,model.extnoise[i],color='blue',label='$\sigma^{2}_{ext}$')
        line3,=ax1.plot(model.t,model.totnoise[i],color='green',label='$\sigma^{2}_{tot}$')
        ax1.legend(handles=[line1,line2,line3],loc='best')
        if model.static==True and i==model.drawnrateindex:
            model.K[model.drawnrateindex]=model.mean
            plt.title('IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font)
        else:
            plt.title('IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font)
        if model.birthdist==None:
            plt.xlabel(r'Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font)
        else:
            plt.xlabel(r'Dis=%s, Obs=%s, Var=%s, ObsVar=%s, Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.birthdist,model.observations,model.var,model.obsvar,model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font)
        fig = plt.gcf()
        plt.show()
        filename=model.reac +'TwoColorNoiseSimulation' + model.ICspec[i] + '.png'
        fig.savefig(results_dir+filename, format='png')

def TotalNoiseComparison(model):
    for i in range(len(model.totnoise)):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        line1,=ax1.plot(model.t,model.testTotnoise[i],color='blue',label='From Formula')
        line2,=ax1.plot(model.t,model.totnoise[i],color='red',label='From $\sigma^{2}_{int} , \sigma^{2}_{ext}$')
        ax1.legend(handles=[line1,line2],loc='best')
        if model.static==True and i==model.drawnrateindex:
            model.K[model.drawnrateindex]=model.mean
            plt.title('IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font)
        else:
            plt.title('IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font)
        if model.birthdist==None:
            plt.xlabel(r'Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font)
        else:
            plt.xlabel(r'Dis=%s, Obs=%s, Var=%s, ObsVar=%s, Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.birthdist,model.observations,model.var,model.obsvar,model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font)
        plt.ylabel('Comparison $\sigma^{2}_{tot}$, Spec=%s'%(model.ICspec[i]),**axis_font)
        fig = plt.gcf()
        plt.show()
        fig.savefig(results_dir + model.reac + 'TotalNoiseComparison' + model.ICspec[i] + '.png', format='png')
    
def ExtNoiseMethodComparison(model):
    string1='Dual Reporter'
    string2='Single Color'
    if model.color==True:
        label1=string1
        label2=string2
    else:
        label1=string2
        label2=string1
    for i in range(len(model.extnoise)):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        line1,=ax1.plot(model.t,model.extnoise[i],color='blue',label=label1)
        line2,=ax1.plot(model.t,model.benchmarkextrinsic[i],color='red',label=label2)
        ax1.legend(handles=[line1,line2],loc='best')
        if model.static==True and i==model.drawnrateindex:
            model.K[model.drawnrateindex]=model.mean
            plt.title('IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font)
        else:
            plt.title('IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font)
        if model.birthdist==None:
            plt.xlabel(r'Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font)
        else:
            plt.xlabel(r'Dis=%s, Obs=%s, Var=%s, ObsVar=%s, Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.birthdist,model.observations,model.var,model.obsvar,model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font)
        plt.ylabel('Dual-Reporter & Single-Color $\sigma^{2}_{ext}$ Spec=%s'%(model.ICspec[i]),**axis_font)
        fig = plt.gcf()
        plt.show()
        fig.savefig(results_dir + model.reac + 'ExtrinsicNoiseMethodComparison' + model.ICspec[i] + '.png', format='png')

def IntNoiseMethodComparison(model):
    string1='Dual Reporter'
    string2='Single Color'
    if model.color==True:
        label1=string1
        label2=string2
    else:
        label1=string2
        label2=string1
    for i in range(len(model.intnoise)):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        line1,=ax1.plot(model.t,model.intnoise[i],color='blue',label=label1)
        line2,=ax1.plot(model.t,model.benchmarkintrinsic[i],color='red',label=label2)
        ax1.legend(handles=[line1,line2],loc='best')
        if model.static==True and i==model.drawnrateindex:
            model.K[model.drawnrateindex]=model.mean
            plt.title('IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font)
        else:
            plt.title('IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font)
        if model.birthdist==None:
            plt.xlabel(r'Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font)
        else:
            plt.xlabel(r'Dis=%s, Obs=%s, Var=%s, ObsVar=%s, Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.birthdist,model.observations,model.var,model.obsvar,model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font)
        plt.ylabel('Dual-Reporter & Single-Color $\sigma^{2}_{int}$ Spec=%s'%(model.ICspec[i]),**axis_font)
        fig = plt.gcf()
        plt.show()
        fig.savefig(results_dir + model.reac + 'IntrinsicNoiseMethodComparison' + model.ICspec[i] + '.png', format='png')

def Benchmark(model):
    for i in range(len(model.extnoise)):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        line1,=ax1.plot(model.t,model.extnoise[i],color='blue',label='$\sigma^{2}_{ext}$')
        line2,=ax1.plot(model.t,model.benchmarkextrinsic[i],color='yellow',label='Extrinsic Benchmark')
        ax1.legend(handles=[line1,line2],loc='best')
        if model.static==True and i==model.drawnrateindex:
            model.K[model.drawnrateindex]=model.mean
            plt.title('IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font)
        else:
            plt.title('IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font)
        if model.birthdist==None:
            plt.xlabel(r'Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font)
        else:
            plt.xlabel(r'Dis=%s, Obs=%s, Var=%s, ObsVar=%s, Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.birthdist,model.observations,model.var,model.obsvar,model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font)
        plt.ylabel('Single Color $\sigma^{2}_{ext}$, Spec=%s'%(model.ICspec[i]),**axis_font)
        fig = plt.gcf()
        plt.show()
        fig.savefig(results_dir + model.reac + 'SingleColorNoiseVsBenchmark' + model.ICspec[i] + '.png', format='png')
        
        
def ExtrinsicAverageOfSpecices(model):
    limitspecies=[np.full((len(model.t),),model.extaveofintaveofp[l][-1]) for l in range(len(model.indexarray))]
    for i in range(len(limitspecies)):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        line1,=ax1.plot(model.t,model.extaveofintaveofp[i],color='blue',label='Data SS')
        line2,=ax1.plot(model.t,limitspecies[i],color='red',label='Data SS t->$\infty$')
        ax1.legend(handles=[line1,line2],loc='best')
        if model.static==True and i==model.drawnrateindex:
            model.K[model.drawnrateindex]=model.mean
            plt.title('IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font)
        else:
            plt.title('IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font)
        if model.birthdist==None:
            plt.xlabel(r'Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font)
        else:
            plt.xlabel(r'Dis=%s, Obs=%s, Var=%s, ObsVar=%s, Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.birthdist,model.observations,model.var,model.obsvar,model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font)
        plt.ylabel('Average Over Birth Of Species, Spec=%s'%(model.ICspec[i]),**axis_font)
        fig = plt.gcf()
        plt.show()
        fig.savefig(results_dir + model.reac + 'ExtAveOfSpecies' + model.ICspec[i] + '.png', format='png')
    
def NoiseRatioTest(model):
    extRQratio=np.divide(model.extnoise[1],model.extnoise[0])
    intRQratio=np.divide(model.intnoise[1],model.intnoise[0])
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    line1,=ax1.plot(model.t,extRQratio,color='blue',label=r'$\frac{\sigma^{2}_{extR}}{\sigma^{2}_{extQ}}$')
    line2,=ax1.plot(model.t,intRQratio,color='red',label=r'$\frac{\sigma^{2}_{intR}}{\sigma^{2}_{intQ}}$')
    ax1.legend(handles=[line1,line2],loc='best')
    model.K[model.drawnrateindex]=model.mean
    plt.title('IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font)
    if model.birthdist==None:
        plt.xlabel(r'Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font)
    else:
        plt.xlabel(r'Dis=%s, Obs=%s, Var=%s, ObsVar=%s, Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.birthdist,model.observations,model.var,model.obsvar,model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font)
    plt.ylabel(' Noise Ratio r_b/r_d=%s, r_d/q_d=%s'%(model.R_b/model.R_d,model.R_d/model.Q_d),**axis_font)
    fig = plt.gcf()
    plt.show()
    fig.savefig(results_dir + model.reac + 'NoiseRatioTest' + '.png', format='png')

def PredictedVsSimulatedIntNoise(model):
    for i in range(len(model.indexarray)):
        index=AutoDetectSteadyState(model.autocorrelation[i])  
        fig1 = plt.figure()
        lines=[]
        ax1 = fig1.add_subplot(211)
        line1,=ax1.plot(model.t[index:],model.intnoise[i][index:],color='red',label='Simulation')
        lines.append(line1)
        if i==1 and model.osc==False and model.static==True:
            difference=np.full((len(model.t[index:]),),model.autocorrelation[1][0]-model.autocorrelation[1][-1])
            line3,=ax1.plot(model.t[index:],difference,color='green',label='Cov[0]-Cov[inf]')
            lines.append(line3)
        if model.predintnoise[i]!=None:
            line2,=ax1.plot(model.t[index:],np.full((len(model.t[index:]),),model.predintnoise[i][index:]),color='blue',label='Prediction')
            ax2=fig1.add_subplot(212)        
            ax2.plot(model.t[index:],np.full((len(model.t[index:]),),np.abs(model.intnoise[i][index:]-model.predintnoise[i][index:])/model.predintnoise[i][index:]),'or')
            ax2.set_ylabel('Relative Residual')
            ax2.grid()
            lines.append(line2)
        if model.static==True and i==model.drawnrateindex:
            model.K[model.drawnrateindex]=model.mean
            fig1.text(.5,.95,'IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font,horizontalalignment='center')
        else:
            fig1.text(.5,.95,'IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font,horizontalalignment='center')
        ax1.legend(handles=lines,loc='best')
        if model.birthdist==None:
            fig1.text(0.5, 0.04, r'Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s $\tau$'%(model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font, ha='center', va='center')
        else:
            fig1.text(0.5, 0.04, r'Dis=%s, Obs=%s, Var=%s, ObsVar=%s, Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.birthdist,model.observations,model.var,model.obsvar,model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font, ha='center', va='center')
        ax1.set_ylabel('Single Color $\sigma^{2}_{int}$ Propagation, Spec=%s'%(model.ICspec[i]),**axis_font)
        fig = plt.gcf()
        plt.show()
        fig.savefig(results_dir + model.reac + 'PredictedVsSimulatedIntNoise' + model.ICspec[i] + '.png', format='png')
    
def PredictedVsSimulatedTotNoise(model):
    for i in range(len(model.indexarray)):
        index=AutoDetectSteadyState(model.autocorrelation[i])
        lines=[]
        lines2=[]
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(211)
        line1,=ax1.plot(model.t[index:],model.totnoise[i][index:],color='red',label='Simulated')
        lines.append(line1)
        if model.predtotnoise[i]!=None: 
            if i==1 and model.osc==False and model.static==True:
                autocortot=np.full((len(model.t[index:]),),model.autocorrelation[i][0])
                func=np.exp(-model.R_d*np.array(model.shift))
                intsimps=sp.integrate.simps(func*model.aveautocorrelation[0],x=model.shift)
                propagatedtotnoiseR=(model.R_b/model.R_d)*model.extaveofintaveofp[0][-1]+(model.R_b**2/model.R_d)*intsimps
                propagatedtotnoiseRvec=np.full((len(model.t[index:]),),propagatedtotnoiseR)
                line2,=ax1.plot(model.t[index:], np.full((len(model.t[index:]),),propagatedtotnoiseR),color='green',label='Noise Prop Formula')
                line3,=ax1.plot(model.t[index:],autocortot,color='blue',label='COV_RR(0) Sim')
                lines.append(line2)
                lines.append(line3)
                line4,=ax1.plot(model.t[index:],np.full((len(model.t[index:]),),model.predtotnoise[i][index:]),color='yellow',label='Predicted')
                lines.append(line4)
                ax2=fig1.add_subplot(212)          
                ax2.set_ylabel('Relative Residual')
                line5,=ax2.plot(model.t[index:],np.abs(model.predtotnoise[i][index:]-propagatedtotnoiseRvec)/model.predtotnoise[i][index:],'or',label='Noise Prop Formula')
                line6,=ax2.plot(model.t[index:],np.abs(model.predtotnoise[i][index:]-autocortot)/model.predtotnoise[i][index:],'og',label='COV_RR(0) Sim)')
                line7,=ax2.plot(model.t[index:],np.abs(model.predtotnoise[i][index:]-model.totnoise[i][index:])/model.predtotnoise[i][index:],'ob',label='Simulation')
                lines2.append(line5)
                lines2.append(line6)
                lines2.append(line7)
                ax2.grid()
                ax2.legend(handles=lines2,loc='best')
            else:
                line4,=ax1.plot(model.t[index:],np.full((len(model.t[index:]),),model.predtotnoise[i][index:]),color='yellow',label='Predicted')
                lines.append(line4)
                ax2=fig1.add_subplot(212)          
                ax2.set_ylabel('Relative Residual')
                line7,=ax2.plot(model.t[index:],np.full((len(model.t[index:]),),np.abs(model.predtotnoise[i][index:]-model.totnoise[i][index:])/model.predtotnoise[i][index:]),'ob',label='Simulation')
                lines2.append(line7)
                ax2.grid()
                ax2.legend(handles=lines2,loc='best')
        ax1.legend(handles=lines,loc='best')
        if model.static==True and i==model.drawnrateindex:
            model.K[model.drawnrateindex]=model.mean
            fig1.text(.5,.95,'IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font,horizontalalignment='center')
        else:
            fig1.text(.5,.95,'IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font,horizontalalignment='center')
        if model.birthdist==None:
            fig1.text(0.5, 0.04, r'Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s $\tau$'%(model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font, ha='center', va='center')
        else:
            fig1.text(0.5, 0.04, r'Dis=%s, Obs=%s, Var=%s, ObsVar=%s, Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.birthdist,model.observations,model.var,model.obsvar,model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font, ha='center', va='center')
        ax1.set_ylabel('Single Color $\sigma^{2}_{tot}$ Propagation, Spec=%s'%(model.ICspec[i]),**axis_font)
        fig = plt.gcf()
        plt.show()
        fig.savefig(results_dir + model.reac + 'PredictedVsSimulatedTotNoise' + model.ICspec[i] + '.png', format='png')
    
def AutoCorrelationVsTau(model):
    extnoiseave=[np.full((len(model.shift),),np.mean(model.extnoise[l])) for l in range(len(model.indexarray))]
    for i in range(len(model.indexarray)):
        index=AutoDetectSteadyState(model.autocorrelation[i])
        index=0
        fig1 = plt.figure()
        lines=[]
        ax1 = fig1.add_subplot(211)
        line1,=ax1.plot(model.shift[index:],model.autocorrelation[i][index:],color='blue',label='Cov t->$\infty$')
        lines.append(line1)
        if i==1 and model.osc==False and model.static==True:
            line2,=ax1.plot(model.shift[index:],np.full((len(model.shift[index:]),),model.predextnoise[i][-1]),color='red',label='$\sigma^{2}_{ext}$ Pred')
            line3,=ax1.plot(model.shift[index:],extnoiseave[i][index:],color='green',label='$\sigma^{2}_{ext}$ Sim')
            lines.append(line2)
            lines.append(line3)
        if model.predautocor[i]!=None:
            line4,=ax1.plot(model.shift[index:],model.predautocor[i],color='yellow',label='Pred Cov t->$\infty$')
            lines.append(line4)
            ax2=fig1.add_subplot(212)
            ax2.set_ylabel('Relative Residual')    
            ax2.plot(model.shift[index:],np.abs(model.predautocor[i]-model.autocorrelation[i][index:])/model.predautocor[i],'or')
            ax2.grid()
        ax1.legend(handles=lines,loc='best')
        if model.static==True and i==model.drawnrateindex:
            model.K[model.drawnrateindex]=model.mean
            fig1.text(.5,.95,'IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font,horizontalalignment='center')
        else:
            fig1.text(.5,.95,'IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font,horizontalalignment='center')
        if model.birthdist==None:
            fig1.text(0.5, 0.04, r'Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s $\tau$'%(model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font, ha='center', va='center')
        else:
            fig1.text(0.5, 0.04, r'Dis=%s, Obs=%s, Var=%s, ObsVar=%s, Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.birthdist,model.observations,model.var,model.obsvar,model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font, ha='center', va='center')
        ax1.set_ylabel('COV, Spec=%s'%(model.ICspec[i]),**axis_font)
        fig = plt.gcf()
        plt.show()
        fig.savefig(results_dir + model.reac + 'AutoCorrelationVsTime' + model.ICspec[i] + '.png', format='png')
    
def DerivativeOfAutoCorrelationVsLag(model):
    slope=np.gradient(model.autocorrelation[0],model.density,edge_order=2)
    for i in range(len(model.indexarray)):
        lines=[]
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        if i==1:
            cons=-(model.R_d**2+model.R_d*model.Q_d-model.R_b*model.Q_d)/(model.R_b+model.R_d+model.Q_d)
            line1,=ax1.plot(model.shift,slope*cons,color='blue',label=r'$\frac{d}{dt} COV*C, \tau$->0')
            lines.append(line1)
        line2,=ax1.plot(model.shift,slope,color='red',label=r'$\frac{d}{dt} COV, \tau$->0')
        lines.append(line2)
        ax1.legend(handles=lines,loc='best')
        if model.static==True and i==model.drawnrateindex:
            model.K[model.drawnrateindex]=model.mean
            plt.title('IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font)
        else:
            plt.title('IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font)
        if model.birthdist==None:
            plt.xlabel(r'Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font)
        else:
            plt.xlabel(r'Dis=%s, Obs=%s, Var=%s, ObsVar=%s, Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.birthdist,model.observations,model.var,model.obsvar,model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font)
        plt.ylabel('Slope of COV Spec=%s'%(model.ICspec[i]),**axis_font)
        fig = plt.gcf()
        plt.show()
        fig.savefig(results_dir + model.reac + 'DerivativeAutoCorrelationVsLag' + model.ICspec[i] + '.png', format='png')
    
def AverageOfAutoCorrelation(model):
    for i in range(len(model.indexarray)):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        line1,=ax1.plot(model.shift,model.aveautocorrelation[i],color='blue',label='<COV>')
        ax1.legend(handles=[line1],loc='best')
        if model.static==True and i==model.drawnrateindex:
            model.K[model.drawnrateindex]=model.mean
            plt.title('IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font)
        else:
            plt.title('IC = %s=%s, %s=%s'%(model.ICspec,model.IC,model.Kspec,model.K),**title_font)
        if model.birthdist==None:
            plt.xlabel(r'Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font)
        else:
            plt.xlabel(r'Dis=%s, Obs=%s, Var=%s, ObsVar=%s, Gil=%s, Step=%s, $\epsilon$=%s, $\omega$=%s, FixT=%s, $\tau$'%(model.birthdist,model.observations,model.var,model.obsvar,model.trials,model.density,model.trigcoef,model.angfreq,model.fixt),**axis_font)
        plt.ylabel('<COV> Spec=%s'%(model.ICspec[i]),**axis_font)
        fig = plt.gcf()
        plt.show()
        fig.savefig(results_dir + model.reac + 'DerivativeAutoCorrelationVsLag' + model.ICspec[i] + '.png', format='png')
    
        

    