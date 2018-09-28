from Parameters import Parameters as P
import numpy as np
import functools as ft
import SingleColorReactions as SR
import TwoColorReactions as TR
import scipy as sp
from tqdm import tqdm
import multiprocessing as mp

class Model(P):
    
    def __init__(self,static,osc):
        self.reactions = [] 
        self.reacrates = []
        self.static=static
        self.osc=osc
        super().__init__()
    
    def Add_Reaction(self,f_reaction,f_rate):
        self.reactions.append(f_reaction)
        self.reacrates.append(f_rate)
        
    def Get_Rates(self,data):
        return np.array([f(data,self.K) for f in self.reacrates])
            
    def Total_Rate(self,rates):
        total=ft.reduce(lambda x,y: x+y,rates)
        return total
        
    def Normalize(self,rates,total_rate):
        norm_rates=[]
        for rate in rates:
            norm=rate/total_rate
            norm_rates.append(norm)
        return norm_rates

    def Fpt_Satisfied(self,fpt_thresh,data):
        if fpt_thresh is None:
            return True
        checks = [fpt_thresh[i] == data[i] for i in range(len(data)) if not fpt_thresh[i] is None]
        return ft.reduce(lambda x,y: x and y, checks)

    def Simulate(self,trials,tmax,fpt_thresh=None):
        self.data_points = [[self.IC] for i in range(trials)]
        self.time_points = [[0.0] for i in range(trials)]
        self.trials=trials
        self.N=len(self.IC)
        #jobs=[]
        if self.static==True:
            #for i in range(trials):
            #    p = mp.Process(target=self.InnerForLoop(i,tmax,fpt_thresh)) #Create a process or job defined by ParallelLoop on the ensemble
            #    jobs.append(p) #append job to a list to keep track of them
            #    p.start() #initialize job
            #for i in range(len(jobs)): #Want to first initialize all jobs (above), then let all finsish loop before moving on
            #    p=jobs[i] #current job
            #    p.join() #wait to finish the job
            for i in range(trials):
                self.InnerForLoop(i,tmax,fpt_thresh)
        else:
            #for i in tqdm(range(trials)):
            #    p = mp.Process(target=self.InnerForLoop(i,tmax,fpt_thresh)) #Create a process or job defined by ParallelLoop on the ensemble
            #    jobs.append(p) #append job to a list to keep track of them
            #    p.start() #initialize job
            #for i in range(len(jobs)): #Want to first initialize all jobs (above), then let all finsish loop before moving on
            #    p=jobs[i] #current job
            #    p.join() #wait to finish the job
            for i in tqdm(range(trials)):
                self.InnerForLoop(i,tmax,fpt_thresh)
        return(self.data_points,self.time_points)
        
    def InnerForLoop(self,i,tmax,fpt_thresh):
        time = 0.0
        data = self.IC
        fpt_sat = self.Fpt_Satisfied(fpt_thresh,data)
        while time < tmax or not fpt_sat:
            if self.osc==False:
                rates = self.Get_Rates(data)
                total_rate = self.Total_Rate(rates)
                time_step=np.random.exponential(1.0/total_rate)
                time += time_step
                norm_rates = self.Normalize(rates,total_rate)
            else:
                rates = self.Get_Rates(data)
                rates = np.delete(rates,self.drawnrateindex)
                propen = self.Total_Rate(rates)
                r=np.random.uniform(0,1)
                #def Integrand(z):
                #    return self.OscPropensity(time)+propen
                #def NumericalFunc(x):
                #    integral,err=sp.integrate.quad(Integrand, time, time+x)
                #    return 1-r-np.exp(-integral)
                def AnalyticalFunc(x):
                    return 1-r-np.exp(-((self.Q_b+propen)*x+((self.Q_b*self.trigcoef)/self.angfreq)*(np.cos(self.angfreq*time)-np.cos(self.angfreq*(time+x)))))
                vfunc=np.vectorize(AnalyticalFunc)
                initialguess=-np.log(1-r)/(self.Q_b*(1+self.trigcoef*np.sin(self.angfreq*time))+propen)
                #roots=sp.optimize.root(vfunc,initialguess,tol=10^-8)
                roots=sp.optimize.fsolve(vfunc, initialguess, xtol=1.49012e-08)
                #time_step=roots.x[0]
                time_step=roots[0]
                time += time_step
                drawnrate=self.K[self.drawnrateindex]*(1+self.trigcoef*np.sin(self.angfreq*(time+time_step)))
                rates = np.insert(rates,self.drawnrateindex,drawnrate)
                total_rate = self.Total_Rate(rates)
                norm_rates = self.Normalize(rates,total_rate)
            choice = np.nonzero(np.random.multinomial(1,norm_rates))[0][0]
            data = self.reactions[choice](data)
            self.data_points[i].append(data)
            self.time_points[i].append(time) 
            if not fpt_sat:
                fpt_sat = self.Fpt_Satisfied(fpt_thresh,time)
    
    def Get_Data(self,t,ensemble):
        data = np.empty((self.trials,len(self.IC)))
        for i in range(self.trials):
            di = 0
            while self.time_points[i][di] <= t:
                di += 1
            data[i,:] = ensemble[i][di-1]
        return data
        
    def Get_Moment(self,times,m,ensemble):
        mvals = np.empty((len(times),len(self.IC)))
        for i in range(len(times)):
            
            d = self.Get_Data(times[i],ensemble)
            mvals[i,:] = np.mean(d ** m, axis = 0)
        return mvals

    def Get_Fpt(self,fpt_thresh):
        fptimes = np.empty(self.trials)
        for i in range(self.trials):
            j = 0
            while not self.Fpt_Satisfied(fpt_thresh,self.data_points[i][j]):
                j += 1
            fptimes[i] = self.time_points[i][j]
        return fptimes
        
    def CopyNumberProbabilityMass(self,val,times):
        histogram = np.empty((self.trials,len(times)))
        for i in range(len(times)):
            results = self.Get_Data(times[i],self.data_points)
            histogram[:,i] = results[:,val]
        return histogram
    
    def MeanScaledProbabilityMass(self,val,times):
        histogram = np.empty((self.trials,len(times)))
        means = self.Get_Moment(times,1,self.data_points)
        for i in range(len(times)):
            results = self.Get_Data(times[i],self.data_points)
            histogram[:,i] = results[:,val]/means[i,val]
        return histogram
    
    def MeanRescaledEnsemble(self,state):
        data = []
        means = self.Get_Moment(self.t,1,state)
        for i in range(len(self.t)):
            results = self.Get_Data(self.t[i],state)
            data.append([])
            for k in range(self.N):
                results[:,k]=results[:,k]/means[i,k]
            for j in range(self.trials):
                obj=results[j]
                data[i].append(obj)
        newdata=[list(i) for i in zip(*data)]
        return newdata
        
    def EnsembleGrid(self,state):
        data = []
        for i in range(len(self.t)):
            results = self.Get_Data(self.t[i],state)
            data.append([])
            for j in range(self.trials):
                obj=results[j]
                data[i].append(obj)
        newdata=[list(i) for i in zip(*data)]
        return newdata
    
    def Correlation(self):
        cor_data = np.zeros((self.N,self.N,len(self.t)))
        for i in range(len(self.t)):
            d = self.Get_Data(self.t[i],self.data_points)
            cor_data[:,:,i] = np.corrcoef(d,rowvar=0)
        return cor_data
    
    def Covariance(self,state):
        cov_data = np.zeros((self.N,self.N,len(self.t)))
        for i in range(len(self.t)):
            d=state[i]
            cov_data[:,:,i] = np.cov(d,rowvar=0)
        return cov_data
        
    def Mean(self):
        mean=self.Get_Moment(self.t,1,self.data_points)
        return mean
            
    def Variance(self):
        variance=self.Get_Moment(self.t,2,self.data_points)-(self.Get_Moment(self.t,1,self.data_points))**2
        return variance
    
    def CoefficientOfVariation(self):
        means = self.Get_Moment(self.t,1,self.data_points)
        var=self.Variance()
        for i in range(len(self.t)):
           for val in range(self.N):
               var[i,val]=np.sqrt(var[i,val])/means[i,val]
        return var
        
    def DrawRateFromPoisson(self,lams,observations): 
        return np.random.poisson(lams, observations)
    
    def TwoColorStatic_Q(self):
        if self.osc==True:
            self.reac='StaticOscillatingQ'
            self.Kspec=['<Q_b(t)>','Q_d']
        else:
            self.reac='StaticQ'
            self.Kspec=['<Q_b>','Q_d']
            while self.static==True:
                self.rates=self.DrawRateFromPoisson(self.mean,self.observations)
                self.obsvar=np.var(self.rates)
                self.obsmean=np.mean(self.rates)
                if np.abs(self.obsvar-self.var)<=self.tolerance and np.abs(self.obsmean-self.mean)<=self.tolerance:
                    break
        self.drawnrateindex=0
        self.drawnstateindex=0
        self.drawnrate=self.Q_b
        self.trials=1
        self.observations=self.trials*self.observations
        self.IC=[self.IC_Q,self.IC_Q]
        self.ICspec=['(Q,Q)']
        self.reac='Q'
        self.K=np.array([self.Q_b,self.Q_d])
        self.Add_Reaction(TR.Birth_Q_1,TR.BirthRate_Q_1)
        self.Add_Reaction(TR.Death_Q_1,TR.DeathRate_Q_1)
        self.Add_Reaction(TR.Birth_Q_2,TR.BirthRate_Q_1)
        self.Add_Reaction(TR.Death_Q_2,TR.DeathRate_Q_2)
        self.indexarray=[[0,1]]
        self.predintnoise=[self.mean/self.Q_d]
        self.predextnoise=[self.var/self.Q_d**2]
        self.predtotnoise=[self.mean/self.Q_d+self.var/self.Q_d**2]
        halfSS=int(len(self.t)*.5)
        self.shift=[tau*self.density for tau in range(0,halfSS)]
        self.predautocor=[self.mean*np.exp(-self.Q_d*self.shift)/self.Q_d+self.var/self.Q_d**2]
   
    def SingleColorStaticOscillating_Q(self):
        if self.static==True and self.osc==False:
            while self.static==True:
                self.rates=self.DrawRateFromPoisson(self.mean,self.observations)
                self.obsvar=np.var(self.rates)
                self.obsmean=np.mean(self.rates)
                if np.abs(self.obsvar-self.var)<=self.tolerance and np.abs(self.obsmean-self.mean)<=self.tolerance:
                    break
            self.Kspec=['<Q_b>','Q_d']
            self.angfreq=0
            self.reac='StaticQ'
        elif self.static==False and self.osc==False:
            self.rates=None
            self.trials=self.trials*self.observations
            self.birthdist=None
            self.var=0
            self.mean=self.Q_b
            self.Kspec=['Q_b','Q_d']
            self.angfreq=0
            self.reac='Q'
        elif self.static==False and self.osc==True:
            self.rates=None
            self.trials=self.trials*self.observations
            self.birthdist=None
            self.var=0
            self.mean=self.Q_b
            self.Kspec=['Q_b(t)','Q_d']
            self.reac='OscillatingQ'
        elif self.static==True and self.osc==True:
            while self.static==True:
                self.rates=self.DrawRateFromPoisson(self.mean,self.observations)
                self.obsvar=np.var(self.rates)
                self.obsmean=np.mean(self.rates)
                if np.abs(self.obsvar-self.var)<=self.tolerance and np.abs(self.obsmean-self.mean)<=self.tolerance:
                    break
            self.Kspec=['<Q_b(t)>','Q_d']
            self.reac='StaticOscillatingQ'
        #self.OscPropensity=self.OscPropensityQ
        self.drawnrateindex=0
        self.drawnstateindex=0
        self.drawnrate=self.Q_b
        self.IC=[self.IC_Q]
        self.ICspec=['Q']
        self.reac='Q'
        self.K=np.array([self.Q_b,self.Q_d])
        self.Add_Reaction(SR.Birth_Q,SR.BirthRate_Q)
        self.Add_Reaction(SR.Death_Q,SR.DeathRate_Q)
        self.indexarray=[0]
        self.predintnoise=[(self.mean/self.Q_d)*(1+self.epsphi*np.sin(self.angfreq*self.t-self.phi))]
        self.predextnoise=[(self.var/self.Q_d**2)*(1+self.epsphi*np.sin(self.angfreq*self.t-self.phi)**2)]
        self.predtotnoise=[self.predintnoise[0]+self.predextnoise[0]]
        self.predautocor=[((self.mean*np.exp(-self.Q_d*self.diff)*(1+self.epsphi*np.sin(self.angfreq*self.fixt-self.phi)))/self.Q_d)+(self.var*(1+self.epsphi*np.sin(self.angfreq*self.fixt-self.phi))*(1+self.epsphi*np.sin(self.angfreq*self.tprime-self.phi)))/self.Q_d**2]
    
    def SingleColorStaticNoiseRatioTest_QR(self):    
        self.R_b=10
        self.R_d=.01
        self.Q_b=15
        self.Q_d=10
        self.mean=6
        self.rates=self.DrawRateFromPoisson(self.mean,self.observations)
        self.obsvar=np.var(self.rates)
        self.obsmean=np.mean(self.rates)
        self.drawnrateindex=0
        self.drawnstateindex=0
        self.drawnrate=self.Q_b
        self.ICspec=['Q','R']
        self.reac='StaticQR'
        self.IC=[self.IC_Q,self.IC_R]
        self.Kspec=['<Q_b>','Q_d','R_b','R_d']
        self.K=np.array([self.Q_b,self.Q_d,self.R_b,self.R_d])
        self.Add_Reaction(SR.Birth_Q,SR.BirthRate_Q)
        self.Add_Reaction(SR.Death_Q,SR.DeathRate_Q)
        self.Add_Reaction(SR.Birth_R,SR.BirthRate_R)
        self.Add_Reaction(SR.Death_R,SR.DeathRate_R)
        self.indexarray=[0,1]
        
    def SingleColorStaticOscillating_QR(self):
        if self.static==True and self.osc==False:
            while self.static==True:
                self.rates=self.DrawRateFromPoisson(self.mean,self.observations)
                self.obsvar=np.var(self.rates)
                self.obsmean=np.mean(self.rates)
                if np.abs(self.obsvar-self.var)<=self.tolerance and np.abs(self.obsmean-self.mean)<=self.tolerance:
                    break
            self.Kspec=['<Q_b>','Q_d','R_b','R_d']
            self.angfreq=0
            self.reac='StaticQR'
        elif self.static==False and self.osc==False:
            self.rates=None
            self.birthdist=None
            self.trials=self.trials*self.observations
            self.var=0
            self.mean=self.Q_b
            self.Kspec=['Q_b','Q_d','R_b','R_d']
            self.angfreq=0
            self.reac='QR'
        elif self.static==False and self.osc==True:
            self.rates=None
            self.trials=self.trials*self.observations
            self.birthdist=None
            self.var=0
            self.mean=self.Q_b
            self.Kspec=['Q_b(t)','Q_d','R_b','R_d']
            self.reac='OscillatingQR'
        elif self.static==True and self.osc==True:
            while self.static==True:
                self.rates=self.DrawRateFromPoisson(self.mean,self.observations)
                self.obsvar=np.var(self.rates)
                self.obsmean=np.mean(self.rates)
                if np.abs(self.obsvar-self.var)<=self.tolerance and np.abs(self.obsmean-self.mean)<=self.tolerance:
                    break
            self.Kspec=['<Q_b(t)>','Q_d','R_b','R_d']
            self.reac='StaticOscillatingQR'
        #self.OscPropensity=self.OscPropensityQ
        self.drawnrateindex=0
        self.drawnstateindex=0
        self.drawnrate=self.Q_b
        self.ICspec=['Q','R']
        self.IC=[self.IC_Q,self.IC_R]
        self.K=np.array([self.Q_b,self.Q_d,self.R_b,self.R_d])
        self.Add_Reaction(SR.Birth_Q,SR.BirthRate_Q)
        self.Add_Reaction(SR.Death_Q,SR.DeathRate_Q)
        self.Add_Reaction(SR.Birth_R,SR.BirthRate_R)
        self.Add_Reaction(SR.Death_R,SR.DeathRate_R)
        self.indexarray=[0,1]
        self.phi=np.arctan(self.angfreq/self.R_d)
        self.thephi=self.trigcoef*np.cos(self.theta)*np.cos(self.phi)
        self.intR=(self.R_b*self.mean*(1+self.thephi*np.sin(self.angfreq*self.t-self.theta-self.phi))/(self.R_d*self.Q_d))+(((self.R_b**2)*self.mean)/self.Q_d)*(((1+self.thepsidel*np.sin(self.angfreq*self.t-self.theta-self.psi-self.delta))/(2*self.R_d*(self.R_d+self.Q_d)))+((1+self.thepsi*np.sin(self.angfreq*self.t-self.theta-self.psi))/(self.R_d**2-self.Q_d**2))-((1+self.thedel*np.sin(self.angfreq*self.t-self.theta-self.delta))/(2*self.R_d*(self.R_d-self.Q_d))))
        self.extR=(((self.R_b**2)*self.var*(1+self.thephi*np.sin(self.angfreq*self.t-self.theta-self.phi))**2)/(self.R_d*self.Q_d)**2)
        self.totR=self.intR+self.extR
        self.predintnoise=[None,self.intR]
        self.predextnoise=[None,self.extR]
        self.predtotnoise=[None,self.totR]
        self.autocorR=(np.exp(-self.R_b*self.diff)*self.R_b*self.mean*(1+self.thephi*np.sin(self.angfreq*self.fixt-self.theta-self.phi))/(self.R_d*self.Q_d))+(((self.R_b**2)*self.mean)/self.Q_d)*((np.exp(-self.R_d*self.diff)*(1+self.thepsidel*np.sin(self.angfreq*self.fixt-self.theta-self.psi-self.delta))/(2*self.R_d*(self.R_d+self.Q_d)))+(np.exp(-self.Q_d*self.diff)*(1+self.thepsi*np.sin(self.angfreq*self.fixt-self.theta-self.psi))/(self.R_d**2-self.Q_d**2))-(np.exp(-self.R_d*self.diff)*(1+self.thedel*np.sin(self.angfreq*self.fixt-self.theta-self.delta))/(2*self.R_d*(self.R_d-self.Q_d))))+(((self.R_b**2)*self.var*(1+self.thephi*np.sin(self.angfreq*self.fixt-self.theta-self.phi))*(1+self.thephi*np.sin(self.angfreq*self.tprime-self.theta-self.phi)))/(self.R_d*self.Q_d)**2)
        self.predautocor=[None,self.autocorR]
  
    def SingleColorStatic_SQR(self):
        if self.static==True:
            while self.static==True:
                self.rates=self.DrawRateFromPoisson(self.mean,self.observations)
                self.obsvar=np.var(self.rates)
                self.obsmean=np.mean(self.rates)
                if np.abs(self.obsvar-self.var)<=self.tolerance and np.abs(self.obsmean-self.mean)<=self.tolerance:
                    break
            self.Kspec=['<Q_b>','Q_d','R_b','R_d','S_b','S_d']
            self.reac='StaticSQR'
        else:
            self.rates=None
            self.birthdist=None
            self.trials=self.trials*self.observations
            self.var=0
            self.mean=self.S_b
            self.Kspec=['Q_b','Q_d','R_b','R_d','S_b','S_d']
            self.reac='SQR'
        self.drawnrateindex=4
        self.drawnstateindex=2
        self.drawnrate=self.S_b
        self.trials=self.observations*self.trials
        self.ICspec=['Q','R','S']
        self.IC=[self.IC_Q,self.IC_R,self.IC_S]
        self.Kspec=['Q_b','Q_d','R_b','R_d','S_b','S_d']
        self.K=np.array([self.Q_b,self.Q_d,self.R_b,self.R_d,self.S_b,self.S_d])
        self.Add_Reaction(SR.Birth_Q,SR.BirthRate_QFromS)
        self.Add_Reaction(SR.Death_Q,SR.DeathRate_Q)
        self.Add_Reaction(SR.Birth_R,SR.BirthRate_R)
        self.Add_Reaction(SR.Death_R,SR.DeathRate_R)
        self.Add_Reaction(SR.Birth_S,SR.BirthRate_S)
        self.Add_Reaction(SR.Death_S,SR.DeathRate_S)
        self.indexarray=[2,0,1]
        self.intR=(self.R_b*self.Q_b*self.mean)/(self.R_d*self.Q_d*self.S_d)
        self.extR=(self.R_b*self.Q_b*self.mean)/(self.R_d*self.Q_d*self.S_d)*(self.R_b/(self.Q_d+self.R_d))+(self.R_b*self.Q_b*self.mean)/(self.R_d*self.Q_d*self.S_d)*(self.R_b*self.Q_b*(self.R_d+self.Q_d+self.S_d)/((self.Q_d+self.R_d)*(self.S_d+self.Q_d)*(self.S_d+self.R_d))) + (self.R_b**2)*(self.Q_b**2)*self.var/(self.R_d*self.Q_d*self.S_d)**2
        self.totR=self.intR+self.extR
        self.predintnoise=[None,None,self.intR]
        self.predextnoise=[None,None,self.extR]
        self.predtotnoise=[None,None,self.totR]
        self.autocorR=(self.R_b*self.Q_b*self.mean)/(self.R_d*self.Q_d*self.S_d)*np.exp(-self.R_d*self.shift)+(self.R_b**2*self.Q_b*self.mean/self.R_d*self.Q_d*self.S_d)*(self.Q_d/self.Q_d**2-self.R_d**2)*(self.Q_d*np.exp(-self.R_d*self.shift)-self.R_d*np.exp(-self.Q_d*self.shift))
        self.autocorR+=(self.R_b**2*self.Q_b*self.mean/self.R_d*self.Q_d*self.S_d)*((self.Q_b*self.S_d/self.S_d**2-self.Q_d**2)*(self.Q_d/(self.Q_d**2-self.R_d**2))*(self.Q_d*np.exp(-self.R_d*self.shift)-self.R_d*np.exp(-self.Q_d*self.shift))-(self.Q_b*self.Q_d/(self.S_d**2-self.Q_d**2))*(self.S_d/(self.S_d**2-self.R_d**2))*(self.S_d*np.exp(-self.R_d*self.shift)-self.R_d*np.exp(-self.S_d*self.shift)))+((self.R_b*2*self.Q_b**2*self.var)/(self.R_d*self.Q_d*self.S_d))**2
        self.predautocor=[None,None, np.full((len(self.shift),),self.autocorR)]
        