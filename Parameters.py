import numpy as np

class Parameters():
    def __init__(self):
        self.singlegill=False
        self.noisescale=False
        self.birthdist='Pois'
        self.R_b=7
        self.R_d=6
        self.Q_b=6
        self.Q_d=5
        self.S_b=3
        self.S_d=1
        self.IC_R=5
        self.IC_Q=3
        self.IC_S=2
        self.mean=5
        self.var=self.mean
        self.trials = 100
        self.tmax = 5
        self.fixt=self.tmax/2
        self.observations= 100
        self.tolerance=.001
        self.density=.01
        self.steps=(self.tmax-0)/self.density
        self.t=np.linspace(0,self.tmax,self.steps)
        self.plot=True
        self.interactive=True
        self.benchmark=True
        self.color=False
        self.angfreq=3
        self.trigcoef=-.5
        self.halfSS=int(len(self.t)*.5)
        self.shift=[tau*self.density for tau in range(0,self.halfSS)]
        self.tprime=np.array(self.shift)+self.fixt
        self.diff=self.tprime-self.fixt
        self.phi=np.arctan(self.angfreq/self.Q_d)
        self.delta=np.arctan(self.angfreq/(2*self.R_d))
        self.psi=np.arctan(self.angfreq/(self.Q_d+self.R_d))
        self.theta=np.arctan(self.angfreq/self.Q_d)
        self.epsphi=self.trigcoef*np.cos(self.phi)
        self.thephi=self.trigcoef*np.cos(self.theta)*np.cos(self.phi)
        self.thepsidel=self.trigcoef*np.cos(self.theta)*np.cos(self.psi)*np.cos(self.delta)
        self.thepsi=self.trigcoef*np.cos(self.theta)*np.cos(self.psi)
        self.thedel=self.trigcoef*np.cos(self.theta)*np.cos(self.delta)
        #self.OscPropensityQ=lambda time:self.Q_b*(1+np.sin(self.angfreq*time))
