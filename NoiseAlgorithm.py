import numpy as np
from tqdm import tqdm

def AverageOfAutoCorrelation(model):
    if model.rates==None:
        model.rates=[model.K[model.drawnrateindex]]
    halfSS=int(len(model.t)*.5)
    Epbirth=np.empty((len(model.indexarray),len(model.rates),halfSS))
    E1birth=np.empty((len(model.indexarray),len(model.rates),halfSS))
    E2birth=np.empty((len(model.indexarray),len(model.rates),halfSS))
    if model.singlegill==True:
        cells=model.trials*len(model.rates)
    else:
        cells=len(model.rates)
    for s in tqdm(range(len(model.indexarray))):
        pos=model.indexarray[s]
        for k in range(cells):
            if model.singlegill==True:
                currentbirthrate=np.random.poisson(model.mean,1)[0] #birth constant
                model.K[model.drawnrateindex]=currentbirthrate #check if this updates the rate properly
                model.Simulate(1,model.tmax) #simulate based on rate constants
                model.state=model.data_points
                gridensemble=model.EnsembleGrid(model.state) #extract data at regular intervals of time (rows are trials,columns are time)
                model.state=np.array(list(map(list, zip(*gridensemble))))
            else:
                grid=model.gilarray[k]
            for tau in range(0,halfSS):
                flat1=grid[halfSS,:,pos]
                flat2=grid[halfSS+tau,:,pos]
                Epbirth[s,k,tau]=np.mean(flat1*flat2)
                E1birth[s,k,tau]=np.mean(flat1)
                E2birth[s,k,tau]=np.mean(flat2)
    pave=[np.mean(Epbirth[l],axis=0) for l in range(len(model.indexarray))]
    p1=[np.mean(E1birth[l],axis=0) for l in range(len(model.indexarray))]
    p2=[np.mean(E2birth[l],axis=0) for l in range(len(model.indexarray))]
    model.aveautocorrelation=[(pave[l]-p1[l]*p2[l]) for l in range(len(model.indexarray))]
    
def AutoCorrelation(model):
    if model.rates==None:
        model.rates=[model.K[model.drawnrateindex]]
    halfSS=int(len(model.t)*.5)
    signal1taubirth=np.empty((len(model.indexarray),len(model.rates),halfSS,model.trials))
    signal2taubirth=np.empty((len(model.indexarray),len(model.rates),halfSS,model.trials))
    for s in tqdm(range(len(model.indexarray))):
        pos=model.indexarray[s]
        for k in range(len(model.rates)): 
            grid=model.gilarray[k]
            for tau in range(0,halfSS):
                signal1taubirth[s,k,tau,:]=grid[halfSS,:,pos]
                signal2taubirth[s,k,tau,:]=grid[halfSS+tau,:,pos]
    model.autocorrelation=np.empty((len(model.indexarray),halfSS))
    for s in range(len(model.indexarray)):
        for tau in range(0,halfSS):
            sig1birth=np.matrix(signal1taubirth[s,:,tau,:])
            sig2birth=np.matrix(signal2taubirth[s,:,tau,:])
            flat1=np.array(sig1birth.flatten())
            flat2=np.array(sig2birth.flatten())
            a = np.mean(flat1[0]*flat2[0])-np.mean(flat1[0])*np.mean(flat2[0])
            model.autocorrelation[s,tau]=a
               
def CheckIfSameNoise(model): #The residual better be <=tolerance at every point
    for i in range(len(model.extnoise)):
        extdiff=np.abs(np.subtract(model.extnoise[i],model.benchmarkextrinsic[i]))
        for element in extdiff:
            if not (element <= model.tolerance):
                if model.interactive==None or model.interactive==True:
                    print('Extrinsic Noise is not the same! ')
                return False,extdiff
        if model.interactive==None or model.interactive==True:
            print('Extrinsic Noise is the same!')
            print('Tolerance =',model.tolerance)
    return (True,extdiff)
    
def GillespieSim(model):
    if model.rates==None:
        model.rates=[model.K[model.drawnrateindex]]
    model.gilarray=[]
    print('Generating Gillepie')
    for k in tqdm(range(len(model.rates))):
        currentbirthrate=model.rates[k] #birth constant
        model.K[model.drawnrateindex]=currentbirthrate #check if this updates the rate properly
        model.Simulate(model.trials,model.tmax) #simulate based on rate constants
        model.state=model.data_points
        gridensemble=model.EnsembleGrid(model.state) #extract data at regular intervals of time (rows are trials,columns are time)
        model.state=np.array(list(map(list, zip(*gridensemble)))) #transpose matrix, for manipulation convience 
        model.gilarray.append(model.state)
    model.gilarray=np.array(model.gilarray)
        
def SingleColorLoop(model): #rates is a list of rate constants, index is the location of species you wish to compute noise on 
    if model.rates==None:
        model.rates=[model.K[model.drawnrateindex]]
    extofintaveofp=np.empty((len(model.indexarray),len(model.rates),len(model.t)))
    extofsqrofintaveofp=np.empty((len(model.indexarray),len(model.rates),len(model.t)))
    extofintnum=np.empty((len(model.indexarray),len(model.rates),len(model.t)))
    extofintaveofsqrofp=np.empty((len(model.indexarray),len(model.rates),len(model.t)))
    model.interactive=True
    for s in tqdm(range(len(model.indexarray))):
        pos=model.indexarray[s]
        for k in range(len(model.rates)): 
            grid=model.gilarray[k]
            for t in range(len(grid)): #For every point in time
                p1=np.mean(grid[t,:,pos])
                p2=np.mean(np.square(grid[t,:,pos]))#mean of all p**2 at time t across every trial
                extofintaveofp[s,k,t]=p1
                extofsqrofintaveofp[s,k,t]=p1**2
                extofintnum[s,k,t]=p2-p1**2
                extofintaveofsqrofp[s,k,t]=p2
        #average over all possible rates drawn from a poisson, aka the extrinsic average
    intrinsicnum=[np.mean(extofintnum[l],axis=0) for l in range(len(model.indexarray))]
    model.extaveofintaveofp=[np.mean(extofintaveofp[l],axis=0) for l in range(len(model.indexarray))]
    sqrofextaveofintaveofp=[np.square(model.extaveofintaveofp[l]) for l in range(len(model.indexarray))]
    extaveofsqrofintaveofp=[np.mean(extofsqrofintaveofp[l],axis=0) for l in range(len(model.indexarray))]
    extaveofintaveofsqrofp=[np.mean(extofintaveofsqrofp[l],axis=0) for l in range(len(model.indexarray))]
    if model.noisescale==True or None:
        model.intnoise=[intrinsicnum[l]/sqrofextaveofintaveofp[l] for l in range(len(model.indexarray))]
        model.extnoise=[(extaveofsqrofintaveofp[l]-sqrofextaveofintaveofp[l])/sqrofextaveofintaveofp[l] for l in range(len(model.indexarray))]
        model.testTotnoise=[(extaveofintaveofsqrofp[l]-sqrofextaveofintaveofp[l])/sqrofextaveofintaveofp[l] for l in range(len(model.indexarray))]
        model.totnoise=[model.intnoise[l]+model.extnoise[l] for l in range(len(model.indexarray))]
    else:
        model.intnoise=[intrinsicnum[l] for l in range(len(model.indexarray))]
        model.extnoise=[(extaveofsqrofintaveofp[l]-sqrofextaveofintaveofp[l]) for l in range(len(model.indexarray))]
        model.testTotnoise=[(extaveofintaveofsqrofp[l]-sqrofextaveofintaveofp[l]) for l in range(len(model.indexarray))]
        model.totnoise=[model.intnoise[l]+model.extnoise[l] for l in range(len(model.indexarray))]
    if model.benchmark==True: 
        return(model.intnoise,model.extnoise)
            
def DifferenceofCopy1Copy2(model,state,index1,index2):
    ensemble=[]
    for i in range(len(state)):
        ensemble.append([])
        trajectory=state[i]
        for j in range(len(trajectory)):
            array=trajectory[j]
            copy1=array[index1]
            copy2=array[index2]
            diff=[copy1-copy2]
            ensemble[i].append(diff)
    return ensemble
    
def TwoColorMean(model,ensemble):
    mvals = np.empty((len(model.t),len(model.IC)))
    for i in range(len(model.t)):
        d = ensemble[i]
        mvals[i,:] = np.mean(d, axis = 0)
    return mvals
        
def TwoColorVariance(model,ensemble):
    mvals = np.empty((len(model.t),len(model.IC)))
    for i in range(len(model.t)):
        d = ensemble[i]
        mvals[i,:] = np.var(d, axis = 0)
    return mvals

def TwoColorLoop(model):
    model.extnoise=np.empty((len(model.indexarray),len(model.t)))
    model.intnoise=np.empty((len(model.indexarray),len(model.t)))
    model.totnoise=np.empty((len(model.indexarray),len(model.t)))
    model.testTotnoise=np.empty((len(model.indexarray),len(model.t)))
    ensemble=[]
    for k in tqdm(range(len(model.rates))): #For a given birth constant
        currentbirthrate=model.rates[k] #birth constant
        model.K[model.drawnrateindex]=currentbirthrate 
        model.Simulate(1,model.tmax) #simulate both copies
        model.state=model.data_points
        gridensemble=model.EnsembleGrid(model.state)
        ensemble.append(gridensemble[0])
    newensemble=np.array(list(map(list, zip(*ensemble))))
    for s in range(len(model.indexarray)):
        pair=model.indexarray[s]
        indexcopy1=pair[0]
        indexcopy2=pair[1]
        copy12covariance=model.Covariance(newensemble)[indexcopy1][indexcopy2]
        meandata=TwoColorMean(model,newensemble) #means of data
        copy1mean=meandata[:,indexcopy1] #copy1 mean
        copy2mean=meandata[:,indexcopy2] #copy2 mean
        diffensemble=DifferenceofCopy1Copy2(model,newensemble,indexcopy1,indexcopy2)
        vardiffdata=TwoColorVariance(model,diffensemble)
        VarianceOfDiffOfCopy1Copy2=vardiffdata[:,0]
        vardata=TwoColorVariance(model,newensemble)
        var1=vardata[:,indexcopy1]
        var2=vardata[:,indexcopy2]
        totalvar=np.add(var1,var2)
        if model.noisescale==True or None:
            extn=copy12covariance/(copy1mean*copy2mean)
            intn=VarianceOfDiffOfCopy1Copy2/(2*copy1mean*copy2mean)
            model.extnoise[s,:]=extn #use the formulas
            model.intnoise[s,:]=intn
            model.totnoise[s,:]=totalvar/(2*copy1mean*copy2mean)
            model.testTotnoise[s,:]=intn+extn
        else:
            extn=copy12covariance
            intn=VarianceOfDiffOfCopy1Copy2/2
            model.extnoise[s,:]=extn #use the formulas
            model.intnoise[s,:]=intn
            model.totnoise[s,:]=totalvar/2
            model.testTotnoise[s,:]=intn+extn
    if model.benchmark==True: 
        return(model.intnoise,model.extnoise)
        
