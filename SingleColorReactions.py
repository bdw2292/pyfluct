def Birth_R(state):
    x = state[:] 
    x[1]+=1
    return x
        
def BirthRate_R(state,KB):
    return KB[2]*state[0]
        
def Death_R(state):
    x = state[:]
    x[1] -= 1
    return x
        
def DeathRate_R(state,KB):
    return KB[3]*state[1]
        
def Birth_S(state):
    x = state[:] 
    x[2]+=1
    return x
        
def BirthRate_S(state,KB):
    return KB[4]
        
def Death_S(state):
    x = state[:]
    x[2] -= 1
    return x
        
def DeathRate_S(state,KB):
    return KB[5]*state[2]
        
def Birth_Q(state):
    x = state[:] 
    x[0]+=1
    return x
    
def BirthRate_QFromS(state,KB):
    return KB[0]*state[2]
        
def BirthRate_Q(state,KB):    
    return KB[0]
        
def Death_Q(state):
    x = state[:]
    x[0] -= 1
    return x
        
def DeathRate_Q(state,KB):
    return KB[1]*state[0]
        

        
