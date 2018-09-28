def Birth_Q_1(state):
    x = state[:] 
    x[0]+=1
    return x
        
def BirthRate_Q_1(state,KB):
    return KB[0]
        
def Death_Q_1(state):
    x = state[:]
    x[0] -= 1
    return x
        
def DeathRate_Q_1(state,KB):
    return KB[1]*state[0]
        
def Birth_Q_2(state):
    x = state[:] 
    x[1]+=1
    return x
        
def BirthRate_Q_2(state,KB):
    return KB[0]
        
def Death_Q_2(state):
    x = state[:]
    x[1] -= 1
    return x
        
def DeathRate_Q_2(state,KB):
    return KB[1]*state[1]
        
def Birth_R_1(state):
    x = state[:] 
    x[2]+=1
    return x
        
def BirthRate_R_1(state,KB):
    return KB[2]*state[0]
        
def Death_R_1(state):
    x = state[:]
    x[2] -= 1
    return x
        
def DeathRate_R_1(state,KB):
    return KB[3]*state[2]
        
def Birth_R_2(state):
    x = state[:] 
    x[3]+=1
    return x
        
def BirthRate_R_2(state,KB):
    return KB[2]*state[1]
        
def Death_R_2(state):
    x = state[:]
    x[3] -= 1
    return x
        
def DeathRate_R_2(state,KB):
    return KB[3]*state[3]
        
  