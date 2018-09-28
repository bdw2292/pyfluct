import matplotlib.pyplot as plt
import NoiseAlgorithm as N
import Models
import Plots as Plots

def SingleColorStaticOscillating_Q(M): 
    N.GillespieSim(M)
    N.SingleColorLoop(M)
    N.AutoCorrelation(M)
    Plots.TotalNoiseComparison(M)
    Plots.SingleColorNoise(M)
    Plots.AutoCorrelationVsTau(M)
    Plots.PredictedVsSimulatedIntNoise(M)
    Plots.PredictedVsSimulatedTotNoise(M)
    
def TwoColor_Q(M):
    N.TwoColorLoop(M)
    Plots.TotalNoiseComparison(M)
    Plots.TwoColorNoise(M)
    
def SingleAndTwoColorComparison_Q(M):   
    M.static=True
    M.benchmark=True
    M.SingleColorStaticOscillating_Q()
    N.GillespieSim(M)
    benchmarkintnoise,benchmarkextnoise=N.SingleColorLoop(M)
    M.benchmark=False
    M.TwoColorStatic_Q()
    M.benchmarkintrinsic=benchmarkintnoise
    M.benchmarkextrinsic=benchmarkextnoise
    N.TwoColorLoop(M)
    N.CheckIfSameNoise(M)
    M.color=True
    Plots.ExtNoiseMethodComparison(M)
    Plots.IntNoiseMethodComparison(M)
    
def SingleColorStaticNoiseRatioTest_QR(M):
    N.GillespieSim(M)
    N.SingleColorLoop(M)
    Plots.NoiseRatioTest(M)
    
def SingleColorStaticOscillating_QR(M):
    N.GillespieSim(M)
    N.SingleColorLoop(M)
    N.AutoCorrelation(M)
    N.AverageOfAutoCorrelation(M)
    Plots.AutoCorrelationVsTau(M)
    Plots.PredictedVsSimulatedIntNoise(M)
    Plots.PredictedVsSimulatedTotNoise(M)

def SingleColorStatic_SQR(M):
    N.GillespieSim(M)
    N.SingleColorLoop(M)
    Plots.TotalNoiseComparison(M)
    Plots.SingleColorNoise(M)
    
def Query():
    plt.close("all")
    question=input('Draw rates from a static distribution?(y/n) ')
    if question=='y':
        static=True
    else:
        static=False
    question=input('Sinusoidal birthrate?(y/n) ')
    if question=='y':
        osc=True
    else:
        osc=False
    return static,osc
    
def main(M): #Options and parameters such as birth-death constants
    print('Enter 1 for Single Color Static Oscillating Simulation Q ')
    print('Enter 2 for Two Color Static Simulation Q ')
    print('Enter 3 for Comparison of Two and Single Color Simulation Q')
    print('Enter 4 for Noise Ratio Test in QR Reaction via Single Color')
    print('Enter 5 for Static Oscillating Noise Propogation in QR Reaction via Single Color')
    print('Enter 0 to exit ')
    sim=int(input('Enter Choice '))
    if sim==1:
        M.SingleColorStaticOscillating_Q()
        SingleColorStaticOscillating_Q(M)
    elif sim==2:
        M.TwoColorStatic_Q()
        TwoColor_Q(M)
    elif sim==3:
        SingleAndTwoColorComparison_Q(M)
    elif sim==4:
        M.SingleColorStaticNoiseRatioTest_QR()
        SingleColorStaticNoiseRatioTest_QR(M)
    elif sim==5:
        M.SingleColorStaticOscillating_QR()
        SingleColorStaticOscillating_QR(M)
    elif sim==0:
        exit 
    
static,osc=Query()
M=Models.Model(static,osc)
main(M)