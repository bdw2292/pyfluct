*************************************************************************************************
PyFluct - a simple python library for investigating stochastic fluctuations in chemical reactions
Brandon Walker
*************************************************************************************************

Models.py - Python class for a model of chemical reactions (coupled or a series of chemical reactions)
          - Input:
               a list of chemical reaction functions, and an array of rate constants, booleans to specify static rate constants or oscilatting rate constants
	  - Contains Gillespie algorithm
          - Contains example models 
            TwoColorStatic_Q, SingleColorStaticOscillating_Q, SingleColorStaticNoiseRatioTest_QR, SingleColorStaticOscillating_QR, SingleColorStatic_SQR

NoiseAlgorithm.py  - Contains functions for computing single color and dual color reporter noise


Parameters.py - input parameters, this  class is inherited by model class


Plots.py - Contains functions for plotting results of gillespie simulations

SingleColorReactions.py - Function definitions for single color chemical reactions

TwoColorReactions.py - Function definitions for dual color single reactions

SpecialCases.py - This python script is used as an example to test stochastic fluctuations in several models of chemical reactions
	Included Options:
    		Single Color Static Oscillating Simulation Q 
    		Two Color Static Simulation Q 
    		Comparison of Two and Single Color Simulation Q
    		Noise Ratio Test in QR Reaction via Single Color
    		Static Oscillating Noise Propogation in QR Reaction via Single Color


Example Plots are in Plots folder


************************************************************************************
Papers used for algorithms 

Intrinsic and extrinsic contributions to stochasticity in gene expression 
*************************************************************************************
