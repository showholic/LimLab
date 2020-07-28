import numpy as np
import math
import sys
import pdb
import subprocess 
##########################################################
# Parameter sets for the calcium and synaptic weight dynamics
class synapticChange2():
    ''' 
        class to calculate the change in synaptic strenght 
    '''
    ###############################################################################
    def __init__(self, tauCa,Cpre,Cpost,thetaD,thetaP,gammaD,gammaP,sigma,tau,rhoStar,D,beta,b,nonlinear=1.):
        
        # chose parameters from predefined set or from file
        self.tauCa = tauCa
        self.Cpre  = Cpre
        self.Cpost = Cpost
        self.thetaD = thetaD
        self.thetaP = thetaP
        self.gammaD = gammaD
        self.gammaP = gammaP 
        self.sigma  = sigma
        self.tau    = tau 
        self.rhoStar= rhoStar
        self.D      = D
        self.beta   = beta
        self.b      = b
        

        
        

    ##########################################################
    # calculate UP and DOWN transition probabilities
    def transitionProbability(self, T_total, rho0, rhoBar,sigmaRhoSquared,tauEff):
        
        # argument for the Error Function
        x1 =  -(self.rhoStar - rhoBar + (rhoBar - rho0)*np.exp(-T_total/tauEff))/(np.sqrt(sigmaRhoSquared*(1.-np.exp(-2.*T_total/tauEff))))
        # transition probability
        if rho0 == 0.:
            return (1. + math.erf(x1))/2.
        else:
            return (1. - math.erf(x1))/2.



    ##########################################################
    # calculate all values for the change in synaptic strength
    def changeInSynapticStrength2(self, T_total,rho0,deltaT,preRate,postRate,ppp,deltaCa):
        if self.Cpre>self.Cpost:
                arguments = str(deltaT) + ' ' + str(self.tauCa) + ' ' + str(self.Cpost) + ' ' + str(self.Cpre) + ' ' + str(self.thetaD) + ' ' + str(self.thetaP) + ' ' + str(preRate) + ' ' + str(postRate) + ' ' + str(ppp) + ' ' + str(deltaCa)
        else:
                arguments = str(deltaT) + ' ' + str(self.tauCa) + ' ' + str(self.Cpre) + ' ' + str(self.Cpost) + ' ' + str(self.thetaD) + ' ' + str(self.thetaP) + ' ' + str(preRate) + ' ' + str(postRate) + ' ' + str(ppp) + ' ' + str(deltaCa)

        #print arguments
        #(out,err) = commands.getstatusoutput('./timeAboveThreshold/poissonPairs_timeAboveThreshold ' + arguments) 
        (out,err) = subprocess.getstatusoutput('~/limlab/CalciumBasedPlasticityModel/timeAboveThreshold/poissonPairs_timeAboveThreshold ' + arguments) 
        alphaD = float(err.split()[0])
        alphaP = float(err.split()[1])
        # fraction of time spent above threshold
        (self.alphaD,self.alphaP) = (alphaD, alphaP)
        
        # average potentiation and depression rates
        self.GammaP = self.gammaP*self.alphaP
        self.GammaD = self.gammaD*self.alphaD
        
        # rhoBar: average value of rho in the limit of a very long protocol equivalent to the minimum of the quadratic potentia
        self.rhoBar = self.GammaP/(self.GammaP + self.GammaD)
        
        # sigmaRhoSquared L standard deviation of rho in the same limit
        self.sigmaRhoSquared = (self.alphaP + self.alphaD)*(self.sigma**2)/(self.GammaP + self.GammaD)
        # tauEff : characteristic time scale of the temporal evolution of the pdf of rho
        self.tauEff = self.tau/(self.GammaP + self.GammaD)
        #
        # UP and the DOWN transition probabilities
        self.UP   = self.transitionProbability(T_total,0.,self.rhoBar,self.sigmaRhoSquared,self.tauEff)
        self.DOWN = self.transitionProbability(T_total,1.,self.rhoBar,self.sigmaRhoSquared,self.tauEff)
        
        # mean value of the synaptic strength right at the end of the stimulation protocol
        self.meanUP   =  self.rhoBar - (self.rhoBar - 0.)*np.exp(-T_total/self.tauEff)
        self.meanDOWN =  self.rhoBar - (self.rhoBar - 1.)*np.exp(-T_total/self.tauEff)
        self.mean     =  self.rhoBar - (self.rhoBar - rho0)*np.exp(-T_total/self.tauEff)
        self.var      =  self.sigmaRhoSquared * (1-np.exp(-2*T_total/self.tauEff))/2
        # change in synaptic strength after/before
        self.synChange = ((self.beta*(1.-self.UP) + (1.-self.beta)*self.DOWN) + (self.beta*self.UP+ (1.-self.beta)*(1.-self.DOWN))*self.b)/(self.beta + (1.-self.beta)*self.b)
        
        #synChange = synapticChange.changeSynapticStrength(synapticChange.beta,UP,DOWN,synapticChange.b)
        
        
    

    

        
