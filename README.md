Hodgkin-Huxley Hippocampal Model
================================

Impliments simulated network of sparsely connected Hodgkin-Huxley type neurons as coupled differential equations reproducing the results in:

### *Wang, X. J., & Buzsáki, G. (1996). Gamma oscillation by synaptic inhibition in a hippocampal interneuronal network model. Journal of Neuroscience, 16(20), 6402–6413. https://doi.org/10.1523/jneurosci.16-20-06402.1996*

Simulations support the hypothesis that synchronous gamma oscillations (20-80 Hz) can emerge from a random network of interconnected GABAergic fast-spiking interneurons.  


**Simulation illustrates that:**

1. Network synchronization occurs only within a frequency band coinciding with the gamma (20–80 Hz) range. 

2. Large-scale network synchronization requires a critical (minimal) average number of synaptic contacts per cell, independent of the network size. 

3. Firing frequencies vary with maximal GABAa synaptic conductance, synaptic decay time constant, and mean external excitatory drive to the network.



Manual
=========

Requires Matlab. 

**Usage:**

Run HHsnetwork.m

To simulate network behavior as a function of number of synaptic inputs per neuron: 

Run MSyn_anal.m 

To simulate network behavior as a function of GABAa synaptic conductance:
 
Run ESyn_anal.m  

parameters:

* N - Number of cells in the network
* Msyn - fixed average number of synaptic inputs per neuron
* p = Msyn/N - probability of two neurons being connected
* Esyn - Synaptic potential (mv)
* alpha - activation rate constant (msec^-1)
* beta - inactivation rate constant (msec^-1)
* Iu - mean applied current (mA/cm^2)
* Isigma - standard deviation of applied current (mA/cm^2)
* T - simulation duration in (ms)


