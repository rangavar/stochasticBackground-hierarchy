# Stochastic GW Background: Direct Model Comparison

This script employs Hierarchical Bayesian statistics to discriminate between different semi-analytic galaxy/supermassive black hole coevolution models, using the stochastic gravitational background emitted by the superposition of the mergers of supermassive black hole binaries (SMBHB). Moreover, it uses the "detections" of individually resolved mergers to improve the constraints on each semi-analytic models (in combination with the constraints from the stochastic GW background). It is optimized for the semi-analytic SMBHB population realizations in Klein et al. 2016, which are based on galaxy/supermassive blackhole coevolution models in Barausse et al. 2012, which are available at: 

https://people.sissa.it/~barausse/catalogs.html

This script follows five steps:

1) Simulating a noisy realization of the stochastic GW background based on a given SMBHB population synthesis model. 
2) Discriminating between different semi-analytic models based on this simulated noisy-background. 
3) Simulating a population of detected individually resolved binaries from one of the SMBHB population synthesis models. 
4) Discrimination between models based on the simulated detection of individually resolved binaries. 
5) Combining the constraints from stages 2 and 4 to get the "joint" constraints. 

The level of noise can be changed at the beginning of the script.

This script returns the constraints from the stochastic background as "LStrain.eps" files, constraints from the individually resolved binaries as "LDistribution.eps" file, and the joint constraints as "LFinal.eps", all in the ./plots directory. 

More details on the methods and the Hierarchical approach is available in the thesis at: 

https://github.com/rangavar/stochasticBackground-hierarchy/raw/master/Stochastic_Gravitational_Wave_Background_from_the_Mergers_of_Supermassive_Black_Hole_Binaries.pdf
