# SchiSTOP model implementation

#### SchiSTOP (Beyond mass drug administration: understanding Schistosomiasis dynamics to STOP transmission), 

a project funded by the European Marie Skłodowska-Curie fellowship.

#

Repository of code for SchiSTOP modelling framework: a stochastic agent-based (ABM) transmission model for the transmission dynamics of *Schistosoma mansoni* in human population. An ODE-based module for the infection dynamics within the snail intermidiate host is integrated into the ABM.

Detailed description of SchiSTOP has been published within the manuscript: *"Revisiting the impact of Schistosoma mansoni regulating mechanisms on transmission dynamics using SchiSTOP , a novel modelling framework"*.

### **Mode of use**

1.  File `02.2_Setting_simulation_scenario.R`: customize scenario-specific parameters. This will prepare the setting to run simulations.

2.  File `05_Run_model.R` : run the model and save results. The script allows to load input data, call previous scripts, customize the parameters for simulations, launch simulations with SchiSTOP from R source (`04_Model_specification.R`). Results are stored in two outputs:

    -   population-level output, i.e. results aggregated at population level (e.g. prevalence)

    -   individual-level output (optional), i.e. results tracked over time for each individuals

### **Author**

Veronica Malizia<sup>1</sup>

### **Supervision**

Sake J. de Vlas<sup>2</sup>, Federica Giardina<sup>1*</sup>

<sup>1</sup> Radboud University Medical Center, Department IQ Health, Biostatistics Research Group, Nijmegen, The Netherlands

<sup>2</sup> Department of Public Health, Erasmus MC, University Medical Center Rotterdam, Rotterdam, The Netherlands

<sup>*</sup> Project leader
