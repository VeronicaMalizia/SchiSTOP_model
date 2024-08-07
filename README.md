# SchiSTOP model implementation

Repository of code for SchiSTOP modelling framework: a stochastic agent-based (ABM) transmission model for the transmission dynamics of *Schistosoma mansoni* in human population. An ODE-based module for the infection dynamics within the snail intermidiate host is integrated into the ABM.

Detailed description of SchiSTOP has been published within the manuscript: *"Revisiting the impact of Schistosoma mansoni regulating mechanisms on transmission dynamics using SchiSTOP , a novel modelling framework"*.

### **Mode of use**

1.  File `02.2_Setting_simulation_scenario.R`: customize scenario-specific parameters. This will prepare the setting to run simulations.

2.  File `05_Run_model.R` : run the model and save results. The script allows to load input data, call previous scripts, customize the parameters for simulations, launched simulations with SchiSTOP from R source (`04_Model_specification.R`). Results are stored in two outputs:

    -   population-level output, i.e. results aggregated at population level (e.g. prevalence)

    -   individual-level output (optional), i.e. results tracked over time for each individuals

### **Author**

Veronica Malizia^1^

### **Supervision**

Sake J. de Vlas^2^, Federica Giardina^1^

^1^ Radboud University Medical Center, Department IQ Health, Biostatistics Research Group, Nijmegen, The Netherlands

^2^ Department of Public Health, Erasmus MC, University Medical Center Rotterdam, Rotterdam, The Netherlands
