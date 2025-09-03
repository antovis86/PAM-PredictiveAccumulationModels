# PAM 

---

This repository contains the MATLAB implementation of PAM framework. The structure and methodologies are inspired from the HGF toolbox, open source code available as part of the TAPAS software collection. For the Python implementation see: https://github.com/francesco-cal98/Pam---Predictive-Accumulator-models
For the R implementation see: https://github.com/Mar-Cald/PAM-PredictiveAccumulationModels
References:

TAPAS - Frässle, S., et al. (2021). TAPAS: An Open-Source Software Package for
Translational Neuromodeling and Computational Psychiatry. Frontiers in
Psychiatry, 12:680811. https://doi.org/10.3389/fpsyt.2021.680811

PAM - https://arxiv.org/abs/2411.13203

The repository is structured in the following way:

```
project-root
├── PAM_master
│   ├── PAM_HGF         # implementation of PAM framework using HGF model
│   ├── PAM_VKF         # implementation of PAM framework using VKF model
│   └── Examples
│       ├── HGF_examples  # usage examples with HGF model
│       └── VKF_examples  # usage examples with VKF model
```


## Dependencies 

| Required | Package           | Remarks         |
| ---------|-------------------|-----------------|
| Yes      | [MATLAB]          |                 |
| Yes      | [TAPAS]           | 4 PAM_HGF       |
| Yes      | [VKF]             | 4 PAM_VKF       |

----

1. Download/clone the latest release and unzip it.
2. Open MATLAB
3. Navigate to the project directory
```
cd '/path/to/PAM'
```
4. Add the project folder and its subfolders to the MATLAB path:
```
addpath(genpath('/path/to/PAM'))
```
5. Run the scripts 
