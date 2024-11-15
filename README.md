down# PAM 

---

This repository contains the MATLAB implementation of PAM framework. The structure and methodologies are inspired from the HGF toolbox, open source code available as part of the TAPAS software collection
Original reference:

Frässle, S., et al. (2021). TAPAS: An Open-Source Software Package for
Translational Neuromodeling and Computational Psychiatry. Frontiers in
Psychiatry, 12:680811. https://doi.org/10.3389/fpsyt.2021.680811

The repository is structured in the following way:

\```
project-root
├── PAM_master
│   ├── PAM_HGF: contains the implementation of PAM framework using HGF model as perceptual model.
│   ├── PAM_VKF: contains the implementation of PAM framework using Volatile Kalman Filter (VKF) model as perceptual model.
│
└── Examples
    ├── HGF_examples: contains different usage examples with HGF model as perceptual model.
    ├── VKF_examples: contains different usage examples with VKF model as perceptual model.
\```

---- 

## Dependencies 

| Required | Package           | Remarks         |
| ---------|-------------------|-----------------|
| Yes      | [MATLAB]          |                 |
| Yes      | [TAPAS]           |                 |

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
