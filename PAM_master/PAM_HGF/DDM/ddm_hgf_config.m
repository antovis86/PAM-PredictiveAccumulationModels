function c = ddm_hgf_config
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Contains the configuration for First passage time for Wiener diffusion model (DDM) according to:
% by Gondan, Blurton, and Kesselmeier (2014)
% https://doi.org/10.1016/j.jmp.2014.05.002.
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The structure and methodologies of this file are inspired
% from the HGF Toolbox, open source code available as part of the TAPAS
% software collection: Frässle, S., et al. (2021). TAPAS: An Open-Source Software Package 
% for Translational Neuromodeling and Computational Psychiatry. Frontiers in Psychiatry, 12:680811. 
% https://www.translationalneuromodeling.org/tapas



% The DDM configuration consists of the priors of parameters and initial values. All priors are
% Gaussian in the space where the quantity they refer to is estimated. They are specified by their
% sufficient statistics: mean and variance (NOT standard deviation).

% The default values for all parameters are pre-defined with the mean set to "0" and the standard deviation set to "4".

% It is possible to modify the specified values and assign custom values to each parameter.
% This config file does not take any input and returns the
% structure "c" as output. "c" contains all the necessary values to 
% define and build a Wiener Diffusion Model.
%
% --------------------------------------------------------------------------------------------------
% 

% Config structure
c = struct;

% Model name
c.model = 'ddm: hgf';

% Intercept of boundary separation (a)
c.a_amu = 0;
c.a_asa = 4;

% Intercept of drift rate (v)
c.a_vmu = 0;
c.a_vsa = 4;

% Muhat slope effect for starting point (w)
c.b_wmu = 0;
c.b_wsa = 4;

% Muhat slope effect for "a"
c.b_amu = 0;
c.b_asa = 4;

% Muhat slope effect for "v"
c.b_vmu = 0;
c.b_vsa = 4;

% Non decision time
c.Termu = 0;
c.Tersa = 4;


% Gather prior settings in vectors
c.priormus = [
    c.a_amu,...
    c.a_vmu,...
    c.b_wmu,...
    c.b_amu,...
    c.b_vmu,...
    c.Termu,...
    ];

c.priorsas = [
    c.a_asa,...
    c.a_vsa,...
    c.b_wsa,...
    c.b_asa,...
    c.b_vsa,...
    c.Tersa,...
    ];

% Model filehandle
c.obs_fun = @ddm_hgf;

% Handle to function that transforms observation parameters to their native space
% from the space they are estimated in
c.transp_obs_fun = @ddm_hgf_transp;

return;
