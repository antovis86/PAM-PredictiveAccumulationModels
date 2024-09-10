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

% The default values for all parameters, except non-decision time (T), 
% are pre-defined with the mean set to "0" and the standard deviation set to "6.25".

% The default values of the parameter "T" are set to "-Inf" for the mean
% and "0" for the standard deviation. With this combination of parameters
% the value of "T" is not estimated. It is possible to estimate "T"
% parameter by changing the values of the aforementioned parameters.

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

% Upper absorbing barriers
c.amu = 0;
c.asa = 4;

% Drift
c.vmu = 0;
c.vsa = 4;

% Influence of muhat on the relative starting point (w)
c.bwmu = 0;
c.bwsa = 4;

% Influence of muhat on a
c.bamu = 0;
c.basa = 4;

% Influence of muhat on v
c.bvmu = 0;
c.bvsa = 4;

% Non decision time
c.Tmu = 0;
c.Tsa = 4;


% Gather prior settings in vectors
c.priormus = [
    c.amu,...
    c.vmu,...
    c.bwmu,...
    c.bamu,...
    c.bvmu,...
    c.Tmu,...
    ];

c.priorsas = [
    c.asa,...
    c.vsa,...
    c.bwsa,...
    c.basa,...
    c.bvsa,...
    c.Tsa,...
    ];

% Model filehandle
c.obs_fun = @ddm_hgf;

% Handle to function that transforms observation parameters to their native space
% from the space they are estimated in
c.transp_obs_fun = @ddm_hgf_transp;

return;
