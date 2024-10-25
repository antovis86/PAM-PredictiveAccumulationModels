function c = RDM_hgf_categorical_config
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Contains the configuration for the Racing Diffusion Model for multiple choices (>=2)according to:
% Tillman, G., Van Zandt, T., & Logan, G. D. (2020). Sequential sampling models without random between-trial variability: the racing diffusion model of speeded decision making. 
% Psychonomic Bulletin & Review, 27(5), 911–936
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The structure and methodologies of this file are inspired
% from the HGF Toolbox, open source code available as part of the TAPAS
% software collection: Frässle, S., et al. (2021). TAPAS: An Open-Source Software Package 
% for Translational Neuromodeling and Computational Psychiatry. Frontiers in Psychiatry, 12:680811. 
% https://www.translationalneuromodeling.org/tapas
% 




% The RDM configuration consists of the priors of parameters and initial values. All priors are
% Gaussian in the space where the quantity they refer to is estimated. They are specified by their
% sufficient statistics: mean and variance (NOT standard deviation).

% The default values for all parameters, except non-decision time (Ter), 
% are pre-defined with the mean set to "0" and the standard deviation set to "6.25".

% The default values of the parameter "Ter" are set to "-Inf" for the mean
% and "0" for the standard deviation. With this combination of parameters
% the value of "Ter" is not estimated. It is possible to estimate "Ter"
% parameter by changing the values of the aforementioned parameters.

% It is possible to modify the specified values and assign custom values to each parameter.


% This config file does not take any input and returns the
% structure "c" as output. "c" contains all the necessary values to 
% define and build a Racing Diffusion Model.
%
% --------------------------------------------------------------------------------------------------
% 


% Config structure
c = struct;

% Model's name
c.model = 'racing diggusion model: hgf_categorical';

% Intercept of decision threshold "a"
c.a_amu = 0;
c.a_asa = 4;

% Muhat slope effect for "a"
c.b_amu = 0;
c.b_asa = 4;

% Intercept of drift rate "v"
c.a_vmu = 0;
c.a_vsa = 4;

% Effect of validity (resp = input) on the drift 
c.b_valmu = 0;
c.b_valsa = 4;

% Muhat slope effect for "v"
c.b_vmu = 0;
c.b_vsa = 4;

% Non decision time
c.Termu = -Inf;
c.Tersa = 0;


% Gather prior settings in vectors
c.priormus = [
    c.a_amu,...
    c.b_amu,...
    c.a_vmu,...
    c.b_valmu,...
    c.b_vmu,...
    c.Termu,...
         ];

c.priorsas = [
    c.a_asa,...
    c.b_asa,...
    c.a_vsa,...
    c.b_valsa,...
    c.b_vsa,...
    c.Tersa,...
     ];

% Model filehandle
c.obs_fun = @RDM_hgf_categorical;

% Handle to function that transforms observation parameters to their native space
% from the space they are estimated in
c.transp_obs_fun = @RDM_hgf_categorical_transp;

return;
