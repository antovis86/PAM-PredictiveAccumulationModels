function c = RDM_hgf_config
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Contains the configuration for the Racing Diffusion Model according to:
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

% The default values for all parameters, except non-decision time (T), 
% are pre-defined with the mean set to "0" and the standard deviation set to "6.25".

% The default values of the parameter "T" are set to "-Inf" for the mean
% and "0" for the standard deviation. With this combination of parameters
% the value of "T" is not estimated. It is possible to estimate "T"
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
c.model = 'inverse gaussian: hgf';

% Intercept of B
c.b0Bmu = 0;
c.b0Bsa = 4;

% Muhat Slope of B
c.b1Bmu = 0;
c.b1Bsa = 4;

% Drift V for valid response (resp == input)
c.Vvmu = 0;
c.Vvsa = 4;

% Drift V for invalid response (resp != input)
c.Vimu = 0;
c.Visa = 4;

% Influence of muhat on v
c.bvmu = 0;
c.bvsa = 4;

% Non decision time
c.Tmu = 0;
c.Tsa = 4;


% Gather prior settings in vectors
c.priormus = [
    c.b0Bmu,...
    c.b1Bmu,...
    c.Vvmu,...
    c.Vimu,...
    c.bvmu,...
    c.Tmu,...
         ];

c.priorsas = [
    c.b0Bsa,...
    c.b1Bsa,...
    c.Vvsa,...
    c.Visa,...
    c.bvsa,...
    c.Tsa,...
     ];

% Model filehandle
c.obs_fun = @RDM_hgf;

% Handle to function that transforms observation parameters to their native space
% from the space they are estimated in
c.transp_obs_fun = @RDM_hgf_transp;

return;
