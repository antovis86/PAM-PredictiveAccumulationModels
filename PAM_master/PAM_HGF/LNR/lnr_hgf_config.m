function c = lnr_hgf_config
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Contains the configuration for the Linear Deterministic Accumulator (LNR) according to:
% Heathcote, A., & Love, J. (2012). Linear Deterministic Accumulator Models of Simple Choice. Frontiers in Psychology, 3.

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The structure and methodologies of this file are inspired
% from the HGF Toolbox, open source code available as part of the TAPAS
% software collection: Frässle, S., et al. (2021). TAPAS: An Open-Source Software Package 
% for Translational Neuromodeling and Computational Psychiatry. Frontiers in Psychiatry, 12:680811. 
% https://www.translationalneuromodeling.org/tapas
% 

% The LNR configuration consists of the priors of parameters and initial values. All priors are
% Gaussian in the space where the quantity they refer to is estimated. They are specified by their
% sufficient statistics: mean and variance (NOT standard deviation).

% The default values for all parameters, except non-decision time (Ter), 
% are pre-defined with the mean set to "0" and the standard deviation set to "4".

% The default values of the parameter "Ter" are set to "-Inf" for the mean
% and "0" for the standard deviation. With this combination of parameters
% the value of "Ter" is not estimated. It is possible to estimate "Ter"
% parameter by changing the values of the aforementioned parameters.

% It is possible to modify the specified values and assign custom values to each parameter.

% This config file does not take any input and returns the
% structure "c" as output. "c" contains all the necessary values to 
% define and build a Linear Deterministic Model.
%
% --------------------------------------------------------------------------------------------------


% Config structure
c = struct;

% Model name
c.model = 'lnr: hgf';

% intercept of lognormal mu (resp == input)
c.amu = 0;
c.asa = 4;

% beta for validity (resp = input)
c.b_valmu = 0;
c.b_valsa = 4;

% Beta for muhat
c.bmu = 0;
c.bsa = 4;

% Standard deviation
c.sigmamu = 0;
c.sigmasa = 4;

% Non decision time
c.Termu = -Inf;
c.Tersa = 0;


% Gather prior settings in vectors
c.priormus = [
    c.amu,...
    c.b_valmu,...
    c.bmu,...
    c.sigmamu,...
    c.Termu,...
    ];

c.priorsas = [
    c.asa,...
    c.b_valsa,...
    c.bsa,...
    c.sigmasa,...
    c.Tersa,...
    ];

% Model filehandle
c.obs_fun = @lnr_hgf;

% Handle to function that transforms observation parameters to their native space
% from the space they are estimated in
c.transp_obs_fun = @lnr_hgf_transp;

return;
