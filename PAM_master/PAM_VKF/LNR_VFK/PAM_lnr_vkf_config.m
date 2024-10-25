function [c] = PAM_lnr_vkf_config
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Contains the configuration for the PAM model combining the following models:
% PERCEPTUAL MODEL: Volatile Kalman Filter (https://doi.org/10.1371/journal.pcbi.1007963)
% DECISION MODEL: Log-normal race model (https://doi.org/10.3389/fpsyg.2012.00292) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Volatile Kalman filter parameters
c.lambda_x0 = 0;
c.lambda_fix = 1;
c.v0_x0 = log(0.5);
c.v0_fix = 1;
c.omega_x0 = 0;
c.omega_fix = 0;

% DDM PARAMETERS:
% Intercept of lognormal mu "a"
c.a_x0 = 0;
c.a_fix = 0;
% beta for balidity "b_val"
c.b_val_x0 = 0;
c.b_val_fix = 0;
% Muhat slope 
c.b_x0 = 0;
c.b_fix = 0;
% Lognormal sigma
c.sigma_x0 = 0;
c.sigma_fix = 0;
% Non decision time (fixed to 0; to estimate Ter, use the commented version)
c.Ter_x0 = -Inf; % c.Ter_x0 = 0
c.Ter_fix = 1; % c.Ter_fix = 0

% Gather initial settings in vectors
c.initpt = [
    c.lambda_x0,...
    c.v0_x0,...
    c.omega_x0,...
    c.a_x0,...
    c.b_val_x0,...
    c.b_x0,...
    c.sigma_x0,...
    c.Ter_x0,...
    ];

c.fix = logical([
    c.lambda_fix,...
    c.v0_fix,...
    c.omega_fix,...
    c.a_fix,...
    c.b_val_fix,...
    c.b_fix,...
    c.sigma_fix,...
    c.Ter_fix,...
    ]);

end