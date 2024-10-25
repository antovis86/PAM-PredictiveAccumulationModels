function [c] = PAM_rdm_vkf_config
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Contains the configuration for the PAM model combining the following models:
% PERCEPTUAL MODEL: Volatile Kalman Filter (https://doi.org/10.1371/journal.pcbi.1007963)
% DECISION MODEL: Racing Diffusion Model (https://doi.org/10.3758/s13423-020-01719-6)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Volatile Kalman filter parameters
c.lambda_x0 = 0;
c.lambda_fix = 1;
c.v0_x0 = log(0.5);
c.v0_fix = 1;
c.omega_x0 = 0;
c.omega_fix = 0;

% RDM PARAMETERS:
% Intercept of decision threshold "a"
c.a_a_x0 = 0;
c.a_a_fix = 0;
% Muhat slope effect for "a"
c.b_a_x0 = 0;
c.b_a_fix = 0;
% Intercept of drift rate "v"
c.a_v_x0 = 0;
c.a_v_fix = 0;
% Effect of validity (resp = input) on the drift 
c.b_val_x0 = 0;
c.b_val_fix = 0;
% Muhat slope effect for "v"
c.b_v_x0 = 0;
c.b_v_fix = 0;
% Non decision time (fixed to 0; to estimate Ter, use the commented version)
c.Ter_x0 = -Inf; % c.Ter_x0 = 0
c.Ter_fix = 1; % c.Ter_fix = 0

% Gather initial settings in vectors
c.initpt = [
    c.lambda_x0,...
    c.v0_x0,...
    c.omega_x0,...
    c.a_a_x0,...
    c.b_a_x0,...
    c.a_v_x0,...
    c.b_val_x0,...
    c.b_v_x0,...
    c.Ter_x0,...
    ];

c.fix = logical([
    c.lambda_fix,...
    c.v0_fix,...
    c.omega_fix,...
    c.a_a_fix,...
    c.b_a_fix,...
    c.a_v_fix,...
    c.b_val_fix,...
    c.b_v_fix,...
    c.Ter_fix,...
    ]);

end