function [c] = PAM_ddm_vkf_config
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Contains the configuration for the PAM model combining the following models:
% PERCEPTUAL MODEL: Volatile Kalman Filter (https://doi.org/10.1371/journal.pcbi.1007963)
% DECISION MODEL: Wiener Diffusion Model (https://doi.org/10.1016/j.jmp.2009.02.003)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Volatile Kalman filter parameters
c.lambda_x0 = 0;
c.lambda_fix = 1;
c.v0_x0 = log(0.5);
c.v0_fix = 1;
c.omega_x0 = 0;
c.omega_fix = 0;

% DDM PARAMETERS:
% Intercept of boundary separation "a"
c.a_a_x0 = 0;
c.a_a_fix = 0;
% Intercept of drift rate "v"
c.a_v_x0 = 0;
c.a_v_fix = 0;
% Muhat slope effect for "w"
c.b_w_x0 = 0;
c.b_w_fix = 0;
% Muhat slope effect for "a"
c.b_a_x0 = 0;
c.b_a_fix = 0;
% Muhat slope effect for "v"
c.b_v_x0 = 0;
c.b_v_fix = 0;
% Non decision time 
c.Ter_x0 = 0; 
c.Ter_fix = 0;

% Gather initial settings in vectors
c.initpt = [
    c.lambda_x0,...
    c.v0_x0,...
    c.omega_x0,...
    c.a_a_x0,...
    c.a_v_x0,...
    c.b_w_x0,...
    c.b_a_x0,...
    c.b_v_x0,...
    c.Ter_x0,...
    ];

c.fix = logical([
    c.lambda_fix,...
    c.v0_fix,...
    c.omega_fix,...
    c.a_a_fix,...
    c.a_v_fix,...
    c.b_w_fix,...
    c.b_a_fix,...
    c.b_v_fix,...
    c.Ter_fix,...
    ]);

end