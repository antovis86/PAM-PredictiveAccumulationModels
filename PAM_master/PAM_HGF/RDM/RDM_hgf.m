function [logp, yhat, res] = RDM_hgf(r, infStates, ptrans)
% RDM_hgf: Calculates the log-probability of response speed y (in units of seconds)
% according to the Racing Diffusion Model:
% Tillman, G., Van Zandt, T., & Logan, G. D. (2020). Sequential sampling models without random between-trial variability: the racing diffusion model of speeded decision making. 
% Psychonomic Bulletin & Review, 27(5), 911–936
%

% [Inputs]
% - r: array of responses
% - infStates: array that contains the parameters of the HGF
% - ptrans: structure that contains the parameters for the response model


% The structure and methodologies of this file are inspired
% from the HGF Toolbox, open source code available as part of the TAPAS
% software collection: Frässle, S., et al. (2021). TAPAS: An Open-Source Software Package 
% for Translational Neuromodeling and Computational Psychiatry. Frontiers in Psychiatry, 12:680811. 
% https://www.translationalneuromodeling.org/tapas
%
% --------------------------------------------------------------------------------------------------

% Transform parameters to their native space
a_a = exp(ptrans(1));
b_a = ptrans(2);
a_v = exp(ptrans(3));
b_val = exp(ptrans(4));
b_v = ptrans(5);

% Initialize returned log-probabilities, predictions,
% and residuals as NaNs so that NaN is returned for all
% irregualar trials
n = size(infStates,1);
logp = NaN(n,1);
yhat = NaN(n,1); % not used
res  = NaN(n,1); % not used

% Weed irregular trials out from inferred states, responses, and inputs
mu1hat = infStates(:,1,1);
mu1hat(r.irr) = [];

% Extract the response times from "y" array
rt = r.y(:,1);
rt(r.irr) = [];

% Extract the response from "y" array
resp = r.y(:,2);
resp(r.irr)  = [];

% Fitting the non-decision time with the minimum value of estimated
% non-decision time
Ter = min(rt)/(1+exp(-ptrans(6)));
rt = max(eps,rt-Ter);


% extract the trial list and remove the irregular trials
u = r.u(:,1);
u(r.irr) = [];

% Calculate trial-wise threshold for both the accumulators
a_c1 = a_a + b_a.*(mu1hat-.5);
a_c0 = a_a + b_a.*((1-mu1hat)-.5);

% Calculate drift rate for both accumulators
drift_c1 = a_v + b_val.*(u==1) +  b_v .* (mu1hat - .5);
drift_c0 = a_v + b_val.*(u==0) +  b_v .* ((1-mu1hat) - .5);

% Calculate predicted log-likelihood
logp_reg = NaN(length(u),1);
for ntrial = 1:length(u)
    if rt(ntrial)>=0
        
        % Extract the drifts for the specific trial
        drift_pdf = (resp(ntrial) == 1)*drift_c1(ntrial) + (resp(ntrial) == 0) * drift_c0(ntrial);
        drift_cdf = (resp(ntrial) == 1)*drift_c0(ntrial) + (resp(ntrial) == 0) * drift_c1(ntrial);

        % Extract the threshold for the specific trial
        a_pdf = (resp(ntrial) == 1)*a_c1(ntrial) + (resp(ntrial) == 0) * a_c0(ntrial);
        a_cdf = (resp(ntrial) == 1)*a_c0(ntrial) + (resp(ntrial) == 0) * a_c1(ntrial);
        
        % Compute the defective distribution
        P = utl_inverse_gaussian_defective(rt(ntrial),drift_pdf,drift_cdf,a_pdf,a_cdf);

        if P>0
            logp_reg(ntrial) = log(P+eps);
        else
            logp_reg(ntrial) = NaN;
        end
    else
        logp_reg(ntrial) = NaN;
    end
end

reg = ~ismember(1:n,r.irr);
logp(reg) = logp_reg;
return;
