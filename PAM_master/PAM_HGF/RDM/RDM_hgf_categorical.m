function [logp, yhat, res] = RDM_hgf_categorical(r, infStates, ptrans)
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
mu1hat = squeeze(infStates(:,1,:,1));
mu1hat(r.irr,:) = []; mu1hat=mu1hat./sum(mu1hat,2);
n_choices = size(mu1hat,2);

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


% Calculate predicted log-likelihood
logp_reg = NaN(length(u),1);
for ntrial = 1:length(u)
    if rt(ntrial)>=0
        resp_yes = resp(ntrial);
        resp_not = setdiff(1:n_choices,resp_yes);
        % Extract the threshold for the specific trial
        a_pdf = a_a + b_a.*(mu1hat(ntrial,resp_yes) - 1/n_choices);
        a_cdf = a_a + b_a.*(mu1hat(ntrial,resp_not) - 1/n_choices);
        % Extract the drifts for the specific trial
        drift_pdf = a_v + b_val*double(u(ntrial) == resp_yes) + b_v*(mu1hat(ntrial,resp_yes) - 1/n_choices);
        drift_cdf = a_v + b_val*double(u(ntrial) == resp_not) + b_v* (mu1hat(ntrial,resp_not) - 1/n_choices);

        % Compute the defective distribution
        P = utl_inverse_gaussian_defective_categorical(rt(ntrial),drift_pdf,drift_cdf,a_pdf,a_cdf);

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
