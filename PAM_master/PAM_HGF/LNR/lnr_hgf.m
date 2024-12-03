function [logp, yhat, res] = lnr_hgf(r, infStates, ptrans)
% Calculates the log-probability of response speed y (in units of ms^-1) 
% 
% lnr_hgf: Calculates the log-probability of response speed y (in units of seconds)
% according to the Linear Deterministic Model:
% Tillman, G., Van Zandt, T., & Logan, G. D. (2020). Sequential sampling models without random between-trial variability: the racing diffusion model of speeded decision making. 
% Psychonomic Bulletin & Review, 27(5), 911–936


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
a = (ptrans(1));
b_val = (ptrans(2)); 
b = (ptrans(3));
sigma = exp(ptrans(4));

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

rt = r.y(:,1);
rt(r.irr) = [];

resp = r.y(:,2);
resp(r.irr)  = [];

% Fitting the non-decision time with the minimum value of estimated
% non-decision time
Ter = min(rt)/(1+exp(-ptrans(5)));
rt = max(eps,rt-Ter);


% extract the trial list and remove the irregular trials
u = r.u(:,1);
u(r.irr) = [];


% Calculate trial-wise drift rates for the two accumulators 
mu_c1 = a + b_val.* double(u == 1) + b .*(mu1hat - .5);
mu_c0 = a + b_val.* double(u == 0) + b .*((1 -mu1hat) - .5);


% Calculate predicted log-likelihood
logp_reg = NaN(length(u),1);
for ntrial = 1:length(u)
    if rt(ntrial)>0
        mu_pdf = (resp(ntrial) == 1)*mu_c1(ntrial) + (resp(ntrial) == 0) * mu_c0(ntrial);
        mu_cdf = (resp(ntrial) == 1)*mu_c0(ntrial) + (resp(ntrial) == 0) * mu_c1(ntrial);
        P = utl_lnr_pdf(rt(ntrial),mu_pdf,mu_cdf,sigma);
        
        % Correct if P<0 
        P(P<0)=0;
   
        % Linear deterministic time distribution
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
