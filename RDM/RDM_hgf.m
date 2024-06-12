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
b0B = exp(ptrans(1));
b1B = 2/(1+exp(-ptrans(2)))-1;
Vv = exp(ptrans(3));
Vi = exp(ptrans(4));
bv = 2/(1+exp(-ptrans(5)))-1;

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
T = min(rt)/(1+exp(-ptrans(6)));
rt = max(eps,rt-T);


% extract the trial list and remove the irregular trials
u = r.u(:,1);
u(r.irr) = [];

% Calculate threshold for both the accumulators
B_dx = b0B + b1B.*(.5-mu1hat).*2.*b0B;
B_sx = b0B + b1B.*(.5-(1-mu1hat)).*2*b0B;


% Calculate drift rate for both accumulators
drift_r1 = Vv .* double(u == 1) + Vi .* double(u == 0);
drift_r1 = drift_r1 + (bv .* (mu1hat - .5).*2.*drift_r1);

drift_r2 = Vv .* double(u == 0) + Vi .* double(u == 1);
drift_r2 = drift_r2 + (bv.*((1-mu1hat) - .5).*2.*drift_r2);


% Calculate predicted log-likelihood
logp_reg = NaN(length(u),1);
for ntrial = 1:length(u)
    if rt(ntrial)>=0
        
        % Extract the drifts for the specific trial
        drift_pdf = (resp(ntrial) == 1)*drift_r1(ntrial) + (resp(ntrial) == 0) * drift_r2(ntrial);
        drift_cdf = (resp(ntrial) == 1)*drift_r2(ntrial) + (resp(ntrial) == 0) * drift_r1(ntrial);

        % Extract the threshold for the specific trial
        B_pdf = (resp(ntrial) == 1)*B_dx(ntrial) + (resp(ntrial) == 0) * B_sx(ntrial);
        B_cdf = (resp(ntrial) == 1)*B_sx(ntrial) + (resp(ntrial) == 0) * B_dx(ntrial);
        
        % Compute the defective distribution
        P = utl_inverse_gaussian_defective(rt(ntrial),drift_pdf,drift_cdf,B_pdf,B_cdf);

        % Correct if P<0
        P(P<0)=0;

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
