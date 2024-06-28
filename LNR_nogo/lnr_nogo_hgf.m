function [logp, yhat, res] = lnr_nogo_hgf(r, infStates, ptrans)
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
bv = (ptrans(1));
bi = (ptrans(2));
b1 = (ptrans(3));
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
rt(rt==7012001)=NaN;

resp = r.y(:,2);
resp(r.irr)  = [];


% Fitting the non-decision time with the minimum value of estimated
% non-decision time
% Tmu = min(rt(not(isnan(rt))))/(1+exp(-ptrans(5)));
% rt(not(isnan(rt))) = max(eps,rt(not(isnan(rt)))-Tmu);


% extract the trial list and remove the irregular trials
u = r.u(:,1);
u(r.irr) = [];
% Calculate drift rate
mu_r1 = bv .* double(u == 1) + bi .* double(u == 0) + b1 .*(.5 -mu1hat);
mu_r2 = bv .* double(u == 0) + bi .* double(u == 1) + b1 .*(.5 - (1 -mu1hat));


% Calculate predicted log-likelihood
logp_reg = NaN(length(u),1);
for ntrial = 1:length(u)

    mu_pdf = (resp(ntrial) == 1)*mu_r1(ntrial) + (resp(ntrial) == 0) * mu_r2(ntrial);
    mu_cdf = (resp(ntrial) == 1)*mu_r2(ntrial) + (resp(ntrial) == 0) * mu_r1(ntrial);
    
    if not(isnan(rt(ntrial)))
        prob = utl_lnr_pdf(rt(ntrial),mu_pdf,mu_cdf,sigma);
        prob = max(prob, 1e-32);
        logp_reg(ntrial) = log(prob);
    elseif isnan(rt(ntrial))
        fun = @(x) utl_lnr_pdf(x,mu_pdf,mu_cdf,sigma);
        prob = integral(fun,0, 3);
        prob = max(prob, 1e-32);
        logp_reg(ntrial) = log(prob);

    end


    % Linear deterministic time distribution

end

reg = ~ismember(1:n,r.irr);
logp(reg) = logp_reg;

return;
