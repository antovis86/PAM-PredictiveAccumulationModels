function [pvec, pstruct] = ddm_hgf_transp(r, ptrans) 
% --------------------------------------------------------------------------------------------------


% ddm_hgf_transp: Function to perform transformation of parameters

% [Inputs]
% - r: array of responses
% - ptrans: structure that contains the parameters


% The structure and methodologies of this file are inspired
% from the HGF Toolbox, open source code available as part of the TAPAS
% software collection: Frässle, S., et al. (2021). TAPAS: An Open-Source Software Package 
% for Translational Neuromodeling and Computational Psychiatry. Frontiers in Psychiatry, 12:680811. 
% https://www.translationalneuromodeling.org/tapas


% extracting times
try
    y = r.y(:,1);
    y(r.irr) = [];
catch
    y =nan;  % for "tapas_bayesian_parameter_average"
end

% Initialize an empty array to store the transformed values
pvec    = NaN(1,length(ptrans));
pstruct = struct;

% Apply exponential transformation to "a_a"
pvec(1)      = exp(ptrans(1));
pstruct.a_a   = pvec(1);

% Apply exponential transformation to "a_v"
pvec(2)      = exp(ptrans(2));
pstruct.a_v   = pvec(2);

% Transform "b_w" values in the range (-1,1)
pvec(3)      = 2/(1+exp(-ptrans(3)))-1;
pstruct.b_w   = pvec(3);

% Untrasformed "b_a" 
pvec(4)      = ptrans(4);
pstruct.b_a   = pvec(4);

% Untrasformed "b_v" 
pvec(5)      = ptrans(5);
pstruct.b_v   = pvec(5);

% Transform "Ter" using a sigmoid centered at the minimum value of y 
if ~isnan(y)
pvec(6)      = min(y)/(1+exp(-ptrans(6)));
pstruct.Ter = pvec(6);
else % for "tapas_bayesian_parameter_average"
pvec(6) = ptrans(6);
pstruct.Ter = pvec(6);
end

return;