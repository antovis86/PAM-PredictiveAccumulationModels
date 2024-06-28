function [pvec, pstruct] = RDM_nogo_hgf_transp(r, ptrans) 
% --------------------------------------------------------------------------------------------------

% RDM_hgf_transp: Function to perform transformation of parameters

% [Inputs]
% - r: array of responses
% - ptrans: structure that contains the parameters


% The structure and methodologies of this file are inspired
% from the HGF Toolbox, open source code available as part of the TAPAS
% software collection: Frässle, S., et al. (2021). TAPAS: An Open-Source Software Package 
% for Translational Neuromodeling and Computational Psychiatry. Frontiers in Psychiatry, 12:680811. 
% https://www.translationalneuromodeling.org/tapas

% extracting times
y = r.y(:,1);
y(r.irr) = [];

% Initialize an empty array to store the transformed values
pvec    = NaN(1,length(ptrans));
pstruct = struct;

% Apply exponential transformation to "B0"
pvec(1)      = exp(ptrans(1));
pstruct.b0B   = pvec(1);

% Transform "b1" values in the range (-1,1)
pvec(2)      = 2/(1+exp(-ptrans(2)))-1;
pstruct.b1B   = pvec(2);

% Apply exponential transformation to "Vv"
pvec(3)      = exp(ptrans(3));
pstruct.Vv   = pvec(3);

% Apply exponential transformation to "Vi"
pvec(4)      = exp(ptrans(4));
pstruct.Vi   = pvec(4);

% Transform "bv" values in the range (-1,1)
pvec(5)      = 2/(1+exp(-ptrans(5)))-1;
pstruct.bv = pvec(5);

% Transform "T" using a sigmoid centered at the minimum value of y 
pvec(6)      = min(y)/(1+exp(-ptrans(6)));
pstruct.Tmu = pvec(6);

return;