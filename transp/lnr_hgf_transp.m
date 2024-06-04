function [pvec, pstruct] = lnr_hgf_transp(r, ptrans) 
% --------------------------------------------------------------------------------------------------
% lnr_hgf_transp: Function to perform transformation of parameters

% [Inputs]
% - r: array of responses
% - ptrans: structure that contains the parameters


% The structure and methodologies of this file are inspired
% from the HGF Toolbox, open source code available as part of the TAPAS
% software collection: Frässle, S., et al. (2021). TAPAS: An Open-Source Software Package 
% for Translational Neuromodeling and Computational Psychiatry. Frontiers in Psychiatry, 12:680811. 
% https://www.translationalneuromodeling.org/tapas

% This script is inspired from  Frässle, S., Aponte, E.A., Bollmann, S., Brodersen, K.H., Do, C.T., Harrison, O.K., Harrison, S.J., Heinzle, J., Iglesias, S., Kasper, L., Lomakina, E.I., Mathys, C., Müller-Schrader, M., Pereira, I., Petzschner, F.H., Raman, S., Schöbi, D., Toussaint, B., Weber, L.A., Yao, Y., Stephan, K.E.: 
% TAPAS: An Open-Source Software Package for Translational Neuromodeling and Computational Psychiatry, 
% Frontiers in Psychiatry 12, 857, 2021.

% extracting times
y = r.y(:,1);
y(r.irr) = [];

% Initialize an empty array to store the transformed values
pvec    = NaN(1,length(ptrans));
pstruct = struct;

% Mantain "bv" in the same space
pvec(1)      = (ptrans(1));
pstruct.bvmu   = pvec(1);

% Mantain "bi" in the same space
pvec(2)      = (ptrans(2));
pstruct.bimu   = pvec(2);

% Mantain "b1" in the same space
pvec(3)      = (ptrans(3));
pstruct.b1mu   = pvec(3);

% Apply exponential transformation 
pvec(4)      = exp(ptrans(4)); 
pstruct.sigmamu = pvec(4);

% Transform "T" using a sigmoid centered at the minimum value of y 
pvec(5)      = min(y)/(1+exp(-ptrans(5)));
pstruct.Tmu = pvec(5);

return;