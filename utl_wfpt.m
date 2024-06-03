function p = utl_wfpt(t, v, a, w,prec)
% First passage time for Wiener diffusion model
% by Gondan, Blurton, and Kesselmeier (2014)
% https://doi.org/10.1016/j.jmp.2014.05.002.
%
% USAGE: p = wfpt(t,v,a,w,prec)
%
% INPUTS:
%   t - hitting time (e.g., response time in milliseconds)
%   v - drift rate
%   a - threshold
%   w - bias (default: 0.5)
%   err - error threshold (default: 1e-4)
%
% OUTPUTS:
%   p - probability density at the lower barrier

if nargin<4; w = .5; end
if nargin<5; prec = 1e-4; end

p = 1/a/a .* exp(-v*a*w - v*v.*t/2);
p = p .* utl_fsw(t/a/a, w, prec./p);

end