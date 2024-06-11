function [pdf] = RDM_pdf(x,drift_pdf,threshold)

pdf = threshold./ sqrt(2*pi*x.^3) .* exp(-0.5 * ((drift_pdf*x - threshold).^2) ./ x);

