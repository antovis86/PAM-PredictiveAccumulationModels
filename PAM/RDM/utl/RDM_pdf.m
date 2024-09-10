function [pdf] = RDM_pdf(x,drift_pdf,threshold)


pdf = (threshold * 1000) ./ sqrt(2*pi*(x.^3)) .* exp(-0.5 * ((drift_pdf*(x/1000) - threshold).^2) ./ (x/1000));
