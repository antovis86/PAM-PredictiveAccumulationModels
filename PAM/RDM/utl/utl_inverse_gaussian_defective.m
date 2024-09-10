function [probs] = utl_inverse_gaussian_defective(x,drift_pdf,drift_cdf,threshold1,threshold2)

% Convert milliseconds to seconds within the formulas
pdf = (threshold1 * 1000) ./ sqrt(2*pi*(x.^3)) .* exp(-0.5 * ((drift_pdf.*(x/1000) - threshold1).^2) ./ (x/1000));

cdf = normcdf(((drift_cdf.*(x/1000)) - threshold2) ./ sqrt(x/1000)) + ...
    exp(2.*drift_cdf.*threshold2) .* normcdf((-(drift_cdf.*(x/1000)) - threshold2) ./ sqrt(x/1000));

probs = pdf.*(1-cdf);