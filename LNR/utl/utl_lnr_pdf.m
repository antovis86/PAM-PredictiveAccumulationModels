function [probs] = utl_lnr_pdf(x,mu1,mu2,sigma)


pdf = lognpdf(x,mu1,sigma);
survival = logncdf(x,mu2,sigma,'upper');

probs = pdf.*survival;