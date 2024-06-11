function f = utl_fsw(t, w, prec)
% Density at the lower barrier - one-parameter form
% https://doi.org/10.1016/j.jmp.2014.05.002.

K = utl_ks(t, w, prec);
f = zeros(1,length(t));
if(K > 0 && isfinite(K))
    for k = K:-1:1
        f = (w+2*k) * exp(-(w+2*k) * (w+2*k)/2./t) + ...
            (w-2*k) * exp(-(w-2*k) * (w-2*k)/2./t) + f;
    end
    f= (1./sqrt(2*pi.*t.*t.*t) .* (f + w .* exp(-w*w/2./t)));
end


end