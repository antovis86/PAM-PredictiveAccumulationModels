function K = utl_ks(t, w, prec)
% Number of terms for the density
% https://doi.org/10.1016/j.jmp.2014.05.002.

K1 = (sqrt(2*t) - w)/2;
K2 = K1;
	u_eps = min(-1, log(2.*pi.*t.*t.*prec.*prec)); 
	arg = -t .* (u_eps - sqrt(-2.*u_eps - 2));
	K2(arg > 0) = 1/2 .* sqrt(arg) - w/2;
	K = ceil(max([K1 K2]));
end