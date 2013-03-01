% implementation of the function that is provided by wp2

function f = wp2inp(x)

x = x.';
f = exp(exp(-2*x)*25.*(1+2*x)).*(abs(x)<=1/2);
f = f/max(f);