function q = percentile_simple(x,p)
%PERCENTILE_SIMPLE Toolbox-free percentile with linear interpolation.
x=sort(x(:));
if isempty(x), q=NaN; return; end
if numel(x)==1, q=x; return; end
r=1+(numel(x)-1)*p/100;
i=floor(r); j=ceil(r);
if i==j,q=x(i);else,q=x(i)+(r-i)*(x(j)-x(i));end
end
