function z = interpolate(xb,zb,x)
%return z = z(x) from the linear interpolation of zb(xb) 

%find the two nearest points in xb to zb:
if x > min(xb) && x < max(xb)
    idx = find((x - xb)>0,1, 'last'); %return the index of largest entry of xb smaller than x
elseif x > max(xb)
    idx = length(xb)-1;
else
    idx = 1;
end

    
x_lo = xb(idx);
x_hi = xb(idx + 1);
z_lo = zb(idx);
z_hi = zb(idx + 1);

z = (z_hi - z_lo)/(x_hi - x_lo) *(x - x_lo) + z_lo;
