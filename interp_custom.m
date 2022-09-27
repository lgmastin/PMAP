function yq = interp_custom(x,y,xq)
% linear 1D interpolation
% holds end values for extrapolation

yq = interp1(x,y,xq,'linear');

if xq<x(1)
    yq = y(1);
elseif xq>x(end)
    yq = y(end);
end

end