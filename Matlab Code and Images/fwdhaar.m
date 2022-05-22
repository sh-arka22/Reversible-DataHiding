function [s, d]=fwdhaar(x)
[m, n] =size(x);
if n==1
    n=m;
end
s0 = x(1:2:n);
d0 = x(2:2:n);
d= d0-s0;
s=s0+floor(d/2);
% if n==m
%     d= (d0-s0)';
%     s=(s0+floor(d/2))';
%     
%     d= (s0-d0)';
%     s=(floor((d0+s0)/2))';
% else
%     d= (s0-d0);
%     s=(floor((d0+s0)/2));
% end
end
