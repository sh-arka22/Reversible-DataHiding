function [z]=invhaar(s,d)
% function [rs, rd,z]=invcdf9by7(s,d,n)
[m,n]=size(s);
if n==1
    n=m;
end
rs0 = s;
rd0 = d;
% n=length(s); 
rs=rs0-floor(rd0/2);
rd=rs+rd0;
z(1:2:2*n)=rs;z(2:2:2*n)=rd;
if n==m
    z=z';
end
end
