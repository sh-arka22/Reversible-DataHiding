function [x]=HaarInv(ss, ds, sd, dd)

m=2*size(ss,1);
n=2*size(ss,2);
for j=1:n/2
    [rs(:,j)]=invhaar(ss(:,j),sd(:,j));
    [rd(:,j)]=invhaar(ds(:,j),dd(:,j));
end

%% Inverse Row Processing
for i=1:m
    [x(i,:)]=invhaar(rs(i,:),rd(i,:));
end

end
