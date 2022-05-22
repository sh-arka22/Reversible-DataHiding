function [xdec,ss, ds, sd, dd]=HaarFwd(x)
% [xdec,ss, ds, sd, dd]=HaarFwd(x)
% Input:    x       of size     m*n
% Output:   xdec    of size     m*n
%           ss      of size     m/2*n/2
%           ds      of size     m/2*n/2
%           sd      of size     m/2*n/2
%           dd      of size     m/2*n/2
% 
% Also see HaarInv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % %    Ahmad Program    % % % %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = double(x);
[m, n]=size(x);


%% Forward Row Processing
for i=1:m
    [s(i,:), d(i,:)]=fwdhaar(x(i,:));
end



%% Forward Column processing
for j=1:n/2
    [ss(:,j), sd(:,j)]=fwdhaar(s(:,j));
    [ds(:,j), dd(:,j)]=fwdhaar(d(:,j));
end
xdec = [ss ds; sd dd];
end
