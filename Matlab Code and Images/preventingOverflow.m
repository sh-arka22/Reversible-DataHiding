function [Y,indx,p,q] = preventingOverflow(V)
% Moving the lower range pixels [0 p] where p=alpha*k here k=number of 
% data embedding levels and alpha (constant) is the embedding depth 
% V=1:1:16; p=3; %For Testing
[m,n]=size(V);
% X=V;
% [indx1, indx2]=find(V<p);
% X(indx1,indx2)=V(indx1,indx2)+p;
% [indx3, indx4]=find(V>255-p);
% X(indx3, indx4)=V(indx3, indx4)-p;
q=max(max(V))-255;p=min(min(V));
if (p<0)||(q>0)
    indx=[];
    for i=1:m
        for j=1:n
            if (abs(p)<=V(i,j))&& (V(i,j)<=255-q)
            Y(i,j)=V(i,j);
            elseif V(i,j)<0
                Y(i,j)=V(i,j)-p;
                indx=[indx 1];
            elseif (p<=V(i,j))&&(V(i,j)<abs(p))
                Y(i,j)=V(i,j);
                indx=[indx 0];
            elseif V(i,j)>255
                Y(i,j)=V(i,j)-q;
                indx=[indx 1];
            elseif (255-q<V(i,j))&&(V(i,j)<=255)
                Y(i,j)=V(i,j);
                indx=[indx 0];
            end
        end
    end
else
    Y=V;
    indx=0;
end
            
end


