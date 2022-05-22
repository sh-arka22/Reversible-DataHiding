function Y=overflowCompensation(V,indx,p,q)
% Compensating or converting the cover image histogram into its original 
% position. This technique is also called as histogram streching.
% Here according to the indx value we can detect the modified pixels 
% and original pixels in the range [p 2p] and [255-2*p 255-p].
[m,n]=size(V);
% Y=V;

if indx==0
    Y=V;
else
    t=1;
    for i=1:m
        for j=1:n
            if (abs(p)<=V(i,j))&& (V(i,j)<=255-q)
                Y(i,j)=V(i,j);
            elseif (0<=V(i,j))&&(V(i,j)<abs(p))
                if indx(t)==1
                    Y(i,j)=V(i,j)+p;
                    t=t+1;
                else
                    Y(i,j)=V(i,j);
                    t=t+1;
                end
            elseif (255-q<V(i,j))&&(V(i,j)<=255)
                if indx(t)==1
                    Y(i,j)=V(i,j)+q;
                    t=t+1;
                else
                    Y(i,j)=V(i,j);
                    t=t+1;
                end
            end
        end
    end
end
end


