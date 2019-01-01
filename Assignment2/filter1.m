
function [I_patch] = filter1(i1,sd)
s=0;
sig=1.5*sd;
for i=-2:2
    for j=-2:2
    K(i+3,j+3) = exp(-(i.^2+j.^2)/(2*sig.^2));
    K(i+3,j+3) = K(i+3,j+3)./(2*pi*(sig.^2));
    s=s+K(i+3,j+3); 
    end
end
K = K./s;

%convolving image with gaussian filter

[r,c]=size(i1);
I_patch = zeros(r,c);
s1=0;
for i = 1:r-4
    for j = 1:c-4
        sum1 = realmin;
        for k = 1:5
           for l = 1:5
               sum1=sum1+i1(i+k-1,j+l-1).*K(k,l);
               s1=s1+sum1;
           end
        end
        I_patch(i+3,j+3)= sum1;
    end
end
if(int32(s1)~=0)
I_patch=I_patch./double(s1);
end
end

