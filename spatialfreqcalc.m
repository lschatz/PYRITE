function [spatialfreq,uniqk,totalvar]=spatialfreqcalc(in,xp,yp,cp,nmodes)

%% Calculate m and n spatial frequency indexing
m=xp-cp;
n=yp-cp;
%% Calculate k which are the spatial frequencies
for i=1:size(m)
    k(i)=sqrt(m(i).^2+n(i).^2); 
end
%% This fixes the indexing. 
%The m index runs from -nmodes to nmodes. 
%The problem is that matlab can't handle negatives and zeros in indexes
%So it maps the m indexs to start from 1 and be positive. so -nmodes now is
%1.
m=m+abs(min(m))+1;
n=n+abs(min(n))+1;

%% Sine and Cosine matrix
for i=1:size(m)
    sinematrix(m(i),n(i))=in(i*2);
    cosinematrix(m(i),n(i))=in(i*2-1);
    totalvar(i)=sqrt(sinematrix(m(i),n(i)).^2+cosinematrix(m(i),n(i)).^2);
end

%figure; imagesc(-nmodes:0,-nmodes:nmodes,quadsinematrix); title('sine matrix 4PWFS'); colorbar; caxis([0, max(max(quadvarspatialmode))]); 
%figure; imagesc(-nmodes:0,-nmodes:nmodes,quadcosinematrix); title('cosine matrix 4PWFS'); colorbar; caxis([0, max(max(quadvarspatialmode))]); 
modes=sqrt(sinematrix.^2+cosinematrix.^2);
figure; imagesc(-nmodes:0,-nmodes:nmodes, modes);title('Variance of fitted Fourier Amplitude'); colorbar; axis equal; caxis([0, max(max(modes))]); xlim([-nmodes-0.5 0.5]); xlabel('n Fourier index'); ylabel('m Fourier index')


%% Combine terms so that there is one value per unique spatial frequency

uniqk=unique(k); %find all the unique values of spatial frequency
ssz=size(uniqk); 
for i=1:ssz(2)
    I=find(k==uniqk(i)); %find indicies of unique values
    sz=size(I);
    for j=1:sz(2)
        values(j)=totalvar(I(j));
    end
    spatialfreq(i)=sqrt(sum(values.^2));
end


end
