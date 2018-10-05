function [uniqk, trispatialfreq, quadspatialfreq]=fbsvarcalc(trival, quadval,xp,yp,cp,npix,nmodes)
%REVISIT THIS
% Get rid of NaN values
trival(isnan(trival))=0;
quadval(isnan(quadval))=0;

%% Calculate variances
% trivar=var(trival,0,2);
% quadvar=var(quadval,0,2);
trivar=rms(trival,2);
quadvar=rms(quadval,2);


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

%% Quad Sine and Cosine matrix
for i=1:size(m)
    quadsinematrix(m(i),n(i))=quadvar(i*2);
    quadcosinematrix(m(i),n(i))=quadvar(i*2-1);
    totalquadvar(i)=sqrt(quadsinematrix(m(i),n(i)).^2+quadcosinematrix(m(i),n(i)).^2);
end

%figure; imagesc(-nmodes:0,-nmodes:nmodes,quadsinematrix); title('sine matrix 4PWFS'); colorbar; caxis([0, max(max(quadvarspatialmode))]); 
%figure; imagesc(-nmodes:0,-nmodes:nmodes,quadcosinematrix); title('cosine matrix 4PWFS'); colorbar; caxis([0, max(max(quadvarspatialmode))]); 

quadmodes=sqrt(quadsinematrix.^2+quadcosinematrix.^2);
figure; imagesc(-nmodes:0,-nmodes:nmodes,quadmodes);title('4PWFS Variance of fitted Fourier Amplitude'); colorbar; axis equal; caxis([0, max(max(quadmodes))]); xlim([-nmodes-0.5 0.5]); xlabel('n Fourier index'); ylabel('m Fourier index')

%% Tri Sine and Cosine matrix
for i=1:size(m)
    trisinematrix(m(i),n(i))=trivar(i*2);
    tricosinematrix(m(i),n(i))=trivar(i*2-1);
    totaltrivar(i)=sqrt(trisinematrix(m(i),n(i)).^2+tricosinematrix(m(i),n(i)).^2);
end

%figure; imagesc(sinematrix); title('sine matrix'); colorbar; caxis([0, max(max(trivarspatialmode))]); axis equal
%figure; imagesc(cosinematrix); title('cosine matrix'); colorbar; caxis([0, max(max(trivarspatialmode))]); axis equal

trimodes=sqrt(trisinematrix.^2+tricosinematrix.^2);
figure; imagesc(-nmodes:0,-nmodes:nmodes, trimodes);title('3PWFS Variance of fitted Fourier Amplitude'); colorbar; axis equal; caxis([0, max(max(trimodes))]); xlim([-nmodes-0.5 0.5]); xlabel('n Fourier index'); ylabel('m Fourier index')


%% Variance vs Spatial Frequency Plots

uniqk=unique(k); %find all the unique values of spatial frequency
ssz=size(uniqk); 
for i=1:ssz(2)
    I=find(k==uniqk(i)); %find indicies of unique values
    sz=size(I);
    for j=1:sz(2)
        trivalues(j)=totaltrivar(I(j));
        quadvalues(j)=totalquadvar(I(j));
    end
    trispatialfreq(i)=sqrt(sum(trivalues.^2));
    quadspatialfreq(i)=sqrt(sum(quadvalues.^2));
end

%% Plot

figure; loglog(uniqk, trispatialfreq)%, 'LineWidth', 4);
hold on
loglog(uniqk, quadspatialfreq)%, 'LineWidth',4); %title('Variance of fitted Amplitude vs Spatial Frequency','fontsize',56); ; xlabel('Spatial Frequency (1/m)', 'fontsize',56), ylabel('Varience of the fitted Fourier Mode Amplitudes (rms radians)','fontsize',56);
legend({'3PWFS','4PWFS'},'fontsize',56)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 50)
set(gca,'LineWidth',5)



sum(sum(abs(trispatialfreq)))
sum(sum(abs(quadspatialfreq)))
sum(sum(abs(trispatialfreq)/abs(quadspatialfreq)))

