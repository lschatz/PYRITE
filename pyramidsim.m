function [rpupil, Rwavefront, phase]=pyramidsim(add_noise, nphot,lambda, error, fbs, rmatrix, kpupil, rdns, npix, Npix, pyramidmask, sampling, tripyramid, MVM, broke)

%% Set up

% Centers of the pyramid pupils
% Sampling= # of pixels across the pyramid pupil 
if broke==false
if sampling==16
    cen1=[81,128];
    cen2=[128,81];
    cen3=[81,33];
    cen4=[33,81];
end

if sampling==32
%     cen1=[38,38];
%     cen2=[92,38];
%     cen3=[38,92];
%     cen4=[92,92];

      cen1=[161,255];
      cen2=[255,161];
      cen3=[161,66];
      cen4=[66,161];
end

if sampling==64
%     cen1=[75,75];
%     cen2=[183,75];
%     cen3=[75,183];
%     cen4=[183,183];

      cen1=[321,510];
      cen2=[510,321];
      cen3=[321,132];
      cen4=[132,321];
end

if sampling==128
%     cen1=[150,150];
%     cen2=[366,150];
%     cen3=[150,366];
%     cen4=[366,366];
      cen1=[641,1019];
      cen2=[1019,641];
      cen3=[641,264];
      cen4=[264,641];
end
end

if broke==true
if tripyramid==false  
    if sampling==16
        cen1=[47,47];
        cen2=[47, 114];
        cen3=[114,47];
        cen4=[114,114];
    end
end

if tripyramid==true
    if sampling==16
        cen1=[40,57];
        %cen1= [57,40];
        cen2=[81,128];
        %cen2= [128,81]:
        cen3=[121,57];
      
        %cen3= [57,121];
        cen4= [57,57];
    end
end
end
szz=size(kpupil);

for i=1:szz(3)
pupil=kpupil(:,:,i);

% Zero pad the pupil
pupil=padarray(pupil,[Npix/2-npix/2, Npix/2-npix/2],0,'both');
% Add in photons
mask=zernike(0,0,npix);
pmask=padarray(mask,[Npix/2-npix/2, Npix/2-npix/2],0,'both');
msize=sum(sum(pmask));
p=sqrt(nphot/msize(1));
pupil=pupil.*p;
    

%% Fourier Transform to focal plane
FTpupil=fftshift(fft2(fftshift(pupil)))/Npix;
%figure, imagesc(abs(FTpupil).^2), title('PSF')

%% Apply Pyramid Phase Mask
pyramid=FTpupil.*pyramidmask; %Apply the OPD mask to simulate the pyramid tip

% figure;imagesc(angle(pyramid));axis equal; title('Pyramid focal plane
% splitting')
%figure; imagesc(angle(pyramid));axis equal

%% Back to Pupil Plane and simulate detection
Pupilpyramids=abs((ifftshift(fft2(ifftshift(pyramid))))/Npix).^2;

%figure; imagesc(Pupilpyramid); axis equal


%% Bin the pupil pixels to get correct detector sampling
% [a,b,c]=size(Pupilpyramids);
% p=npix/sampling;
% 
% 
% for i= 1:c
% Pupilpy=Pupilpyramids(:,:,i);    
% Pupilpy=sum(reshape(Pupilpy, p, []),1);
% Pupilpy=reshape(Pupilpy, a/p, []).';
% Pupilpy=sum(reshape(Pupilpy, p, []),1);
% Pupilpy=reshape(Pupilpy, b/p, []).';
% %figure;imagesc(Pupilpyramid); axis equal; title('After binning');
% Pupilpyramid(:,:,i)=Pupilpy;
% end

[a,b]=size(Pupilpyramids);
p=npix/sampling;

Pupilpyramids=sum(reshape(Pupilpyramids, p, []),1);
Pupilpyramids=reshape(Pupilpyramids, a/p, []).';
Pupilpyramids=sum(reshape(Pupilpyramids, p, []),1);
Pupilpyramids=reshape(Pupilpyramids, b/p, []).';

Pupilpyramid(:,:,i)=Pupilpyramids;


end



%% Intensity Centroid Calculation MVM
if MVM==true
PupilOne=Pupilpyramid(cen1(1)-sampling/2-1:cen1(1)+sampling/2+1, cen1(2)-sampling/2-1:cen1(2)+sampling/2+1,:);
PupilTwo=Pupilpyramid(cen2(1)-sampling/2-1:cen2(1)+sampling/2+1, cen2(2)-sampling/2-1:cen2(2)+sampling/2+1,:);
PupilThree=Pupilpyramid(cen3(1)-sampling/2-1:cen3(1)+sampling/2+1, cen3(2)-sampling/2-1:cen3(2)+sampling/2+1,:);
PupilFour=Pupilpyramid(cen4(1)-sampling/2-1:cen4(1)+sampling/2+1, cen4(2)-sampling/2-1:cen4(2)+sampling/2+1,:);

%% Add noise
if add_noise==true
PupilOne=noise(PupilOne,-1,rdns);
PupilTwo=noise(PupilTwo,-1,rdns);
PupilThree=noise(PupilThree,-1,rdns);
PupilFour=noise(PupilFour,-1,rdns);
end
%% Mask region so only pupils remain
sz=size(PupilOne);
mask=circle(sz(1), sampling+1);
PupilOne=PupilOne.*mask;
PupilTwo= PupilTwo.*mask;
PupilThree=PupilThree.*mask;
PupilFour=PupilFour.*mask;
%%
for i=1:szz(3)
   
P1=PupilOne(:,:,i);
P2=PupilTwo(:,:,i);
P3=PupilThree(:,:,i);
P4=PupilFour(:,:,i);
%%  Pull out the pixels only in the pupils (no 0s from mask included)
    
mask=reshape(mask',[size(mask,1)*size(mask,2) 1]);
P1=reshape(P1',[size(P1,1)*size(P1,2) 1]);
P1=P1(mask~=0);
P2=reshape(P2',[size(P2,1)*size(P2,2) 1]);
P2=P2(mask~=0);
P3=reshape(P3',[size(P3,1)*size(P3,2) 1]);
P3=P3(mask~=0);
P4=reshape(P4',[size(P4,1)*size(P4,2) 1]);
P4=P4(mask~=0);


P1dummy(:,i)=P1;
P2dummy(:,i)=P2;
P3dummy(:,i)=P3;
P4dummy(:,i)=P4;

end

PupilOne=P1dummy;
PupilTwo=P2dummy;
PupilThree=P3dummy;
PupilFour=P4dummy; 

%%

% figure
% subplot(2,2,1)
% imagesc(PupilOne);
% axis equal
% subplot(2,2,2)
% imagesc(PupilTwo);
% axis equal
% subplot(2,2,3)
% imagesc(PupilThree);
% axis equal
% subplot(2,2,4);
% imagesc(PupilFour);
% axis equal

%% Calculate Centroids 

sz=size(PupilOne);

%%

if MVM==true

% 4PWFS quad-cell equation
if tripyramid==false
    Sx= (PupilOne+PupilTwo-PupilThree-PupilFour)./(PupilOne+PupilTwo+PupilThree+PupilFour);
    Sy= (PupilOne-PupilTwo-PupilThree+PupilFour)./(PupilOne+PupilTwo+PupilThree+PupilFour);
end

% 3PWFS tri-cell equation
if tripyramid==true
    Sx= (PupilTwo*(sqrt(3)/2)-PupilThree*(sqrt(3)/2))./(PupilOne+PupilTwo+PupilThree);
    Sy= (PupilOne-PupilTwo*(0.5)-PupilThree*(0.5))./(PupilOne+PupilTwo+PupilThree); 
end   
%%
Centroids=[];
for k=1:sz(2)
    
    dummy=[Sx(:,k)' Sy(:,k)'];
    Slopes_mat(:,k)=dummy';
end
end     
            
end



%% Full Frame Method
if MVM==false
    Slopes_mat=Pupilpyramid;
end


%% G-Matrix Generation
ppupil=Pupilpyramid;
Slopes_mat(isnan(Slopes_mat))=0;   


Rwavefront= rmatrix*Slopes_mat;
Rwavefront(abs(Rwavefront)<1e-6)=0;

% % Normalize
% sz=size(kval);
% for i=1:sz(2)
%  Rwavefront(:,i)=Rwavefront(:,i)./kval(:,i);
% end


%% Reconstruct


% Generate the pupil with WF error      
sz=size(Rwavefront);
phase=zeros(npix,npix, sz(2));
for i=1:sz(2)
for j=1:sz(1)
  phase(:,:,i)= phase(:,:,i)+fbs(:,:,j).*Rwavefront(j,i);
end
end
    %pupil = complex(zeros(Npix))
%    success=success+1;
    
%%  
phase=phase.*error.*(2*pi/lambda);
%%
rpupil=zernike(0,0,npix).*exp(1i*phase);

end