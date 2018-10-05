function [rpupil,Rwavefront, Pupilpyramid, phase, mcount, ncount]=zernpupil(nn, mm, nphot,rdns,error,npix, Npix, sampling, rmatrix, lambda, tripyramid, MVM, pyramidmask)

%% Zernike Generation

%Scanning through Z -35 to Z 35
counter=1;
success=0;
for n=0:nn
    for m=-mm:mm
ncount(counter)=n;
mcount(counter)=m;
counter=counter+1;
ma = abs(m);
if m==0 & n==0
    continue
elseif mod(n-ma,2)~=0
        continue
    elseif n<ma
        continue
    else  
% Generate the pupil with WF error        
    ef= zernike(0,0,npix).*exp(1i*((2*pi)/lambda)*error*zernike(n,m, npix));        
    pupil = complex(zeros(Npix));
    pupil(Npix/2-npix/2:Npix/2+npix/2-1,Npix/2-npix/2:Npix/2+npix/2-1) =ef;
    success=success+1;
   %ef(:,:,success)=ef;  

%%
%call pyramidsim
[rpupil(:,:,success),Rwavefront(:,:,success), Pupilpyramid(:,:,success), phase(:,:,success)]= pyramidsim(nn,mm,nphot,rdns,error,npix, Npix, sampling, rmatrix, lambda, tripyramid, MVM, pupil, pyramidmask);

    
    end
end
end
    




%%end