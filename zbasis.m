function [zpupil,zbs]=zbasis(nn, mm, npix, Npix, nmodes,lambda,error)

%% Zernike Generation

G_mat=[];
success=1;
%Scanning through Z -35 to Z 35
 for n=0:nn
    for m=-mm:mm
ma = abs(m);
if m==0 & n==0
    continue
elseif mod(n-ma,2)~=0
        continue
    elseif n<ma
        continue
    else  
% Generate the pupil with WF error   

    z=zernike(n,m,npix).*sqrt(1/sum(sum(zernike(n,m,npix).*zernike(n,m,npix)))); %normalize
    zbs(:,:,success)=z;
    zpupil(:,:,success)= zernike(0,0,npix).*exp(1i*((2*pi)/lambda)*error.*z); 
    success=success+1;
end
    end
 end
 

end
