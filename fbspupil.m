function [fbpupil,fbs,xp,yp,cp]= fbspupil(npix, nmodes, lambda, error, rnphot)

%% Fourier Mode Generator
[fbs,xp,yp,cp]=makefbs(npix,nmodes);
mask=zernike(0,0,npix);
fbs=fbs.*mask;


%% Normalize fourier modes

fsz=size(fbs);
for i=1:fsz(3)
fbs(:,:,i)=fbs(:,:,i).*sqrt(1/sum(sum(fbs(:,:,i).*fbs(:,:,i))));% normalize
end

%% Make orthonormal with Stabilized Gram Schmidt

fbs=stabilizedgs(fbs);


%% Make into complex Pupil
fbpupil=zernike(0,0,npix).*exp(1i*((2*pi)/lambda)*error*fbs);








%pupil=padarray(fphase,[Npix/2-npix/2, Npix/2-npix/2],0,'both');

% %% Add in photons
% 
% %pmask=padarray(mask,[Npix/2-npix/2, Npix/2-npix/2],0,'both');
% msize=sum(sum(mask));
% p=sqrt(rnphot/msize(1));
% 
% pupil=pupil.*p;

%% saving
save('fbs.mat','fbs')
save('fbpupil.mat','fbpupil', '-v7.3')
save('xp.mat','xp')
save('yp.mat','yp')
save('cp.mat','cp')

end