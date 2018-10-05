function [pupil,fbs]= fbspupil(npix, nmodes, Npix, lambda, error, rnphot)

%% Fourier Mode Generator
[fbs]=makefbs(npix,nmodes);
mask=zernike(0,0,npix);
fbs=fbs.*mask;


%% Normalize fourier modes

fsz=size(fbs);
for i=1:fsz(3)
fbs(:,:,i)=fbs(:,:,i).*sqrt(1/sum(sum(fbs(:,:,i).*fbs(:,:,i))));
end

%% Make into complex Pupil
fphase=zernike(0,0,npix).*exp(1i*((2*pi)/lambda)*error*fbs);
pupil=padarray(fphase,[Npix/2-npix/2, Npix/2-npix/2],0,'both');

%% Add in photons

pmask=padarray(mask,[Npix/2-npix/2, Npix/2-npix/2],0,'both');
msize=sum(sum(pmask));
p=sqrt(rnphot/msize(1));

pupil=pupil.*p;

%% saving
save('fbs.mat','fbs')
%save('fbpupil.mat','pupil')

end