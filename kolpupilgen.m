function [kscreen,kpupil]=kolpupilgen(error,npix, Npix, lambda,nphot)

%% Set Up Masks
mask=zernike(0,0,npix);
msize=sum(sum(mask)); 
p=sqrt(nphot/msize);
%msize=sum(sum(mask));
%p=sqrt(nphot/msize); % number of photons in the pupil plane, so that you get the correct number in the focal plane.
%% Generate Masks

    screen=kolphase(npix,'roddier');
    kscreen=screen*((2*pi)/lambda).*error;

    kscreen=subtilt(kscreen); %subtract out tip/tilt so we don't break the pyramid
    kscreen=kscreen.*mask;
    
    piston=mean(mean(kscreen)); 
    kscreen=kscreen-piston; %subtract out piston.
    kscreen=kscreen.*mask;
%% Put phase masks into a complex pupil

kpupil=zernike(0,0,npix).*exp(1i*kscreen);


%    ef= zernike(0,0,npix).*exp(1i*kscreen);        
%     pupil = complex(zeros(Npix));
%     pupil(Npix/2-npix/2:Npix/2+npix/2-1,Npix/2-npix/2:Npix/2+npix/2-1)=ef;
%     kpupil=pupil.*p;
%   


end

