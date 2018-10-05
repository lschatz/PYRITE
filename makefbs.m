 function [fbs,xp,yp,cp] = makefbs(s,n)

% fbs = makefbs(s,n)
%     Makes a Fourier basis set defined on a square array of
%     linear dimension s out to n cycles/aperture. The number
%     of modes returned in fbs is (2n+1)^2-1.  The constant
%     mode is ignored.  Modes are organized by pairs in increasing
%     spatial frequency with the cosine term followed by the sine
%     term.

nmode = (2*n+1)^2-1;
fbs = zeros(s,s,nmode+2);
cp = floor(s/2)+1;

[x,y] = meshgrid([-n:n],[-n:n]);
[theta,r] = cart2pol(x,y);
[rs,idx] = sort(r(:),'ascend');
thetas = theta(idx);
idx = idx((rs>0)&(thetas<=0)); %chops off DC r=0 component, and half the plane to get rid of double count
xp = x(idx)+cp;
yp = y(idx)+cp;

% Make the first two modes be tip-tilt since they're not
% well modeled by the Fourier basis set.


[fbs(:,:,1),fbs(:,:,2)] = meshgrid([-s/2+1:s/2]/(s/2),[-s/2+1:s/2]/(s/2));
%c=3;

c = 3;


for k=1:nmode/2
    ft = zeros(s);
    ft(yp(k),xp(k)) = 1;
    sw = fftshift(fft2(fftshift(ft)));
    fbs(:,:,c) = real(sw); %cosines
    fbs(:,:,c+1) = imag(sw); %Sines
    c = c+2;
end


