function ph=kolphase(s,method)

%% Lauren Edit:
% Only returning the phase screen computed from the real part

% ph = kolphase(s,'method')
%	Computes a pair of Kolmogorov phase screens using the FT of
%	random complex numbers with appropriate amplitudes.
%	Screens are computed on a grid of size 2s, with s x s sized
%	pieces cut out from the center and returned. This helps overcome
%	the problem of under-representing tilt.
%
% ph = kolphase(s,'lane')
%	Uses the subharmonic method of Lane, Glindemann, & Dainty:
%	"Simulation of a Kolmogorov phase screen," Waves in Random Media 2,
%	209-224 (1992)
%	This explicitly adds in undersampled low-frequency components with
%	spatial scales larger than s to prevent a lack of low-frequency
%	power. More accurate than 'roddier'.
%
% ph = kolphase(s,'roddier')
%	Uses the method of Roddier: "The Effects of Atmospheric Turbulence
%	in Optical Astronomy," Progress in Optics, 19, 281-376 (1981)
%	Random tilts are explicitly added to the phase screens to give
%	a reasonable approximation of the overall Kolmogorov power spectrum
%	at low frequency. Faster than 'lane'.
%
%	size(ph) = [s s 2]

%ph=zeros(s,s,2);
ph = zeros(s,s);

[x y] = meshgrid([-s:s-1],[-s:s-1]);
r = sqrt(x.*x+y.*y);

f = 2*pi*rand(2*s);
pconv = r.^(-11/6);
pconv(s+1,s+1) = 0;
psub = zeros(2*s);
p1 = floor(s/2);
p2 = p1+s-1;

if strcmp(method,'lane')
   for n=1:5
      w = 3^-(2*n);
      for k=-1:1
         ks = k*3^-n;
         sx = Sinc(pi*(x-ks));
         if k==0; t=2; else; t=1; end
         for l=-1:t:1
            ls = l*3^-n;
            sy = Sinc(pi*(y-ls));
            psub = psub+w*sx.*sy.*exp(1i*rand*2*pi);
         end
      end
   end
   xt = [0 0];
   yt = [0 0];
elseif strcmp(method,'roddier')
   scale = 15.2;
   xt = randn([1,2])*scale/s;
   yt = randn([1,2])*scale/s;
else
   error('Unknown method specified');
end

sc = fft2(fftshift((pconv+psub).*exp(1i*f)));
ph(:,:) = real(sc(p1:p2,p1:p2))+xt(1)*x(p1:p2,p1:p2)+yt(1)*y(p1:p2,p1:p2);
%ph(:,:,2) = imag(sc(p1:p2,p1:p2))+xt(2)*x(p1:p2,p1:p2)+yt(2)*y(p1:p2,p1:p2);
ph(:,:) = ph(:,:,1)-mean(mean(ph(:,:,1)));
%ph=ph.*sqrt(1/sum(sum(ph.*ph)));

%ph(:,:,2) = ph(:,:,2)-mean(mean(ph(:,:,2)));
end
