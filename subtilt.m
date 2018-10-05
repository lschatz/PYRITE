function notilt = subtilt(phases,varargin);

% notilt = subtilt(phases,varargin);
%
% Removes global tip-tilt and mean value from each frame of a cube of phase screens
% notilt = subtilt(phases)		Removes tilt over non-zero pixels of phases
% notilt = subtilt(phases,mask)		Removes tilt over non-zero pixels of mask

nframe = size(phases,3);
notilt = zeros(size(phases));

if nargin>1
   pmask = varargin{1};
else
   pmask = sum(abs(phases),3);
end
pmask(abs(pmask)>0) = 1;

[x,y] = meshgrid([1:size(pmask,2)],[1:size(pmask,1)]);
for k=1:nframe
    p = double(phases(:,:,k));
    %fo = fit([x(pmask>0),y(pmask>0)],p(pmask>0),'poly11');
    fo = planefit(x(pmask>0),y(pmask>0),p(pmask>0));
    %ft = fo.p00+x*fo.p10+y*fo.p01;
    ft = x*fo(1)+y*fo(2)+fo(3);
    notilt(:,:,k) = (p-ft);
    %notilt(:,:,k) = (p-ft).*pmask;
end
