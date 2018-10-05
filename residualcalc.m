function [trivariance, quadvariance, quadresidual, triresidual]=residualcalc(triphase, quadphase, kscreen,npix)

%% Calculate Residual Wavefront

triresidual=kscreen-triphase;
quadresidual=kscreen-quadphase;


%% Calculate Variances of Residual Pupil
%calculates the variances of across the pupil, (inside where the mask has
%a value of 1). 

mask=circle(npix,npix);

for i=1:size(triresidual,3)
    triresid=triresidual(:,:,i);
    tripix(:,i)=triresid(mask>0);
    
    quadresid=quadresidual(:,:,i);
    quadpix(:,i)=quadresid(mask>0);
end

trivariance=var(tripix);
quadvariance=var(quadpix);




