function rp = randp(mean)

% randp(mean)
% Returns an array of Poisson random numbers of the same size as mean

rp = zeros(size(mean));

p2g = 20;

% Flag number big enough that we can just use a Gaussian distribution
pmask = logical(ones(size(mean)));
pmask(mean>p2g) = false;
pg = randgauss(mean,mean);
rp(~pmask) = round(pg(~pmask));
rp(rp<0) = 0;

mean(~pmask) = 0;

p = rand(size(mean));
e = exp(-mean);
pn = zeros(size(mean));

done = false;
while ~done
    pn(pmask) = pn(pmask)+mean(pmask).^rp(pmask)./factorial(round(rp(pmask))).*e(pmask);
    pmask(pn>p) = false;
    rp(pmask) = rp(pmask)+1;
    done = ~logical(sum(pmask(:)));
end
