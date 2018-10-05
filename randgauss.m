function rg = randgauss(mean,variance)

% rg = randgauss(mean,variance)
%    Returns an array of Gaussian random numbers of the same size as mean and variance
%    mean and variance can be scalars, vectors, or 2D arrays.
%    If they are both not scalars, they must be the same size.

if isequal(size(mean),size(variance))
    rg = zeros(size(mean));
    v = variance;
    m = mean;
else
    if isequal(size(mean),[1 1])
        rg = zeros(size(variance));
        m = ones(size(variance))*mean;
        v = variance;
    elseif isequal(size(variance),[1 1])
        rg = zeros(size(mean));
        v = ones(size(mean))*variance;
        m = mean;
    else
        error('randgauss: mean and variance must be equal sized arrays or scalars\n');
    end
end

randnum = 1000;

for k = 1:randnum; rg = rg+rand(size(rg)); end
rg = (rg-randnum/2)*3.4459.*sqrt(v/randnum)+m;

