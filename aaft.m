function surr = aaft(s)
% Amplitude Adjusted Fourier Transfer, multivariate version


[m,n] = size(s);
if mod(m,2) == 0
    s = s(1:m-1,:);
    m = m-1;
end

% Rank order data according to a set of Gaussian random numbers
sgauss = zeros(m,n);
for ii = 1:n
    [~,indrank] = sort(s(:,ii));
    randgauss = randn(m,1);
    sortgauss = sort(randgauss);
    sgauss(indrank,ii) = sortgauss;
end

% Get parameters for phase randomisazion
len_ser = (m-1)/2;
interv1 = 2:len_ser+1;
interv2 = len_ser+2:m;

% Fourier transform of the gaussian dataset
fft_recblk = fft(sgauss);
ph_rnd = rand([len_ser,1]);

% Create the random phases for all the time series
ph_interv1 = repmat(exp( 2*pi*1i*ph_rnd),1,n);
ph_interv2 = conj( flipud( ph_interv1));

% Randomize all the time series simultaneously
fft_recblk_surr = fft_recblk;
fft_recblk_surr(interv1,:) = fft_recblk(interv1,:).*ph_interv1;
fft_recblk_surr(interv2,:) = fft_recblk(interv2,:).*ph_interv2;

% Inverse transform
surrblk= real(ifft(fft_recblk_surr));

% Invert gaussian transformaation by rank ordering phase
% randomised data according to original data VE

surr = zeros(m,n);
for ii = 1:n
    [~,indrank] = sort(surrblk(:,ii));
    invgauss = sort(s(:,ii));
    surr(indrank,ii) = invgauss;
end