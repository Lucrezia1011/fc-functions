function [surr,it] = iaaft(s)
% Surrogata data generation via Iterative amplitude adjusted Fourier
% transform (IAAFT) algorithm, based on Schreiber and Schmitz 1996.
%%

[m,n] = size(s);
if mod(m,2) == 0
    s = s(1:m-1,:);
    m = m-1;
end

% sortrows? to mainting linear corr
[ssort,~] = sort(s);

S = fft(s);
% Amplitude of Fourier transform
Samp = abs(S);
% To have a reference for convergence with iterations

% % Random shuffle of the data
% rshuffle = randperm(size(s,1));
% 
% sshuffle0 = s(rshuffle,1);
% Sshuffle0 = fft(sshuffle0);
% phaseshuffle = angle(Sshuffle0);
% 

% Get parameters for phase randomisazion
len_ser = (m-1)/2;
interv1 = 2:len_ser+1;
interv2 = len_ser+2:m;

% Fourier transform of the original dataset
fft_recblk = fft(s);

ph_rnd = rand([len_ser,1]);
% ph_rnd = rand([len_ser,n]);

% Create the random phases for all the time series
ph_interv1 = repmat(exp( 2*pi*1i*ph_rnd),1,n);
% ph_interv1 = exp( 2*pi*1i*ph_rnd);
ph_interv2 = conj( flipud( ph_interv1));

% Randomize all the time series simultaneously
fft_recblk_surr = fft_recblk;
fft_recblk_surr(interv1,:) = fft_recblk(interv1,:).*ph_interv1;
fft_recblk_surr(interv2,:) = fft_recblk(interv2,:).*ph_interv2;

% Inverse transform
s0= real(ifft(fft_recblk_surr));

% rank order to match initial data
s1 = zeros(size(s));
[~,indrank] = sort(s0);
for nn = 1:n    
    s1(indrank(:,nn),nn) = ssort(:,nn);    
end

S1 = fft(s1);
% epsi2 = norm(S1-S);
% epsi1 = norm(S1-S);

% epsi1 =1;

% Stamdard criteria is that seuqence will remain the same, not feasible for
% multivariate data. 


% epsi1 = epsi2;

% Run forst iteration to check if series converges
phaseshuffle = angle(S1);
%  phaseshuffle = atan(bsxfun(@rdivide, imag(S1),real(S1)));
S2 = Samp.*exp(1i*phaseshuffle);
s2 = ifft(S2);

% rank order to match initial data
s3 = zeros(size(s));

for nn = 1:n
    [~,indrank] = sort(s2(:,nn));
    s3(indrank,nn) = ssort(:,nn);
    
end
S3 = fft(s3);

epsi2 = norm(S3-S);
fprintf('Norm difference: %e\n',epsi2)

S1 = S3;


it = 1;

% while epsi2 - epsi1 < 0 && it < 20
while  it < 20
%     s1 = s3;
%     epsi1 = epsi2;
    
    phaseshuffle = angle(S1); 
%     phaseshuffle = atan(bsxfun(@rdivide, imag(S1),real(S1)));
    
    S2 = Samp.*exp(1i*phaseshuffle);
    s2 = ifft(S2);
    
    % rank order to match initial data
    s3 = zeros(size(s));
    
    [~,indrank] = sort(s2);
    for nn = 1:n       
        s3(indrank(:,nn),nn) = ssort(:,nn);
       
    end
    S3 = fft(s3);
    
%     epsi2 = norm(S3-S);
    fprintf('Norm difference: %e\n',epsi2)
    
    S1 = S3;    
    it =it+1; 

end

surr = s3;
% plot(s1(:,1))



% for ii = 1:100
% 
% Sshuffle1 = Samp.*exp(1i*phaseshuffle);
% sshuffle1 = real(ifft(Sshuffle1));
% 
% % rank order to match initial data
% s2 = zeros(size(s));
% [~,indrank] = sort(sshuffle1);
% s2(indrank,1) = ssort;
% 
% S2 = fft(s2);
% epsi(ii) = norm(S2-S);
% phaseshuffle = angle(S2);
% end
% 
% plot(1:ii,epsi)