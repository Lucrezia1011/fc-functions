function [staticcorr, dynamiccorr] = envelopecorr_averaged(data,freq,Tav,Twind)
% [staticcorr, dynamiccorr] = envelopecorr_averaged(data,freq,Tav,Twind)
% freq = sampling frequency of the signal (in Hz)
% Tav = window for averaging (in seconds)
% Twind = window for dynamic connectivity measure (in seconds) 
% data needs to be leakage corrected

Nw = freq*Tav;

[m,n] = size(data);
K1 = fix((m-Nw/2)/(Nw/2));
env = zeros(K1,n);

% Real data
for j = 1:K1
    ibeg = (Nw/2)*(j-1) + 1;
    ht = hilbert(data(ibeg: ibeg+Nw-1,:));
    ht(1:5,:) = [];
    ht(end-4:end,:) = [];
    ht = bsxfun(@minus,ht,sum(ht,1)/(Nw-10));
    env(j,:) = mean(abs(ht),1);
end

staticcorr = corr(env).*~eye(n,n);

N = (round(Twind/Tav*2)); % finds equivalent time 
N = N-mod(N,2);
K = fix((size(env,1)-N/2)/(N/2)) ;
dynamiccorr = zeros(n,n,K);
for j = 1:K
    ibeg = (N/2)*(j-1) + 1;
    dynamiccorr(:,:,j)  = corr(env(ibeg: ibeg+N-1,:)).*~eye(n);
end

end