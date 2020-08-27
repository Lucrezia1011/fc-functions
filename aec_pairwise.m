function AEC_wind = aec_pairwise(nodeData,N)
% Envelope correlation with 50% overlapping windows and pairwise
% orthogonalization for leakage correction.
%
% Needs script "leakage_reduction.mexa64"
%
% AEC_wind = aec_pairwise(nodeData,N)
% AEC_wind : 3D matrix, nodes x nodes x windows
% nodeData : beamformed timecourses, time x nodes
% N : length of window in data points, e.g. for a window of 1.5s and data
% acquired at 600Hz => N = 1.5 * 600 = 900.

[m,R] = size(nodeData);
N = round(N/2)*2; % The points in the window need to be even for 50% overlapping windows
cut = 10;  % cuts points from hilbert transform

K = fix((m-N/2-cut*2)/(N/2)); % Number of windows
Neff = N + cut*2;

xm = zeros(Neff*K, R*(R-1)); % Corrected timecourses
for j =1:K
    ibeg = (N/2)*(j-1) + 1;
    data = nodeData(ibeg:ibeg+Neff-1,:);
    for jj = 1:R
        y = data(:,jj);
        ii =  [1:jj-1,jj+1:R];
        for iii =  1:R-1
            x = data(:,ii(iii));
            xm((j-1)*Neff+1:j*Neff, (jj-1)*(R-1)+iii) = leakage_reduction(x,y);
            
        end
    end
end
htm = hilbert(xm);

ht = hilbert(nodeData);
ht([1:cut,end-cut+1:end],:)=[];
ht = bsxfun(@minus, ht, mean(ht,1));
env = abs(ht);

AEC_wind = zeros(R,R,K);
aec = zeros(R,R);

for j =1:K
    
    ibeg = (N/2)*(j-1) + 1;
    imbeg = (j-1)*Neff +1 +cut;
    htx = htm(imbeg:imbeg+N-1,:);
    htx = bsxfun(@minus,htx,sum(htx,1)/N); % normalize
    
    envx = abs(htx);    
    
    A=corr(envx,env(ibeg:ibeg+N-1,:));
    
    for jj = 1:R
        ii =  [1:jj-1,jj+1:R];       
        aec(ii,jj) = A((jj-1)*(R-1)+(1:R-1),jj);
        
    end
    AEC_wind(:,:,j) = (aec + aec')/2;
    
    
end
