function PLV = phaselockvaluedyn(data,mw)
% compute the phase lock value for timeseries x and y
% PLV = phaselockvaluestat(data, mw)
% mw = number of data points for the window

% cuts points from hilbert transform
cut = 25;
[m, n]= size(data);
mw = floor(mw/2)*2;

% step 1, compute the hilbert transform and bring him to the origin
ht = hilbert(data);   
ht(m-cut+1:m,:) = [];
ht(1:cut,:) = []; 
ht = bsxfun(@minus,ht,mean(ht,1));

% step 2, compute the instantenous phase
theta = angle(ht);

K = fix((m-mw/2-cut*2)/(mw/2));
PLV = zeros(n,n,K);
PLVw = zeros(n,n);
for ii = 1:K
    ibeg = (mw/2)*(ii-1) + 1 ;
    thetaw = theta(ibeg:ibeg+mw-1,:);   
    
    for jj = 2:n             
        RP = bsxfun(@minus,thetaw(:, jj),thetaw(:, 1:jj-1));        
        PLVw(jj,1:jj-1) = abs(sum(exp(1i*RP),1)/mw);       
    end
    PLV(:,:,ii) = PLVw+PLVw';
end

end