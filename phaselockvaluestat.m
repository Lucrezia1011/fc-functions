function PLV = phaselockvaluestat(data,mw)
% compute the phase lock value for timeseries x and y
% PLV = phaselockvaluestat(data, mw)
% mw = number of data points for the window

% cuts points from hilbert transform
cut = 25;
mweff = mw + cut*2;
[m, n]= size(data);
mw = floor(mw/2)*2;

K = fix((m-mw/2-cut*2)/(mw/2));
PLVw = zeros(n,n,K);
for i = 1:K
    
    ibeg = (mw/2)*(i-1) + 1 ;
    datawind = data(ibeg:ibeg+mweff-1,:);
    % step 1, compute the hilbert transform and bring him to the origin
%     datawind = econleakagecorr(datawind);
    ht = hilbert(datawind);   
    ht = bsxfun(@minus,ht,sum(ht,1)/mweff);
    ht(mweff-cut+1:mweff,:) = [];
    ht(1:cut,:) = []; 
    
    % step 2, compute the instantenous phase
    theta = angle(ht);
    
    for j = 2:n        
      
        RP = bsxfun(@minus,theta(:, j),theta(:, 1:j-1));        
        PLVw(j,1:j-1,i) = abs(sum(exp(1i*RP),1)/mw);
        
    end
end
PLV = mean(PLVw,3);

PLV = PLV + PLV';

end