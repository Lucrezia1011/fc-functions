function PLI = phaselagindexmlc(data,mw)
% compute the phase lock value for timeseries x and y
% PLV = phaselockvaluestat(data, mw)
% mw = number of data points for the window

% cuts points from hilbert transform
cut = 25;
mweff = mw + cut*2;
[m, n]= size(data);
mw = floor(mw/2)*2;

K = fix((m-mw/2-cut*2)/(mw/2));
PLI = zeros(n,n,K);
PLIw = zeros(n,n);
for ii = 1:K
    ibeg = (mw/2)*(ii-1) + 1 ;
    datawind = data(ibeg:ibeg+mweff-1,:);    
    datawind = econleakagecorr(datawind);  
    ht = hilbert(datawind);   
    ht = bsxfun(@minus,ht,sum(ht,1)/mweff);
    ht(mweff-cut+1:mweff,:) = [];
    ht(1:cut,:) = [];
    theta = angle(ht);
    
    for jj = 2:n        
      
        RP = bsxfun(@minus,theta(:, jj),theta(:, 1:jj-1));  
        RPe = exp(1i*RP);
        RPs = sign(imag(RPe));
        PLIw(jj,1:jj-1) = abs(sum(RPs,1)/mw);               
    end
    PLI(:,:,ii) = PLIw+PLIw';
end


end