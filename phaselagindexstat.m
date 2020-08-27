function PLI = phaselagindexstat(data,mw)
% compute the phase lag index for timeseries x and y
% mw = number of data points for the window
% PLI = phaselagindexstat(data, mw)
% cuts points from hilbert transform

cut = 25;
mweff = mw + cut*2;

[m, n]= size(data);
m = floor(m/2)*2;
K = fix((m-mw/2-cut*2)/(mw/2));

PLIw = zeros(n,n,K);
for i = 1:K
    
    ibeg = (mw/2)*(i-1) + 1;
    datawind = data(ibeg:ibeg+mweff-1,:);    
    % step 1, compute the hilbert transform and bring him to the origin
    datawind = econleakagecorr(datawind);
    ht = hilbert(datawind);
    ht = bsxfun(@minus,ht,sum(ht,1)/mweff);   
    ht(mweff-cut+1:mweff,:) = [];
    ht(1:cut,:) = []; 
    % step 2, compute the instantenous phase
    theta = angle(ht);
    
    for j = 2:n
        
      
        RP = bsxfun(@minus,theta(:, j),theta(:, 1:j-1));        
        % step 3, phase difference between -pi and pi
        dtheta = mod(RP,2*pi)-pi;
        
        % step 4, is this greater or smaller than zero?
        dtheta = sign(dtheta);
        PLIw(j,1:j-1,i) = abs(sum(dtheta,1)/mw);
        
    end
end
PLI = mean(PLIw,3);

PLI = PLI + PLI';

end