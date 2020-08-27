function PLT = phaselagtimestat(data,samplefreq, meanfreq,mw)
% compute the phase lag time for timeseries x and y
% PLT = phaselagtimestat(data,samplefreq, meanfreq, mw)
%
% samplefreq = sampling frequency in Hz
% meanfreq = mean frequnecy of oscillation in Hz
% mw = number of data points for the window
thresh = 1/meanfreq;
[m, n]= size(data);
m = floor(m/2)*2;

K = fix((m-mw/2)/(mw/2)) ;
PLTw = zeros(n,n,K);
for i = 1:K
    
    ibeg = (mw/2)*(i-1) + 1;
    datawind = data(ibeg:ibeg+mw-1,:);
    % step 1, compute the hilbert transform and bring him to the origin
    ht = hilbert(datawind);
    ht = bsxfun(@minus,ht,sum(ht,1)/mw);
    
    % step 2, compute the instantenous phase
    theta = angle(ht);
    
    for j = 2:n
        
        %     x1 = repmat(theta(:, j), 1, j-1);
        %     x2 = theta(:, 1:j-1);
        %     RP = x1 - x2;
        RP = bsxfun(@minus,theta(:, j),theta(:, 1:j-1));
        
        % step 3, phase difference between -pi and pi
        dtheta = mod(RP,2*pi)-pi;
        
        % step 4, is this greater or smaller than zero?
        dtheta = sign(dtheta);
        
        % step 5, find crossings
        cross = dtheta(2:end,:)-dtheta(1:end-1,:);
        for k = 1:j-1
            % step 6, compute time between crossings
            crossings  = find(cross(:,k)) +1;
            time = (crossings(2:end) - crossings(1:end-1))/samplefreq;
            
            % step 7, if time between crossings too small, no synchronization
            time = time((time-thresh)>0);
            % compute phase lag time
            if ~isempty(time)
                PLTw(j,k,i) = mean(1-exp(-time),1);
            end
        end
    end
end
PLT = mean(PLTw,3);

PLT = PLT + PLT';

end