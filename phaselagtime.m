function PLT = phaselagtime(data, samplefreq, meanfreq)
% PLT = phaselagtime(data, samplefreq, meanfreq)
% compute the phase lag time for timeseries x and y
% samplefreq is samplefrequency
% meanfreq=f, synchronization/phase difference time should be larger than 1/f 

thresh = 1/meanfreq;
[m, n ]= size(data);
PLT = zeros(n,n);

% step 1, compute the hilbert transform and bring him to the origin
ht = hilbert(data);
ht = bsxfun(@minus,ht,sum(ht,1)/m); 

ht(m-24:m,:)=[];
ht(1:25,:)=[];

% step 2, compute the instantenous phase
theta = angle(ht);

for j = 2:n
    
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
            PLT(j,k) = mean(1-exp(-time),1);
        end
    end
end
PLT = PLT + PLT';

end

