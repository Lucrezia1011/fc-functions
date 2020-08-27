function PLI = phaselagindex(data)
% compute the phase lag index for timeseries x and y

[m, n ]= size(data);

PLI = zeros(n,n);

% step 1, compute the hilbert transform and bring him to the origin
ht = hilbert(data);
ht = bsxfun(@minus,ht,sum(ht,1)/m); 

ht(m-24:m,:)=[];
ht(1:25,:)=[];

m = m-50;
% step 2, compute the instantenous phase
theta = angle(ht);

for jj = 2:n    

    RP = bsxfun(@minus,theta(:, jj),theta(:, 1:jj-1));  
    RPe = exp(1i*RP);
    RPs = sign(imag(RPe));
    PLI(jj,1:jj-1) = abs(sum(RPs,1)/m);           
        
end
   
PLI = PLI + PLI';

end