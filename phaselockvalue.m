function PLV = phaselockvalue(data)
% compute the phaselocking value for timeseries

[m, n ]= size(data);

PLV = zeros(n,n);

% step 1, compute the hilbert transform and bring him to the origin
ht = hilbert(data);
ht = bsxfun(@minus,ht,sum(ht,1)/m); 

ht(m-24:m,:)=[];
ht(1:25,:)=[];

m = m-50;
% step 2, compute the instantenous phase
theta = angle(ht);

for j = 2:n    
    RP = bsxfun(@minus,theta(:, j),theta(:, 1:j-1)); 
    PLV(j,1:j-1) = abs(sum(exp(1i*RP),1)/m);    
end
   
PLV = PLV + PLV';

end
