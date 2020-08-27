function MI = mutualinformation_Peng(data)
% kth neighbour, default k =1
% Mutual information based on method 1, only 1 trial
ht = hilbert(data);
ht = bsxfun(@minus,ht,mean(ht,1));
ht(1:5,:) = [];
ht(end-4:end,:) = [];
phi = angle(ht);
m = size(phi,2);

I = zeros(m,m);
H = zeros(m,m);
for nodex = 1:m
    x = phi(:,nodex);    
    for nodey = nodex+1:m
            y = phi(:,nodey);  
            [p12, p1, p2] = estpab(x,y);
            I(nodex,nodey) = estmutualinfo(p12,p1,p2);            
            H(nodex,nodey) = estjointentropy(p12);
    end
end

I = I+I';
H = H+H';
MI = I./H;

end
