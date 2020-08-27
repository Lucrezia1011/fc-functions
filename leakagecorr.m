function P = leakagecorr(Z)

R = size(Z,2);

D = eye(R,R);
% economy size decomposition

[U,~,V] = svd(Z*D,0);
O = U*V';
d = diag(Z'*O);
D = eye(R).*repmat(d,1,R);

e0 = trace(Z'*Z) - 2*trace(Z'*O*D) + trace(D^2);

[U,~,V] = svd(Z*D,0);
O = U*V';
d = diag(Z'*O);
D = eye(R).*repmat(d,1,R);
e1 = trace(Z'*Z) - 2*trace(Z'*O*D) + trace(D^2);

% semilogy([1,2],[e0,e1],'*'); hold on
while e0-e1 > 1e-6
    

    [U,~,V] = svd(Z*D,0);
    O = U*V';
    d = diag(Z'*O);
    D = eye(R).*repmat(d,1,R);
    e0 = e1;
    e1 = trace(Z'*Z) - 2*trace(Z'*O*D) + trace(D^2);
    % semilogy(j,e1,'*'); hold on
end   
P = O*D;
% figure
% imagesc(corr(Z)); colorbar
% 
% figure
% imagesc(corr(P)); colorbar




