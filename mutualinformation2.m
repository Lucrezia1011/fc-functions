function MI = mutualinformation2(phi,k)
% Mutual information of the phases,  Lucrezia Liuzzi
% 
% Following rectangular algorithm (2) for estimating mutual information 
% based on k-nearest neighbour distances: 
% Kraskov A. (2004) 'Estimating mutual information'
% DOI: 10.1103/PhysRevE.69.066138
%
% MI = mutualinformation1(phi,k)
% 
% phi = signal angles,  i.e.   angle(hilbert (data))
% dimensions: time x voxels
% 
% k can be omitted: 
% kth neighbour, default k from Wilmer 2012, PLoS ONE 7(9): e44633.
% doi:10.1371/journal.pone.0044633


[n,m] = size(phi);
psiterm = zeros(m,m);

if nargin == 1
    k = round(0.032*n + 1.5);
end


for nodex = 1:m-1
    x = phi(:,nodex);
    % Adjust angle difference to consider circualrity of the phase
    xdiff = pi - abs(abs(bsxfun(@minus,x ,x'))-pi); 
    
    for nodey = nodex+1:m        
            
            y = phi(:,nodey);  
            % Adjust angle difference to consider circualrity of the phase
            ydiff = pi - abs(abs(bsxfun(@minus,y ,y'))-pi);        
                   
            zdiffmax = max(cat(3,xdiff,ydiff),[],3);
            
            [~,ind] = sort(zdiffmax,1);                        
            ind = ind + repmat(0:n:n*(n-1),n,1);           
            
            epsix =  max(xdiff(ind(1:k+1,:)),[],1);
            epsiy =  max(ydiff(ind(1:k+1,:)),[],1);
            
            % Finds points in X which have a distance from x(i) < epsi(i)
            nx = sum(bsxfun(@lt,xdiff,epsix),1) ;
            ny = sum(bsxfun(@lt,ydiff,epsiy),1) ;            
            
           
                %%
%                 clf
%                 i = 3;
%                 scatter(x,y)
%                 hold on
%                 scatter(x(i),y(i),'r')
%                 plot([x(i)-epsix(i),x(i)+epsix(i),x(i)+epsix(i),x(i)-epsix(i),x(i)-epsix(i)],...
%                     [y(i)-epsiy(i),y(i)-epsiy(i),y(i)+epsiy(i),y(i)+epsiy(i),y(i)-epsiy(i)],'k--')
%                 axis equal
                %%
               
                      
            psiterm(nodex,nodey) = mean(psi(nx)+psi(ny));
        
    end
end

nodey = nodex;
y = phi(:,nodey);
ydiff = pi - abs(abs(bsxfun(@minus,y ,y'))-pi);
zdiffmax = max(cat(3,xdiff,ydiff),[],3);
[~,ind] = sort(zdiffmax,1);
ind = ind + repmat(0:n:n*(n-1),n,1);
epsix =  max(xdiff(ind(1:k+1,:)),[],1);
epsiy =  max(ydiff(ind(1:k+1,:)),[],1);
nx = sum(bsxfun(@lt,xdiff,epsix),1) ;
ny = sum(bsxfun(@lt,ydiff,epsiy),1) ;
psiterm_auto = mean(psi(nx)+psi(ny));
Imax =psi(k) -1/k -psiterm_auto + psi(n);

psiterm = psiterm+psiterm';
I = psi(k) -1/k -psiterm + psi(n);

% Normalize I diving by its maximum value
MI = I/Imax;
MI = bsxfun(@times,MI,~eye(m));

end
