function MI = mutualinformation1(phi,k)
% Mutual information of the phases,  Lucrezia Liuzzi
% 
% Following square algorithm (1) for estimating mutual information based on
% k-nearest neighbour distances: 
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

if nargin == 1
    k = round(0.032*n + 1.5);
end

psiterm = zeros(m,m);
% hterm =  zeros(m,m);
for nodex = 1:m-1
    x = phi(:,nodex);
    % Adjust angle difference to consider circularity of the phase
    xdiff = pi - abs(abs(bsxfun(@minus,x ,x'))-pi);
    
    for nodey = nodex+1:m

            y = phi(:,nodey); 
            % Adjust angle difference to consider circularity of the phase
            ydiff = pi - abs(abs(bsxfun(@minus,y ,y'))-pi);                
            zdiffmax = max(cat(3,xdiff,ydiff),[],3);            
            
            zdiffmax = sort(zdiffmax,1);          
            epsi = zdiffmax(k+1,:);
                       
%%            
%                 clf
% %                 i = 1490;
%                 scatter(x,y)
%                 hold on
%                 scatter(x(i),y(i),'r')
%                 plot([x(i)-epsi(i),x(i)+epsi(i),x(i)+epsi(i),x(i)-epsi(i),x(i)-epsi(i)],...
%                     [y(i)-epsi(i),y(i)-epsi(i),y(i)+epsi(i),y(i)+epsi(i),y(i)-epsi(i)],'k--')
%                 axis equal
            %%
%             hterm(nodex,nodey) = sum(log(epsi));
            
            % Finds points in X which have a distance from x(i) < epsi(i)
            nx = sum(bsxfun(@lt,xdiff,epsi),1);
            ny = sum(bsxfun(@lt,ydiff,epsi),1);           
                       
            psiterm(nodex,nodey) = sum(psi(nx+1)+psi(ny+1));
    end
end

nodey = nodex;
y = phi(:,nodey);
ydiff = pi - abs(abs(bsxfun(@minus,y ,y'))-pi);         
zdiffmax = max(cat(3,xdiff,ydiff),[],3);
zdiffmax = sort(zdiffmax,1);
epsi = zdiffmax(k+1,:);
% Finds points in X which have a distance from x(i) < epsi(i)
nx = sum(bsxfun(@lt,xdiff,epsi),1);
ny = sum(bsxfun(@lt,ydiff,epsi),1);
psiterm_auto = sum(psi(nx+1)+psi(ny+1));
Imax = psi(k) -psiterm_auto/n + psi(n);

psiterm = psiterm+psiterm';
% hterm = hterm + hterm';

I = psi(k) -psiterm/n + psi(n);
% H = -psi(k) + psi(n) + 2*hterm/n;

% Normalize I diving by its maximum value
MI = I/Imax;
MI = bsxfun(@times,MI,~eye(m));


end
