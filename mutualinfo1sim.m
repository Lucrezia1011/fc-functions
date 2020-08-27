clear all
load uncorrected/111514_MEG_4-Restin_rmegpreproc_ROInets_correction-none_13-30Hz_ROI_timecourses

n = 360;
m = 39;

nodeData = nodeData';
l = size(nodeData,1);
l = floor((l-1)/2)*2 +1; % makes m odd
nodeData = nodeData(1:l,:);

cut = 25;
Neff = n + cut*2;
K = fix((l-n/2-cut*2)/(n/2));
ht = hilbert(nodeData);
ht = bsxfun(@minus,ht,mean(ht,1));
ht(end-cut+1:end,:) = [];
ht(1:cut,:) = [];

phiall = angle(ht);

I = zeros(m,m,K);
k = round(0.032*n + 1.5);
 
for j = 1:K
    
    ibeg = (n/2)*(j-1) + 1;        
    phi = phiall(ibeg: ibeg+n-1,:);

    psiterm = zeros(m,m);
    hterm =  zeros(m,m);
    for nodex = 1:m
        x = phi(:,nodex);
        xdiff = pi - abs(abs(bsxfun(@minus,x ,x'))-pi);
        xdiff = xdiff + pi*eye(n);
        for nodey = nodex+1:m
            
            y = phi(:,nodey);
            ydiff = pi - abs(abs(bsxfun(@minus,y ,y'))-pi);
            
            ydiff = ydiff + pi*eye(n);
            zdiffmax = max(cat(3,xdiff,ydiff),[],3);
            
            
            zdiffmax = sort(zdiffmax,1);
            epsi = zdiffmax(k,:);
                                 
            hterm(nodex,nodey) = sum(log(epsi));
            
            nx = sum(bsxfun(@lt,xdiff,epsi),1);
            ny = sum(bsxfun(@lt,ydiff,epsi),1);
            
            psiterm(nodex,nodey) = sum(psi(nx+1)+psi(ny+1));
        end
    end
    
    psiterm = psiterm+psiterm';
    hterm = hterm + hterm';
    
    I(:,:,j) = psi(k) -psiterm/n + psi(n);
%     H(:,:,j) = -psi(k) + psi(n) + 2*hterm/n;
    
%     MI(:,:,j) = ~eye(m).*I(:,:,j)./H(:,:,j);
    
    
end

MI = I/I(1,1,1);
MI = bsxfun(@times,MI,~eye(m));


%%

sim =phaseran2(nodeData,1);

ht = hilbert(sim);
ht = bsxfun(@minus,ht,mean(ht,1));
ht(end-cut+1:end,:) = [];
ht(1:cut,:) = [];

phiall = angle(ht);

Isim = zeros(m,m,K);


for j = 117:K
    
    ibeg = (n/2)*(j-1) + 1;   
    
    
    phi = phiall(ibeg: ibeg+n-1,:);   
    
    psiterm = zeros(m,m);
    hterm =  zeros(m,m);
    for nodex = 1:m
        x = phi(:,nodex);
        xdiff = pi - abs(abs(bsxfun(@minus,x ,x'))-pi);
        xdiff = xdiff + pi*eye(n);
        for nodey = nodex+1:m
            
            y = phi(:,nodey);
            ydiff = pi - abs(abs(bsxfun(@minus,y ,y'))-pi);
            
            ydiff = ydiff + pi*eye(n);
            zdiffmax = max(cat(3,xdiff,ydiff),[],3);           
           
                   
            zdiffmax = sort(zdiffmax,1);
            epsi = zdiffmax(k,:);
          
                      
            hterm(nodex,nodey) = sum(log(epsi));
            
            nx = sum(bsxfun(@lt,xdiff,epsi),1);
            ny = sum(bsxfun(@lt,ydiff,epsi),1);
            
            psiterm(nodex,nodey) = sum(psi(nx+1)+psi(ny+1));
        end
    end
    
    psiterm = psiterm+psiterm';
    hterm = hterm + hterm';
    
    Isim(:,:,j) = psi(k) -psiterm/n + psi(n);

    
    
end

MIsim = Isim/Isim(1,1,1);
MIsim = bsxfun(@times,MIsim,~eye(m));