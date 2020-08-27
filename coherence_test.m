% Calculates dynamic connectivity metrics and surrogates.
% version 3: adjusted the hilbert transform to cut off 50 points at the
% ends while manteining the correct window and overlap.

clear all
close all
clc
addpath /net/carados/data_local/Lucrezia/ROI-timecourses/uncorrected/
addpath /home/ppyll1/matlab/svdandpca
addpath /home/ppyll1/matlab/Repeatibility/

acq_freq = 300; % frequency of signal (Hz)
R = 39; % number of ROIs

% minimun common multiple of the bandwidths is 3060
% N needs to be a multiple of 2 >> T multiple of 6120


subnum = dlmread('names_uncorrected2.txt');
s_tot = length(subnum);


freq =1;
for effp = 30
    
    if freq == 1
        highf = 13;
        lowf = 8;
        folder = 'alpha';
        
    elseif freq == 2
        highf = 30;
        lowf = 13;
        folder = 'beta';
        
    elseif freq == 3
        highf = 48;
        lowf = 30;
        folder = 'gamma';
        
    elseif freq == 4
        highf = 4;
        lowf = 1;
        folder = 'delta';
        
    elseif freq == 5
        highf = 8;
        lowf = 4;
        folder = 'theta';
    end
    
    
    f = [num2str(lowf),'-',num2str(highf)];
    B = highf-lowf;
    T = effp/(2*B);
    % N needs to be a multiple of 2
    N = round(T*acq_freq/2)*2;
    
    % cuts points from hilbert transform
    cut = 10;
    Neff = N + cut*2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters for coherence calculation
    Kw = fix(7*N/acq_freq);
    if Kw < 2
        Kw = 2;
    end
    p = .75; % overlap of windows for coherence
    if p == .5
        M = 2*fix(fix((1/(1-p))*N / (Kw+1))/4);
    elseif p == .75
        M = 4*fix(fix(4*N / (Kw+3)) /4);
    end
    
    nfft =  2^(ceil(log2(M)));
    if nfft < 256
        nfft = 256;
    end
    fx = linspace(0,acq_freq/2,nfft/2+1);
    fband = fx>= lowf & fx <= highf;
    % makes hamming envelope for each window
    w = hamming(M);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subfolder =[folder,'_allwindows_N',num2str(N)];
    
   
    
    warning( 'off','MATLAB:MKDIR:DirectoryExists')
    
    disp([folder,' f. band, window ',num2str(T),'s'])
    aecc = [];
    plvc = [];
    plic = [];
    cohc = [];
    icohc = [];
    iplvc = [];

    for s = 56:60
        
        %% Metrics calculations
        for t = 1:3
            
            name = [num2str(subnum(s)),'_MEG_',num2str(t+2),...
                '-Restin_rmegpreproc_ROInets_correction-none_',f,'Hz_ROI_timecourses.mat'];
            if exist(name, 'file') == 2
                
                
                load(name);
                nodeData = nodeData';
                
                m = size(nodeData,1);
                m = floor((m-1)/2)*2 +1; % makes m odd
                nodeData = nodeData(1:m,:);
                
                fprintf('loaded %6.0f subject %2.0f, session %2.0f\n',subnum(s),s,t)
                
                % Calculation optimised for small data points in each window
                K = fix((m-N/2-cut*2)/(N/2));
                
                AEC_wind = zeros(R,R,K);
                PLI_wind = zeros(R,R,K);
                PLV_wind = zeros(R,R,K);
                COH_wind = zeros(R,R,K);
                iCOH_wind = zeros(R,R,K);
                iPLV_wind = zeros(R,R,K);
                
                aec = zeros(R,R);
                plv = zeros(R,R);
                iplv = zeros(R,R);
                pli = zeros(R,R);
                coh = zeros(R,R);
                icoh = zeros(R,R);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                               
                % makes hamming envelope for each window
                
                X = zeros(nfft,R,Kw);
                Xm = zeros(nfft,R*(R-1),Kw);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % no correction
                ht = hilbert(nodeData);
                ht([1:cut,end-cut+1:end],:)=[];
                ht = bsxfun(@minus, ht, mean(ht,1));
                phi = angle(ht);
                env = abs(ht);
                clear ht
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Orthoghonalization
                xm = zeros(Neff*K, R*(R-1)); % Corrected timecourses
                for j =1:K
                    ibeg = (N/2)*(j-1) + 1;
                    data = nodeData(ibeg:ibeg+Neff-1,:);
                    for jj = 1:R
                        y = data(:,jj);
                        ii =  [1:jj-1,jj+1:R];
                        for iii =  1:R-1
                            x = data(:,ii(iii));
                            xm((j-1)*Neff+1:j*Neff, (jj-1)*(R-1)+iii) = leakage_reduction(x,y); 
                        end
                    end
                end
                htm = hilbert(xm);
                
                
                for j =1:K
                    
                    ibeg = (N/2)*(j-1) + 1;
                    imbeg = (j-1)*Neff +1 +cut;
                    htx = htm(imbeg:imbeg+N-1,:);
                    htx = bsxfun(@minus,htx,sum(htx,1)/N); % normalize
                    phix = angle(htx);
                    envx = abs(htx);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Non corrected imaginary coherence
                  
                    for kk = 1:Kw
                        iibeg = (M*(1-p))*(kk-1);
                        % selects 50% overlapping windows
                        % All fourier transforms for time segments
                        X(:,:,kk) = fft(bsxfun( @times,nodeData(iibeg+cut+ibeg:iibeg+cut+ibeg+M-1,:),w), nfft);
                        Xm(:,:,kk) = fft(bsxfun( @times,xm(iibeg+imbeg:iibeg+imbeg+M-1,:),w), nfft);
                    end
                    Xc = X(fband,:,:);
                    Xmc = Xm(fband,:,:);
                    % Takes sum of the overalpping hamming windows
                    % Could also use mean, but the denominator cancels with Pxy
                    P = sum(abs(Xc).^2,3);
                    Pm = sum(abs(Xmc).^2,3);
                    
                    for jj = 2:R
                        RP = bsxfun(@minus,phi(ibeg:ibeg+N-1,jj),phi(ibeg:ibeg+N-1, 1:jj-1));
                        RPe = imag(exp(1i*RP));
                        iplv(jj,1:jj-1)=abs(sum(RPe,1)/N);                                            
                        pli(jj,1:jj-1) = abs(sum(sign(RPe),1)/N);
                        Pxy = sum(bsxfun(@times,conj(Xc(:,jj,:)),Xc(:, 1:jj-1,:)),3);
                        Sxy = imag(Pxy).^2 ./ (bsxfun( @times,P(:,jj),P(:,1:jj-1)));
                        icoh(jj,1:jj-1)  = mean(Sxy,1);
                    end
                    
                    iPLV_wind(:,:,j) = iplv + iplv';
                    PLI_wind(:,:,j) = pli + pli';
                    iCOH_wind(:,:,j) = icoh + icoh';
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    A=corr(envx,env(ibeg:ibeg+N-1,:));
                    
                    for jj = 1:R                       
                                              
                        ii =  [1:jj-1,jj+1:R];                        
                        aec(ii,jj) = A((jj-1)*(R-1)+(1:R-1),jj);
                        RP = bsxfun(@minus,phix(:,(jj-1)*(R-1)+(1:R-1)),phi(ibeg:ibeg+N-1,jj));
                        plv(ii,jj)=abs(sum(exp(1i*RP),1)/N);
                        
                        Pxy = sum(bsxfun(@times,conj(Xmc(:, (jj-1)*(R-1)+(1:R-1),:)),Xc(:,jj,:)),3);
                        Sxy = abs(Pxy).^2 ./ (bsxfun( @times,Pm(:,(jj-1)*(R-1)+(1:R-1)),P(:,jj)));
                        coh(ii,jj)  = mean(Sxy,1);
                        
                    end
                    
                    AEC_wind(:,:,j) = (aec + aec')/2;
                    PLV_wind(:,:,j) = (plv + plv')/2;
                    COH_wind(:,:,j) = (coh + coh')/2;
                end
                
  
%                 save(['/net/carados/data_local/Lucrezia/ROI-timecourses/bandwidth_window/sim_allwindows_pairwise_test/'...
%                     ,subfolder,'/simZstats_sub',num2str(s),'_',...
%                     num2str(t+2),'_',f,'Hz_N',num2str(N),'_pairwise_test.mat'],...
%                     'AEC_wind','PLV_wind','iPLV_wind','PLI_wind','COH_wind','iCOH_wind');%                 
                
                aecc = cat(3,aecc,AEC_wind);
                plvc = cat(3,plvc,PLV_wind);
                iplvc = cat(3,iplvc,iPLV_wind);
                plic = cat(3,plic,PLI_wind);
                cohc = cat(3,cohc,COH_wind);
                icohc = cat(3,icohc,iCOH_wind);
              
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
            else
                disp(['Skipped sub',num2str(s),', trial ',num2str(t+2),', ',f,'Hz'])
                
                
            end
            

        end
        
    end
end


%% Find network of FC edges correlated in time
addpath (genpath('/home/ppyll1/matlab/brain_map'))

fc = plvc; 
Ktot = size(fc,3);
c = repmat(triu(ones(R),1),[1,1,Ktot]);
ll =find(triu(ones(R),1));
% v = reshape(fc,[R^2, Ktot]);
v = reshape(fc(c==1),[741, Ktot]);
Cv = corr(v').*~eye(741);


figure; imagesc(Cv)


excursion = zeros(R);
for ii = 1:R-1
    for jj = ii+1:R
        % Difference of FC from the median
        diff = squeeze(fc(ii,jj,:)) - median(fc(ii,jj,:),3);
        % Positive or negative difference
        cross = diff>0;
        % Find crossing points
        dexcurs = cross(2:end) - cross(1:end-1);
        crossings  = find(dexcurs) +1;
        % adds first crossing
        crossings = cat(1,1,crossings);
        % Finds length of excursion
        In = crossings(2:end) - crossings(1:end-1);
        % Find height of excursion
        Hn = zeros(length(crossings)-1,1);
        for kk = 1:length(crossings)-1
            Hn(kk) = max(abs(diff(crossings(kk):crossings(kk+1)-1)));
        end
        
        excursion(ii,jj) = sum(In.^.9.*Hn);
    end
end
excursion = excursion + excursion';

figure; go_3Dbrain(excursion,.9)


corrpar = .2;
% Find edges pairs
C = Cv > corrpar;
figure; imagesc(C)
colormap gray
[nComponents,sizes,members] = networkComponents(C);


colrange = [min(excursion(~eye(R))), max(excursion(~eye(R)))];
%%
% Order based on correlation 
% ind = zeros(1,nnz(sizes>2));
% for ii = 1:nnz(sizes>2)
%     lc = nchoosek(members{ii},2);
%     llc = lc(:,1) + (lc(:,2)-1)*741;
%     cllc =Cv(llc);
%     ind(ii) = mean(cllc(cllc>corrpar));
% %     ind(ii) = mean(cllc);
% end
% 
% [~,indord] = sort(ind,'descend');
% 
% for ii = indord(1:3)
%     Z = NaN(R);
%     Z(ll(members{ii}))= excursion(ll(members{ii}));      
%     
%     
%     figure;   
%     go_3Dbrain(Z,0,colrange);
% end

% Order based on variance
exll = zeros(1,nnz(sizes>2));
for ii = 1:nnz(sizes>1)
    exll(ii) = mean(excursion(ll(members{ii})));    
end

[~,indord] = sort(exll,'descend');


for ii = indord(1:3)
    Z = NaN(R);
    Z(ll(members{ii}))= excursion(ll(members{ii}));      
    
    
    figure;   
    go_3Dbrain(Z,0,colrange);
end

%%
for ii = 1:5
    Z = NaN(R);
    Z(ll(members{ii}))= excursion(ll(members{ii}));      
    
    
    figure;   
    go_3Dbrain(Z,0,colrange);
end

