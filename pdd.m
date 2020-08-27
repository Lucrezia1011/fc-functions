%% PDD of data, time points x ROIs
function PDDmat = pdd(data,lowpass,highpass,f)

 minosc =  f/((lowpass+highpass)/2); % minimum number of oscillations to consider FC

PDDmat = zeros(size(data,2)); % nr = number of ROIs
cut = 10; % number of points to cut at data edges after hilbert transform
for ii= 1:78
    datax = data(:,ii);
    for jj = [1:ii-1,ii+1:78]
        datay = data(:,jj);
        datay = leakage_reduction(datay,datax);
        ht = hilbert([datax,datay]);
        phi = unwrap(angle(ht));
        
        RP = phi(:,1)-phi(:,2); % units [radians]
        dRP = diff(RP);   % units [radians/sample]
        
        pdd = exp(-abs(dRP)*minosc);  % units [ radians/sample  * sample]
        
        %% For discontinuities in the data
        %         disc_p = np:np:np*nt-1;
        %                 pdd0 = pdd;
        %                 crossi =diff(pdd> 0.05); % Consider a crossing a pdd greater than 1deg
        %                 crossi(disc_p,:) = 1; % Adds crossings at discontinuities between trials
        %                 crossi(disc_p-1,:) = 1;
        %                 crossi = [1 ;crossi; 1];
        %
        %                 crossf = find(crossi) ;
        %                 timep = diff(crossf); % timepoints between each crossing
        %                 ff = find(timep<minosc); % if crossing is less than 1 oscillation
        %
        %                 cmtimep = cumsum(timep);
        %                 if ff(1) == 1
        %                     pdd(1:cmtimep(1),1) = 0;
        %
        %                     for jjj = 2:length(ff)
        %                         n = ff(jjj);
        %                         pdd(cmtimep(n-1)+1: cmtimep(n), 1) = 0;
        %                     end
        %                 else
        %
        %                     for jjj = 1:length(ff)
        %                         n = ff(jjj);
        %                         pdd(cmtimep(n-1)+1: cmtimep(n), 1) = 0;
        %                     end
        %                 end
        %%
        
        pdd = smooth(pdd,5,'moving');
        %                 pdd0 = smooth(pdd0,5,'moving');
        %                 PDD0mat(ii,jj) = sum(pdd0(cut+1:end-cut,:),1)/(nt*np -1 - cut*2);
        PDDmat(ii,jj) = sum(pdd(cut+1:end-cut,:),1)/(size(data,1) -1 - cut*2);
    end
    
end