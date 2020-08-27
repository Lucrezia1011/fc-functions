function coh = coherence(data, nfft, Fs, filt, n_overlap )
% computing coherence
%
% All input parameters are equivalent to csd or the other spectrum 
% related function in MATLAB.
%
% This program use csd.m in Signal Processing Toolbox.
%
% Input
%	f1 and f2	input data
%	nfft		number of data for FFT
%	Fs		sampling frequency
%	filt		filter vector, hanning(nfft/2)
%	n_overlap	number of overlap for smoothing
%
% Output
%	coh		coherence
%========================================================================

coh = zeros(size(data,2),size(data,2));
for i = 1: size(data,2)
    for j = 1: size(data,2)
        
        f1 = data(:,i);
        f2 = data(:,j);
        if i<j            
            
            [b,bint,f2] = regress(f2,f1);
            
            [Pxx,freq] = pwelch( f1,filt, n_overlap, nfft, Fs);
            [Pyy,freq] = pwelch( f2,filt, n_overlap, nfft, Fs);            
            [Pxy,freq] = cpsd( f2, f1, filt, n_overlap, nfft, Fs );
%             Pxx(1:2)=[];
%             Pyy(1:2)=[];
%             Pxy(1:2)=[];
            
            coh(i,j)  =mean(abs(Pxy.^2)./(Pxx.*Pyy));
            
        elseif j>i
                       
            
            [b,bint,f1] = regress(f1,f2);
            
            [Pxx,freq] = pwelch( f1,filt, n_overlap, nfft, Fs);
            [Pyy,freq] = pwelch( f2,filt, n_overlap, nfft, Fs);            
            [Pxy,freq] = cpsd( f1, f2, filt, n_overlap, nfft, Fs );
%             Pxx(1:2)=[];
%             Pyy(1:2)=[];
%             Pxy(1:2)=[];
            
            coh(i,j)  = mean(abs(Pxy.^2)./(Pxx.*Pyy));
        end
    end
end


coh = (coh + coh')/2;

end