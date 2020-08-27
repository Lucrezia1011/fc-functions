function  icoh = imag_coherence(data,M,f,highf,lowf)
% icoh = imag_coherence(data, M,f,highf,lowf)
%
% Measure of imaginary coherence based on power spectrum by Welch's method
% with 50% overlapping Hamming windows 
%
% data (time x voxels)
% M = Hamming window length in number of points (e.g. M = 10*600, for a 10s window at 600Hz)
% f = frequency of the signal (e.g. 600Hz)
% highf, lowf = frequency limits (e.g. highf = 30, lowf = 13, for beta band)
% N = FFT length
%
% Power Spectrum Estimation by Welch's method based on:
%---------------------------------------------------------------
% copyright 1994, by C.S. Burrus, J.H. McClellan, A.V. Oppenheim,
% T.W. Parks, R.W. Schafer, & H.W. Schussler.  For use with the book
% "Computer-Based Exercises for Signal Processing Using MATLAB"
% (Prentice-Hall, 1994).

% http://www.rpi.edu/dept/ecse/rta/burrus/Func-v4/welch.m
% accessed 06/10/2015
%---------------------------------------------------------------

[n1,n2]= size(data);

if ~isempty(M)
    M = M - mod(M,2); % needs to be an even number for 50% overlap windows
    K = fix((n1-M/2)/(M/2)) ;
else
    K = 8; 
    M = fix(2*n1 / 9); % length of fft section (and Hamming window) with Welch's method, set to use 8 windows
end

%     
if M <=256
    N = 256;
else
%  set N to the next power of 2 greater than M (from matlab default settings of pwelch)
    N = 2^(ceil(log2(M)));    
end

% Default frequency 600Hz
if isempty(f)
    f= 600;
end

if isempty(lowf)
    lowf = 0;
end

if isempty(highf)
    highf = f/2; % Nyquist frequency 
end

freq = linspace(0,f/2,N/2+1);

icoh = zeros(n2,n2);

X = zeros(N,n2,K);
% makes hamming envelope for each window
w = hamming(M);

for ii = 1:K
    ibeg = (M/2)*(ii-1) + 1;   
    % selects 50% overlapping windows
    % All fourier transforms for time segments
    X(:,:,ii) = fft(bsxfun( @times,data(ibeg:ibeg+M-1,:),w), N);
end
X = X(freq>= lowf & freq <= highf,:,:);
% Takes sum of the overalpping hamming windows
% Could also use mean, but the denominator cancels with Pxy
P = sum(abs(X).^2,3); 
% Second half of fourier transform is just a mirror of the data
% P = P(1:N/2+1,:);

for jj = 2:n2
    Pxy = sum(bsxfun( @times,conj(X(:,jj,:)),X(:, 1:jj-1,:)),3);
%     Pxy = Pxy(1:N/2+1,:);     
    icoh(jj,1:jj-1)  = mean(imag(Pxy./sqrt(bsxfun( @times,P(:,jj),P(:,1:jj-1)))),1);
end
% Imaginary coherence can be positive or negative depending on the order of
% conj(X)*Y and is antisymmetric along the diagonal -> Take absolute value
icoh = abs(icoh + icoh');

end


%% If working with matlab 2015
% See help file for pwelch or cpsd for more information on input parameters

% function  icoh = icoherence_2015(data,M,NOVERLAP,N)%     
% %     data (time x voxels)
% 
% %     M = WINDOW.  If WINDOW is omitted or specified as
% %     empty, a Hamming window is used to obtain eight sections of X and Y
% % 
% %     NOVERLAP samples of overlap from
% %     section to section. NOVERLAP must be an integer smaller than the
% %     length of WINDOW if WINDOW is a vector, or smaller than WINDOW if
% %     WINDOW is an integer. If NOVERLAP is omitted or specified as empty, it
% %     is set to obtain a 50% overlap.
% % 
% %     N = FFT length
% 
% 
% n = size(data,2);
% icoh = zeros(n,n);
% 
% % P = pwelch(data, M,NOVERLAP,N);
% 
% for j = 2:n
%     data1 = repmat(data(:,j),1,j-1);
%     data2 = data(:,1:j-1);
%     Pxy = cpsd(data1, data2, M,NOVERLAP,N);  
%     icoh(1:jj-1,jj)  = mean(imag(Pxy./sqrt(bsxfun( @times,P(:,jj),P(:,1:jj-1)))),1);     
% end
% 
% icoh = icoh + icoh';
% 
% end