function  [coh,icoh] = coherence_dyn(data,f,highf,lowf,filtyn)
% [coh,icoh] = coherence_LL(data,f,highf,lowf)
% Measure of coherence based on power spectrum by Welch's method
% with 75% overlapping Hamming windows, optimized for small windows 
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

p = 0.9; % overlap of windows for coherence
p1= 0.1;

% if p == 0.5
%     M = 2*fix(fix((1/(1-p))*N / (Kw+1))/4);
% elseif p==.75
%     M = 4*fix(fix(4*N / (Kw+3)) /4);
% else
    M = round(round((n1/4)*(p1))/ (p1));
    K = fix( (n1- (M*p) )  / (M*(p1)) ) ;
% end


% K = ceil(n1/f*7); 
% if K <2
%     K = 2;
% end
% 
% M = 4*fix(fix(4*n1 / (K+3)) /4);

% if isempty(M) && ~isempty(K)
%     M = (1/(1-p))*n1 / (K+1); 
%     M = round(M*(1-p))/(1-p);   
%         
% elseif  ~isempty(M) && isempty(K)
%     M = round(M*(1-p))/(1-p);
%     M = round(round(M/(1-p))*(1-p));
%     K = fix(n1/((1-p)*M) -1);
% elseif isempty(M) && isempty(K)
%     K = 8; % sets 8 windows to default
%     M = (1/(1-p))*n1 / (K+1);
%     M = round(M*(1-p))/(1-p); 
% end

 
%  set N to the next power of 2 greater than M (from matlab default settings of pwelch)
N = 2^(ceil(log2(M)));   
if N < 256
    N = 256;
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

coh = zeros(n2,n2);
icoh = zeros(n2,n2);

X = zeros(N,n2,K);
% makes hamming envelope for each window
w = hamming(M);

for ii = 1:K
    ibeg = (M*(p1))*(ii-1) + 1;   
    % selects 50% overlapping windows
    % All fourier transforms for time segments
    X(:,:,ii) = fft(bsxfun( @times,data(ibeg:ibeg+M-1,:),w), N);
end
% Second half of fourier transform is just a mirror of the data
if filtyn == 0
    X = X(freq>= lowf & freq <= highf,:,:);    
else
    B = highf - lowf; % Bandwidth
    X = X(freq>= (lowf-B/2) & freq <= (highf+B/2),:,:);    
end

% Takes sum of the overalpping hamming windows
% Could also use mean, but the denominator cancels with Pxy
P = sum(abs(X).^2,3); 

for jj = 2:n2
    Pxy = sum(bsxfun(@times,conj(X(:,jj,:)),X(:, 1:jj-1,:)),3);
    Sxyden = sqrt(bsxfun( @times,P(:,jj),P(:,1:jj-1)));    
    coh(jj,1:jj-1)  = mean(abs(Pxy) ./ Sxyden,1);    
    icoh(jj,1:jj-1)  = mean(abs(imag(Pxy)) ./ Sxyden, 1);
end

coh = coh + coh';
icoh = icoh + icoh';
end
