function AEC = envelopecorr(data)
ht = hilbert(data);
[m,n]= size(data);

ht(m-24:m,:)=[];
ht(1:25,:)=[];
ht = bsxfun(@minus,ht,mean(ht,1)); 

envelope = abs(ht);
AEC = corr(envelope).*~eye(n,n);
end