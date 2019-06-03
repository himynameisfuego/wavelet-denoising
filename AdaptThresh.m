function w_new = AdaptThresh(c,w,order,n_MT,R,N)

% w_new = AdaptThresh(c,w,order,n_MT,R,N)
cnt = 0;
qj = zeros(size(c));
qnoise = mean(c);
qj(1) = ((sum(c(1:1+R)))/(R+1) - qnoise)/(2^order); 
qj(2) = ((sum(c(1:1+R+1)))/(R+2) - qnoise)/(2^order);
    
for index = (1+R):(length(c)-R)
    qj(index) = ((sum(c(index-R:index+R)))/(2*R+1) - qnoise)/(2^order);
end

qj(end-1) = ((sum(c(end-R-1:end)))/(R+2) - qnoise)/(2^order); 
qj(end) = ((sum(c(end-R:end)))/(R+1) - qnoise)/(2^order);

for q = qj
    cnt = cnt+1;
    sigmaZ(cnt) = sqrt(var(log10(n_MT+q)));
end

threshold = sqrt(2*log(N))*sigmaZ;
for index = 1:length(w)
    w_new(index) = thresholding(w,index,threshold(1:length(w)));
end

end


function a = thresholding(c,index,threshold)

if c(index) > threshold(index) 
       c(index) = c(index)-threshold(index);
elseif c(index) < -threshold(index)
       c(index) = c(index)+threshold(index);
else c(index) = 0;
end
a = c(index);
end