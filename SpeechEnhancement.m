clear all; close all; clc;

% Audio acquisition
[file,path] = uigetfile('./Databases/*.wav', 'Select the speech files', 'MultiSelect', 'on'); 
[audio_signal,fs] = audioread([path,file]);
Fs = fs; N = 1024; seqAxis = 0:1:N-1; step = 0;

audio_signal = audio_signal.'/max(abs(audio_signal));
audio_signal = cat(2,zeros(1,N),audio_signal);

% Additive noise + filtering @8KHz
noise = randn(1,length(audio_signal));
nsVar = 0.05;
noise = nsVar * noise/max(abs(noise));

Hd = mylowpass;
audio_noisy = audio_signal + noise;
audio_noisy = filter(Hd,audio_noisy);

% SNR values 
Pa = bandpower(audio_signal);
Ps = bandpower(audio_noisy);
Pn = bandpower(noise);
SNR = 10*log10(Pa/Pn);

% Window choice and preparation
h = hann(N); % Hanning
%h = hamming(N); % Hamming
h = h'; weight = sum(abs(h.^2));

% Wiener Filter parameters
umax = 10; u0 = (1+4*umax/5); z = 25/(umax-1);

%% MTS: Sine Tapers + Noise estimation

m = 0:N-1; L = 5; l = 1:L;
a = sqrt(2/(N+1)).*sin((pi * l' * (m+1))/(N+1));

figure(1); clf
set(gcf,'Name','Sine tapers')
S_k_n = 0.0;

for currentTaper = 1:L
    S_k_n = S_k_n + abs(fft(a(currentTaper,:) .* audio_noisy(1:N))).^2;
    plot(a(currentTaper,:));
    hold on
end

S_n = S_k_n/L; % Noise estimation from speech-absent frame

title('Sine tapers: up to order L = 5')
legend('a = 1','a = 2','a = 3','a = 4','a = 5')
axis tight

%% Standard deviation plot

n = randn(1,N); n = n - mean(n); n_MT = 0.0;

for currentTaper = 1:L
    s = abs(sum(a(currentTaper,:) .* n .* exp(-1i*seqAxis'.*m))).^2;
    n_MT = n_MT + s;
end

n_MT = n_MT/L;

Q = logspace(-4,4,1000);
cnt = 0;

for q = Q
    cnt = cnt+1;
    sigmaZ(cnt) = sqrt(var(log10(n_MT+q)));
end

figure(2); clf
set(gcf,'Name','Relation of q to the standard deviation σz');
loglog(Q,sigmaZ)
title('σ_z vs q for σ_n = 1')
xlabel('q'); ylabel('σ_z')

%% PSD estimation w/ Multitapers

freq = linspace(0,1, N+1);
freq = freq(1:end-1);
step = step + 1;

% Windowing and overlap
audio_noisy_ol = buffer(audio_noisy,N,N/2);
audio_noisy_ol = audio_noisy_ol.';

audio_signal_ol = buffer(audio_signal,N,N/2);
audio_signal_ol = audio_signal_ol.';

nwind = size(audio_noisy_ol,1); 

% MTS estimation
for i = 1:nwind
    S_k = 0.0; 
    for currentTaper = 1:L
    S_k = S_k + abs(fft(a(currentTaper,:) .* audio_noisy_ol(i,:))).^2;
    end
    S_MT(i,:) = S_k/L;
    
    % Verification
    S_x_true(i,:) = abs(sqrt(1/N) * fft(audio_signal_ol(i,:))).^2;
end

wndw = 10;
figure(3); clf
set(gcf, 'name', ['Multitaper log PSD (Window ',num2str(wndw),')']);
title('log PSD comparison')
plot((0:N/2-1)*Fs/N,log(S_x_true(wndw, 1:N/2)),'r');
hold on
plot((0:N/2-1)*Fs/N,log(S_MT(wndw,1:N/2)),'b');
xlabel('Frequency')

%% DWT + Thresholding + Wiener Filter

logSMT = log(S_MT); R = 2;

for window = 1:nwind
    
    % Wavelet coefficient extraction
    [c_log,l_log] = wavedec(logSMT(window,:),L,'db4');
    pos = l_log(1); w5 = c_log(pos+1:pos+l_log(2));
    pos = l_log(2); w4 = c_log(pos+1:pos+l_log(3));
    pos = l_log(3); w3 = c_log(pos+1:pos+l_log(4));
    pos = l_log(4); w2 = c_log(pos+1:pos+l_log(5));
    pos = l_log(5); w1 = c_log(pos+1:pos+l_log(6));
    
    % Scale coefficient extraction 
    [c,l] = wavedec(S_MT(window,:),1,'db4');
    c1 = c(1:l(1));
    [c,l] = wavedec(S_MT(window,:),2,'db4');
    c2 = c(1:l(1));
    [c,l] = wavedec(S_MT(window,:),3,'db4');
    c3 = c(1:l(1));
    [c,l] = wavedec(S_MT(window,:),4,'db4');
    c4 = c(1:l(1));
    [c,l] = wavedec(S_MT(window,:),L,'db4');
    c5 = c(1:l(1));
    
    % Adaptive Thresholding
    w1_new = AdaptThresh(c1,w1,1,n_MT,R,N);
    w2_new = AdaptThresh(c2,w2,2,n_MT,R,N);
    w3_new = AdaptThresh(c3,w3,3,n_MT,R,N);
    w4_new = AdaptThresh(c4,w4,4,n_MT,R,N);
    w5_new = AdaptThresh(c5,w5,5,n_MT,R,N);
    
    % Coefficients update
    pos = l_log(1); c_log(pos+1:pos+l_log(2)) = w5_new;
    pos = l_log(2); c_log(pos+1:pos+l_log(3)) = w4_new;
    pos = l_log(3); c_log(pos+1:pos+l_log(4)) = w3_new;
    pos = l_log(4); c_log(pos+1:pos+l_log(5)) = w2_new;
    pos = l_log(5); c_log(pos+1:pos+l_log(6)) = w1_new;
    
    % IDWT and exponentiation
    logSMT(window,:) = waverec(c_log,l_log,'db4');
    logSMT(window,N/2+1:end) = fliplr(logSMT(window,1:N/2));
    S_wlet(window,:) = exp(logSMT(window,:));
    
    % Wiener filter + Reconstruction
    localSNR = sum(S_wlet(window,:))/sum(S_n);
    
    if localSNR >= 20
        u = 1;
    elseif localSNR <= -5
        u = umax;
    else 
        u = u0 - 10*log10(localSNR)/z;
    end
    
    g_k(window,:) = S_wlet(window,:)./(S_wlet(window,:)+u*S_n);
    
    xdft(window,:) = g_k(window,:).*sqrt(1/weight) .* ...
                    fft(h.*audio_noisy_ol(window,:));
    
    x_reconstructed(window,:) = ifft(sqrt(weight)*xdft(window,:));
end

%% Reconstruction and results

% Overlap and add (50%) - Reconstruction phase
x_rec = x_reconstructed(1,:);
scaling = 1.0;
for window = 2:nwind
    x_rec((window-1)*N/2+1:window*N/2) = scaling .* ...
        (x_rec((window-1)*N/2+1:window*N/2) + x_reconstructed(window,1:N/2));
    x_rec = cat(2,x_rec,x_reconstructed(window,N/2+1:N));
end

% Sound adaptation and play
x_rec = real(x_rec.'/abs(max(x_rec)));
padding = length(x_rec)-length(audio_noisy);
x_rec = x_rec(padding/2+2:end-padding/2+1);
%soundsc(audio_noisy,Fs);
%soundsc(x_rec,Fs);

figure(3);
title('Wavelet thresholding of the logPSD')
plot((0:N/2-1)*Fs/N,logSMT(wndw, 1:N/2),'k');
legend('Clean','Noisy','Thresholded');
axis tight

figure(4); clf
set(gcf, 'name', 'Audio waveform comparison');
subplot(2,1,1)
plot(audio_noisy,'k')
title(['Speech enhancement, N = ',int2str(N),', 50% overlap, Step ',int2str(step)]);
hold on
plot(x_rec,'r--')
legend('Noisy speech','Enhanced speech')
axis tight
xlabel('Samples')
subplot(2,1,2);
plot(audio_signal,'b')
title('Original "clean" data')
legend('Original clean speech')
axis tight
xlabel('Samples')

%% To repeat the denoising procedure:

% Suggested when noise is highly relevant before first denoising
audio_noisy_start = audio_noisy;
audio_noisy = x_rec;
% Execute the code again, starting from the MTS estimation
% (section PSD estimation w/ Multitapers)
