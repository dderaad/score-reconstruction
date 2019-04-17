%% (a)
close all; clear variables; clc
load handel

v = y'/2; v(end) = [];
t = (1:length(v))/Fs;

figure(1)
subplot(3,1,1)
plot(t,v);
hold on;
xlabel('Time [sec]');
ylabel('Amplitude');
title('Signal of Interest, v(n)');

% p8 = audioplayer(v,Fs);
% playblocking(p8);

a = .05; % Filter width in seconds
samples = 500; % Number of samples
samples = 100;
tauvec = linspace(t(1), t(end), samples); % Sample points (filter centers)

n = length(v); % fourier modes
L = t(end) - t(1); % length of recording
k=(2*pi/L)*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);

for filtertype = 1:3
spec = zeros([length(ks)/2, length(tauvec)]);

for j = 1:length(tauvec)
    switch filtertype
        case 1
            subplot(3,1,1)
            filter = Shannonfilter(t,a,tauvec(j));%
            plot(t,v);
            hold on
            plot(t,filter/(2*max(filter))); axis([0 t(end) -1 1])
            hold off
            subplot(3,1,2)
            vf = filter.*v;
            plot(t,filter.*v); axis([0 t(end) -1 1])
            subplot(3,1,3)
            vft = fftshift(fft(vf));
            spec(:,j) = abs(vft(end/2+1:end)');
            plot(ks/(2*pi), abs(vft));
            drawnow
        case 2
            filter = Gaussianfilter(t,a,tauvec(j));
            vf = filter.*v;
            vft = fftshift(fft(vf));
            spec(:,j) = abs(vft(end/2+1:end)');
        case 3
            filter = MexiHatfilter(t,a,tauvec(j));%
            vf = filter.*v;
            vft = fftshift(fft(vf));
            spec(:,j) = abs(vft(end/2+1:end)');
    end
end

figure(1+filtertype)
colormap(bone)
pcolor(tauvec, ks(end/2+1:end)/(2*pi), log(1+spec)), shading interp
switch filtertype
    case 1
        title(sprintf('Logscale Spectrogram of Handel''s Messiah \n using Shannon Filter of Filter Width %.3f with %.0f Samples', a, samples))
    case 2
        title(sprintf('Logscale Spectrogram of Handel''s Messiah \n using Gaussian Filter of Filter Width %.3f with %.0f Samples', a, samples))
    case 3
        title(sprintf('Logscale Spectrogram of Handel''s Messiah \n using Mexican Hat Filter of Filter Width %.3f with %.0f Samples', a, samples))
end
ylabel('frequency (Hz)')
xlabel('t (sec)')
end

%% (b)
%close all; clear variables; clc

figure(5)
subplot(2,1,1)
[y1, Fs] = audioread('music1.wav');
tm1 = (1:length(y1))/Fs;
plot(tm1,y1);
ylabel('Amplitude');
title('Mary had a little lamb (piano)');
axis([tm1(1) tm1(end) -1 1])
%p8 = audioplayer(y1,Fs); playblocking(p8);
figure(6)
subplot(2,1,1)
[y2, Fs] = audioread('music2.wav');
tm2 = (1:length(y2))/Fs;
plot(tm2,y2);
ylabel('Amplitude');
title('Mary had a little lamb (recorder)');
axis([tm2(1) tm2(end) -1 1])
%p8 = audioplayer(y2,Fs); playblocking(p8);

a = .0625; % Filter width in seconds
samples = 120; % Number of samples
tauvec1 = linspace(tm1(1), tm1(end), samples); % Sample points (filter centers)
tauvec2 = linspace(tm2(1), tm2(end), samples); % Sample points (filter centers)

n1 = length(y1); % fourier modes
n2 = length(y2);
L1 = tm1(end) - tm1(1); % length of recording
L2 = tm2(end) - tm2(1);
k1=(2*pi/L1)*[0:(n1/2-1) -n1/2:-1]; ks1=fftshift(k1);
k2=(2*pi/L2)*[0:(n2/2-1) -n2/2:-1]; ks2=fftshift(k2);

spec1 = zeros([length(ks1)/2, length(tauvec1)]);

for j = 1:length(tauvec1)
    % filter = Shannonfilter(tm1,a,tauvec1(j));%
    filter = Gaussianfilter(tm1,a,tauvec1(j));%
    % filter = MexiHatfilter(t,a,tauvec(j));%
    y1f = filter.*y1';
    y1ft = fftshift(fft(y1f));
    spec1(:,j) = abs(y1ft(end/2+1:end)');
end

spec2 = zeros([length(ks2)/2, length(tauvec2)]);

for j = 1:length(tauvec2)
    filter = Gaussianfilter(tm2,a,tauvec2(j));%
    
    y2f = filter.*y2';
    y2ft = fftshift(fft(y2f));
    spec2(:,j) = abs(y2ft(end/2+1:end)'); 
end

figure(5)
subplot(2,1,2)
colormap(hot)
pcolor(tauvec1, ks1(end/2+1:end)/(2*pi), spec1), shading interp
axis([tm1(1) tm1(end) .9*261.63 1.1*329.63])
hold on
line([tm1(1) tm1(end) tm1(end) tm1(1) tm1(1) tm1(end)], [261.63 261.63 293.66 293.66 329.63 329.63])
xlabel('Time (sec)'); 
ylabel('Frequency (Hz)');
yticks([261.63 293.66 329.63]);
yticklabels({'C: 261.63','D: 293.66','E: 329.63'});

figure(6)
subplot(2,1,2)
colormap(hot)
pcolor(tauvec2, ks2(end/2+1:end)/(2*pi), spec2), shading interp
axis([tm2(1) tm2(end) .9*783.99 1.1*987.77])
hold on
line([tm2(1) tm2(end) tm2(end) tm2(1) tm2(1) tm2(end)], [783.99 783.99 880 880 987.77 987.77])
xlabel('Time (sec)');
ylabel('Frequency (Hz)');
yticks([783.99 880 987.77]);
yticklabels({'G: 783.99','A: 880','B: 987.77'});
