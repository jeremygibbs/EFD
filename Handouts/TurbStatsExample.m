%%% homework #1 mfile
close all;clear all;clc;warning off

%%% problem #3
U=load('v1obsdata.txt');

%%% plot the pdf with 25 bins using function hist %%%
[Uh,bins] = hist(U,25);
Updf = Uh/(sum(Uh)*(bins(2)-bins(1)));

figure;plot(bins,Updf);
xlabel('velocity [m/s]');ylabel('pdf(U)');

%%% stats (using matlab functions var and mean, could also just use
%%% function sum)

display(['mean =',num2str(mean(U),4)])
display(['var  =',num2str(var(U),4)])
display(['skew =',num2str(mean((U-mean(U)).^3)/var(U)^(3/2),4)])
display(['kurt =',num2str(mean((U-mean(U)).^4)/var(U)^(2),4)])
% 
%%% plot the autocorrelation
Uprime = U-mean(U); %define primes

for k=0:length(U)/2-1 %lags to 1/2 the timeseries length (could do fewer)
    autocorr(k+1) = mean(Uprime(1:length(U)/2).*Uprime(k+1:k+length(U)/2));
end
autocorr = autocorr/var(U);

%%% plot the autocorrelation vs time lag (divide by 20 for aquisition
%%% speed of 20Hz)
figure;plot([0:4000]/20,autocorr(1:4001));
hold on; plot([0,4000/20],[0,0],'--k');
xlabel('time [s]');ylabel('autocorrelation');

clear U;
%%% problem #4

data=load('v1simdata.txt');
U = reshape(data(:,3),256,256);

Uprime = U-mean(mean(U)); %define primes (not required)

for i=1:256
    uk = fft(Uprime(:,i),[],1);
    psp(:,i) = uk.*conj(uk);
end

psp = squeeze(mean(psp,2));
wn = 2*pi*[1:1:128]/(256.0*3.125);
figure;loglog(wn,psp(1:128),wn,2.5*wn.^(-5/3.0));set(gca,'Xlim',[wn(1) wn(128)]);
xlabel('wavenumber k');ylabel('Power E(k)');

figure;plot(wn,psp(1:128));set(gca,'Xlim',[wn(1) wn(128)]);
xlabel('wavenumber k');ylabel('Power E(k)');
