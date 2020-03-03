%   Simulating a PSK transmission
%   Nb bits are randomly generated and modulated to PSK
%   The random bits and the signal waveform are plotted
%   In this case each bit is one symbol. So the symbols will be omitted
%   By Marcelino Silva, March 2nd, 2020

%Optional
clear all;
close all;

%Parameters
Tb = 1e-3;  %Bit duration (s) and Symbol duration
t = [0:Tb/99:Tb]; %Time interval to plot one waveform period

nc = 2;     %Cycles number
fc = nc/Tb;  %Frequency

Eb = 1e-6;  %Bit energy

Nb = 100;   %Number of bits to be transmitted

J = 5;   %Number of symbols to plot the signal waveform (symbols being transmitted)

%Base function
phit = sqrt(2/Tb).*cos(2*pi*fc.*t);

%Signals: s1(t) and s2(t)
s1t = sqrt(Eb).*phit;    %Bit 1
s2t = -1*sqrt(Eb).*phit; %Bit 0

%Random bits to be transmitted
bits = randn(Nb,1);
bits(find(bits<=0)) = -1; % -1 means bit 0
bits(find(bits>0)) = +1; % +1 means bit 1

%Plot the simulation

%Index of bits being trasmitted
beginI = 1;
endI = J;

%plotting bits
subplot(2,1,1);
hold on;
axis([1 Nb -1.1 1.1]);
bar((1:beginI),bits(1:beginI),'b');
bar((beginI:endI),bits(beginI:endI),'r');
bar((endI:Nb),bits(endI:Nb),'k');
title('Bits being transmitted');
legend('Transmitted Bits','Bits being transmitted','Bits to be transmitted');

%vector to signal waveform
lt = length(t);
s = zeros(1,lt*J);

while(endI<=Nb),
    
    %plotting bits
    b = bits(beginI:endI);   
    subplot(2,1,1);
    bar((1:beginI),bits(1:beginI),'b'); %Transmitted Bits
    bar((beginI:endI),bits(beginI:endI),'r'); %Bits being transmitted
    bar((endI:Nb),bits(endI:Nb),'k'); %Bits to be transmitted
    
    
    %signal wavefor for the bits being transitted
    for i=1:J,
        if(b(i)==1),
            s(lt*(i-1)+1:lt*i) = s1t;
        else
            s(lt*(i-1)+1:lt*i) = s2t;
        end;
    end;
    
    %plotting signal
    subplot(2,1,2);
    ts = [(Tb/100)+(beginI-1)*Tb:Tb/100:endI*Tb];
    plot(ts,s,'r');
    axis([(Tb/100)+(beginI-1)*Tb endI*Tb -1*sqrt(2*Eb/Tb) sqrt(2*Eb/Tb)]);
    title('PSK Signal being transmitted');
    xlabel('Time (s)');
    ylabel('Signal');
    
    %update index
    beginI = beginI + 1;
    endI = endI + 1;
    
    pause(0.5);
    
end;

%Plotting the costellation and the decision regions

v = [-sqrt(Eb) 0 sqrt(Eb) 0];

decision_regions(v,2*sqrt(Eb),1,'PSK Constellation','Phi1(t)');

%It is all