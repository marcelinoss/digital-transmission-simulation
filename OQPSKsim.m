%   Simulating a OQPSK transmission
%   Nb bits are randomly generated and modulated to OQPSK
%   The random bits, the symbols and the OQPSK signal waveform are plotted
%   Also, the QPSK signal waveform is plotted to be compare with OQPSK
%   By Marcelino Silva, March 2nd, 2020

%Optional
clear all;
close all;

%Parameters
Tb = 1e-3;  %Bit duration (s)
T = 2*Tb;   %Symbol duration
t = [T/100:T/100:T]; %Time interval to plot one waveform period

nc = 2;     %Cycles number
fc = nc/T;  %Frequency

Eb = 1e-6;  %Bit energy
E = Eb*2;   %Symbol energy

Nb = 100;   %Number of bits to be transmitted
NM = Nb/2;  %Number of symbols to be transmitted

J = 6;   %Number of symbols to plot the signal waveform (symbols being transmitted)

%Base functions
phi1t = sqrt(2/T).*cos(2*pi*fc.*t);
phi2t = sqrt(2/T).*sin(2*pi*fc.*t);

%Signals: s1(t), s2(t), s3(t) and s4(t) for QPSK
i=1; s1t = +sqrt(E)*cos((2*i-1)*pi/4)*phi1t -sqrt(E)*sin((2*i-1)*pi/4)*phi2t;
i=2; s2t = +sqrt(E)*cos((2*i-1)*pi/4)*phi1t -sqrt(E)*sin((2*i-1)*pi/4)*phi2t;
i=3; s3t = +sqrt(E)*cos((2*i-1)*pi/4)*phi1t -sqrt(E)*sin((2*i-1)*pi/4)*phi2t;
i=4; s4t = +sqrt(E)*cos((2*i-1)*pi/4)*phi1t -sqrt(E)*sin((2*i-1)*pi/4)*phi2t;

%Random bits to be transmitted
bits = randn(Nb,1);
bits(find(bits<=0)) = -1; % -1 means bit 0
bits(find(bits>0)) = +1; % +1 means bit 1

%Symbols Codification
m = zeros(1,NM);
for i=1:2:Nb,
    if(bits(i)==+1),
        if(bits(i+1)==+1)  %symbol m4 = [1 1]
            m((i+1)/2) = 4;
        else            %symbo ml = [1 0]
            m((i+1)/2) = 1;
        end
    else
        if(bits(i+1)==+1)  %symbol m3 = [0 1]
            m((i+1)/2) = 3;
        else            %symbol m2 = [0 0]
            m((i+1)/2) = 2;
        end
    end
end

%offset
lt = length(t);
delay = lt/2;

%QPSK Signal
s = zeros(1,lt*NM);
for i=1:NM,
    if (m(i)==1)
        s(lt*(i-1)+1:lt*i) = s1t;
    elseif (m(i)==2)
        s(lt*(i-1)+1:lt*i) = s2t;
    elseif (m(i)==3)
        s(lt*(i-1)+1:lt*i) = s3t;
    elseif (m(i)==4)
        s(lt*(i-1)+1:lt*i) = s4t;
    end
end

%OQPSK signal
soff = zeros(1,lt*NM);
OFFSET = zeros(1,delay);

for i=1:NM,
    
    if (m(i)==1)
        inp  = +sqrt(E/2)*phi1t;
        quad = -sqrt(E/2)*phi2t;
        soff(lt*(i-1)+1:lt*i) = inp + [OFFSET quad(1:delay)];
        
    elseif (m(i)==2)
        inp  = -sqrt(E/2)*phi1t;
        quad = -sqrt(E/2)*phi2t;
        soff(lt*(i-1)+1:lt*i) = inp + [OFFSET quad(1:delay)];
    
    elseif (m(i)==3)
        inp  = -sqrt(E/2)*phi1t;
        quad = +sqrt(E/2)*phi2t;
        soff(lt*(i-1)+1:lt*i) = inp + [OFFSET quad(1:delay)];
        
    elseif (m(i)==4)
        inp  = +sqrt(E/2)*phi1t;
        quad = +sqrt(E/2)*phi2t;
        soff(lt*(i-1)+1:lt*i) = inp + [OFFSET quad(1:delay)];
        
    end
    
    OFFSET = quad(delay+1:lt);
end

%Plotting the simulation

%Index of symbols being trasmitted
beginM = 1;
endM= J;

%Index of bits being trasmitted
beginB = beginM*2-1;
endB = endM*2;

%Index of signal being trasmitted
beginS = 1;
endS = endM*lt;

%plotting bits
subplot(4,1,1);
hold on;
axis([1 Nb -1.1 1.1]);
bar((1:(beginB)),bits(1:beginB),'b');
bar((beginB:endB),bits(beginB:endB),'r');
bar((endB:Nb),bits(endB:Nb),'k');
title('Bits');
legend('Transmitted Bits','Bits being transmitted','Bits to be transmitted');

%plotting symbols
subplot(4,1,2);
hold on;
axis([1 NM 0 4.1]);
bar((1:(beginM)),m(1:beginM),'b');
bar((beginM:endM),m(beginM:endM),'r');
bar((endM:NM),m(endM:NM),'k');
title('Symbols');
legend('Transmitted Symbols','Symbols being transmitted','Symbols to be transmitted')

%plotting QPSK signal
subplot(4,1,3);
hold on;
title('QPSK Signal being transmitted');
xlabel('Time (s)');
ylabel('Signal');
    
%plotting QPSK signal
subplot(4,1,4);
hold on;
title('OQPSK Signal being transmitted');
xlabel('Time (s)');
ylabel('Signal');

while(endM<=NM),
    
    %plotting bits  
    subplot(4,1,1);
    bar((1:beginB),bits(1:beginB),'b'); %Transmitted Bits
    bar((beginB:endB),bits(beginB:endB),'r'); %Bits being transmitted 
    bar((endB:Nb),bits(endB:Nb),'k'); %Bits to be transmitted
    
    %plotting symbols  
    subplot(4,1,2);
    bar((1:beginM),m(1:beginM),'b'); %Transmitted Symbols
    bar((beginM:endM),m(beginM:endM),'r'); %Symbols being transmitted
    bar((endM:NM),m(endM:NM),'k'); %Symbols to be transmitted
    
    %plotting QPSK signal
    subplot(4,1,3);
    ts = [(T/100)+(beginM-1)*T:T/100:endM*T];
    plot(ts,s(beginS:endS),'r');
    axis([((T/100)+(beginM-1)*T) endM*T -1*sqrt(2*E/T) sqrt(2*E/T)]);
    
    %plotting OQPSK signal
    subplot(4,1,4);
    %ts = [(T/100)+(beginM-1)*T:T/100:endM*T];
    plot(ts,soff(beginS:endS),'r');
    axis([((T/100)+(beginM-1)*T) endM*T -1*sqrt(2*E/T) sqrt(2*E/T)]);
    
    %update indexes
    beginM = beginM + 1;
    endM = endM + 1;
    
    beginB = beginM*2-1;
    endB = endM*2;
    
    beginS = beginS + lt;
    endS = endS + lt;
    
    pause(0.7);
    
end;

%Plotting the costellation and the decision regions

v = [+sqrt(E)*cos((2*1-1)*pi/4) -sqrt(E)*sin((2*1-1)*pi/4) +sqrt(E)*cos((2*2-1)*pi/4) -sqrt(E)*sin((2*2-1)*pi/4) ...
    +sqrt(E)*cos((2*3-1)*pi/4) -sqrt(E)*sin((2*3-1)*pi/4) +sqrt(E)*cos((2*4-1)*pi/4) -sqrt(E)*sin((2*4-1)*pi/4)];

decision_regions(v,2*sqrt(E/2),2*sqrt(E/2),'OQPSK Constellation','Phi1(t)','Phi2(t)');

%It is all