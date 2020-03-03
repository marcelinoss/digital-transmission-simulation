function m = decision_regions(v,nx,ny,tstr,xl,yl);
%
%   Calc and plot the costellation
%
%
%   INPUTS
%   v = constellation points
%   v = [s11 s12 s21 s22 s31 s32 ... sn1 sn2]; 
%   nx = axis x limits
%   ny = axis y limits
%
%               +ny
%                |
%                |
%    -nx ---------------- +nx
%                |
%                |
%               -ny
%
%   tstr = string(optional)
%   xl = axis x label (optional)
%   yl = axis y label (optional)
%
%   RETURN
%   m = decision regions matrix


np = 500; % plot dimension: np x np

n = length(v);

if(mod(n,2)==0),
    inphase = v(1:2:n);
    quad = v(2:2:n);
else
    error('Number of elements in v must be even!');
end

x = linspace(-nx,nx,np);
y = linspace(-ny,ny,np);

figure;
hold on;
axis([x(1) x(length(x)) y(1) y(length(y))]);

if nargin >= 4,
    %hold on;
    title(tstr);
end

if nargin >= 5,
    %hold on;
    xlabel(xl);
end

if nargin == 6,
    %hold on;
    ylabel(yl);
end

cor = linspace(1,20,n/2);

m = zeros(length(y),length(x));

for i=1:length(x),
    for j=1:length(y),
        
        d = sqrt( (x(i) - inphase(1))^2 + (y(j) - quad(1))^2 );
        p = 1;
        
        for k=2:n/2,
            dk = sqrt( (x(i) - inphase(k))^2 + (y(j) - quad(k))^2 );
            if(dk<d)
                d = dk;
                p = k;
            end
        end
        
        m(j,i) = cor(p);
        
    end
end

image(x,y,m);

for i=1:n/2,
    plot(inphase(i),quad(i),'w*');
end

grid;