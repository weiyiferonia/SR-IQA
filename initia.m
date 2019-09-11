function beta0 = initia(OS, MOS)

%Divide data points into bands according to their OS
T = 60; %number of bands
x = zeros(T,1);
y = zeros(T,1);
dev = zeros(T-1,1);

OS_min = min(OS);
OS_max = max(OS);
interval = (OS_max - OS_min) / T;

for i=1:T
    index = find(OS>OS_min + interval*(i-1) & OS<OS_min + interval*i);
    x(i,1) = mean(OS(index));
    y(i,1) = mean(MOS(index));
end


for j=1:T-1
    dev(j,1) = abs((y(j+1,1) - y(j,1))/(x(j+1,1) - x(j,1)));
end

[G, j_M]= max(dev);
y1 = 2*(rand()-0.5);
y2 = max(MOS) - min(MOS);
y3 = 4 * G / y2;
y4 = OS(j_M);
y5 = max(MOS);

%convert y to beta0
beta0 = zeros(5,1);
beta0(1) = y2;
beta0(2) = y3;
beta0(3) = y4;
beta0(4) = y1;
beta0(5) = y5 - y2/2;







