%% 间接法求功率谱估计(维纳欣钦关系)
clc;
clearvars;

%% 设定信号
A1=5;           A2=3;           A3=4;
f1=0.1;         f2=0.298;       f3=0.3;
phi1=pi/3;      phi2=pi/4;      phi3=pi/2;      %三个正弦波的参数

Fs=1;                                           %采样频率定为1Hz
N=1024;                                         %设置1024个采样点
n=0:N-1;                                        %从0开始到1023截止共1024个点（1到1024也行）
x=A1*sin(2*pi*f1*n+phi1)+A2*sin(2*pi*f2*n+phi2)+A3*sin(2*pi*f3*n+phi3)+10*randn(size(n));      %其中的rand()表示返回一个和n有同样维数大小的随机正态分布数组

%% 绘制信号图像
subplot(2,2,1);
stem(n,x);                                      %绘制原始信号的点图
title('原始信号x的图像')
grid;

%% 计算信号的自相关估计 并绘制点图
Rx=xcorr(x,'unbiased');                         %求其无偏自相关估计
subplot(2,2,2);
stem(Rx);                                       %绘制自相关估计的图像
title('自相关估计的点图')
grid;

%% 由自相关估计计算功率谱估计
nfft=1024;                                      %设定fft的点数为1024
xfft=fft(Rx,nfft);                              %对无偏自相关估计进行1024点fft
Pxx=abs(xfft);                                  %原fft求得的CXk存在负值，对其求绝对值
spot=1:nfft/2;                                  %index为1到512的512点(观测单边的功率谱）                                
PSD=10*log(Pxx(spot));                          %计算功率谱密度
subplot(2,2,3);
f=spot*Fs/nfft;                                 %将点位转换为对应的频率f

%% 绘制功率谱估计的图像
plot(f,PSD);
title('间接法得到的功率谱密度图像');
xlabel('频率(Hz)');
ylabel('幅值(dB)');
grid;

%% 取功率谱函数的最大值(没什么必要，故注释掉)
%{

P1=P.*(k>99&k<101);                             % 分段处理，分别求最大功率密度值
P1hat=max(P1);
disp(P1hat);

P2=P.*(k>249&k<251);
P2hat=max(P2);
disp(P2hat);

P3=P.*(k>299&k<301);
P3hat=max(P3);
disp(P3hat);

%}

%% 寻找f1 f2 f3的估计值
[pks,locs]=findpeaks(PSD,'SortStr','descend');  %寻找峰值 并依次降序排列 
locs=locs(1:3);
locs=sort(locs);

for i=1:3                                       %显示最大值所在的坐标（即fhat）
    fprintf("F%d=%f\n",i,locs(i)/nfft)
end