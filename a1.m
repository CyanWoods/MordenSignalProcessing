%% 直接法求功率谱估计(帕塞瓦尔定理)
clc;
clearvars;

%% 设定信号
A1=5;           A2=3;           A3=4;
f1=0.1;         f2=0.298;       f3=0.3;
phi1=pi/3;      phi2=pi/4;      phi3=pi/2;      %三个正弦波的参数

Fs=1;                                           %采样频率定为1Hz
N=1024;                                         %设置1024个采样点
n=0:N-1;                                        %从0开始到1023截止共1024个点（1到1024也行）
x=A1*sin(2*pi*f1*n+phi1)+A2*sin(2*pi*f2*n+phi2)+A3*sin(2*pi*f3*n+phi3)+1*randn(size(n));      %其中的rand()表示返回一个和n有同样维数大小的随机正态分布数组

%% 绘制信号图像
subplot(2,2,1);
stem(n,x);                                      %绘制原始信号的点图
title('原始信号x的图像')
grid;

%% 对信号进行FFT
nfft=1024;                                      %选择1024个点用来fft
xfft=fft(x,nfft);                               %对信号进行1024点fft
y=abs(xfft);                                    %取绝对值（运算结果有负有正，全部转换为正值方便计算功率谱）
spot=1:nfft/2;                                  
Y=y(spot);                                      %取一半的采样点进行分析
subplot(2,2,2);                                 
plot(spot,Y);                                   %将输出的fft表示出来
title('经FFT后取绝对值的图像')
grid;

%% 功率谱密度
P=Y.^2/nfft;                                    %功率谱密度为傅立叶变换的平方除以点数
f=spot*Fs/nfft;                                 %将横坐标从点数形式归一化为频率形式
PSD=10*log(abs(P));                             %变成分贝形式
subplot(2,2,3);
plot(f,PSD);                                    %绘制功率谱图像
title('直接法求得的功率谱图像');
xlabel('频率(Hz)');
ylabel('幅值(dB)');
grid;

%% 寻找f1 f2 f3的估计值
[pks,locs]=findpeaks(PSD,'SortStr','descend');  %寻找峰值 并依次降序排列 
locs=locs(1:3);
locs=sort(locs);                                %按照横坐标升序排列峰值点

for i=1:3                                       %显示最大值所在的坐标（即fhat）
    fprintf("F%d=%f\n",i,locs(i)/nfft)
end
