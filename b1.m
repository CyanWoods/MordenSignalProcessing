%% 最小二乘法
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

%% 求自相关函数 并绘制图像
Rx=xcorr(x);                                    %求x的自相关序列Rx
subplot(1,2,1);
plot(Rx);
title("输入信号的自相关函数");
grid;

%% 设定矩阵 
p=1000;                                         %指定AR的阶数p=1000
for i=1:p
   R(i,:)=Rx(N+i-p:N+i-1);                      %按行将求得的自相关函数导入到矩阵中，获得一个p阶方阵R
end
for i=1:p
   r(:,i)=Rx(N+i);                              %将自相关函数导入一个p列1行的行向量中
end
a=-(R'*R)^(-1)*R'*r';                           %使用最小二乘法求解参数a —— '表示转置 a是一个列向量（100x1）顺序为ap到a1

%%
[h,f]=freqz(1,[1 fliplr(a')],N,Fs);             %利用freqz求出频响 所以对a的转置进行左右易位 按照a1到ap顺序排列
PSD=10*log(abs(h).^2);                          %转换为分贝形式
subplot(1,2,2);
plot(f,PSD);
title('p为100情况下的最小二乘法');
xlabel('频率(Hz)');
ylabel('幅值(dB)');
grid;

%% 寻找f1 f2 f3的估计值
[pks,locs]=findpeaks(PSD,'SortStr','descend');  %寻找峰值 并依次降序排列 
locs=locs(1:3);
locs=sort(locs);

for i=1:3                                       %显示最大值所在的坐标（即fhat）
    fprintf("F%d=%f\n",i,locs(i)/N/2)
end
