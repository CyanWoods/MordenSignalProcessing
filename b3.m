%% Burg法求功率谱密度
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

%% 设定初始值
p=100;                                              %设定阶数p为100
fm=zeros(p,N);                                  
fm(1,:)=x;                                          %设定一个pxN的零矩阵并将第一行填充为x的数值 作为前向预测误差
gm=zeros(p,N);                                  
gm(1,:)=x;                                          %设定一个pxN的零矩阵并将第一行填充为x的数值 作为后向误差
P0=x*x'/N;                                          %预测误差功率的初始值
a=zeros(p,p);                                       %定义预测滤波器系数（二维p阶0矩阵）

%% 计算反射系数以及前向滤波器系数
k(1)=0;                                             %预定义反射系数
for m=2:p+1                                         %m从2到50开始循环（k1已经被定义过了，不需要再循环）
    x1=0;                                           %预定义变量x1
    x2=0;                                           %预定义变量x2
    for n=m:N
        x1=-1*fm(m-1,n)*gm(m-1,n-1)+x1;             %分子的表达式
        x2=0.5*(fm(m-1,n)^2+gm(m-1,n-1)^2)+x2;      %分母的表达式
    end
    K(m)=x1/x2;                                     %计算反射系数
    a(m,m)=K(m);                                    %a为从1到51的51阶对角阵 （第一行第一列为初始设定的0）
    for n=m:N
        fm(m,n)=fm(m-1,n)+K(m)*gm(m-1,n-1);         %当前阶次前向预测误差
        gm(m,n)=conj(K(m))*fm(m-1,n)+gm(m-1,n-1);   %当前阶次后向预测误差 (K(m)在这里是只有实数，所以加不加conj函数意义不大）
    end
end

%% 计算前向预测误差系数
K=K(1,2:end);                                   %去除第一列的0 转换为1x50的矩阵
a=a(2:p+1,2:p+1);                               %去除第一行第一列的0，转换为50x50的矩阵
P(1)=(1-K(1)^2)*P0;                             %计算初始误差功率P(1)
for m=2:p
    for i=1:m-1
        a(m,i)=a(m-1,i)+K(m)*conj(a(m-1,m-i));  %计算前向预测误差系数
        P(m)=(1-K(m)^2)*P(m-1);                 %预测误差功率
    end
end

%% 计算功率谱 并绘制图像
[h,f]=freqz(1,[1,a(p,:)],N,Fs);                 %求解频率响应
PSD=10*log10(abs(h).^2);                        %转换为分贝形式
plot(f,PSD);                                    %绘制图像
title('采用burg算法');
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