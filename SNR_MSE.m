%% 绘制直接法下SNR与MSE的关系
clc;
clearvars;

%% 设定信号
A1=5;           A2=3;           A3=4;
f1=0.1;         f2=0.298;       f3=0.3;
phi1=pi/3;      phi2=pi/4;      phi3=pi/2;                      %三个正弦波的参数

Fs=1;                                                           %采样频率定为1Hz
N=1024;                                                         %设置1024个采样点
n=0:N-1;                                                        %从0开始到1023截止共1024个点（1到1024也行）
x=A1*sin(2*pi*f1*n+phi1)+A2*sin(2*pi*f2*n+phi2)+A3*sin(2*pi*f3*n+phi3);      %其中的rand()表示返回一个和n有同样维数大小的随机正态分布数组

%% 设定初始值
M1=0;
M2=0;
M3=0;
t=0;
xxx=1000;

%% 开始循环计算MSE
for c=-25:0.1:10                                                %设定信噪比从-25dB开始，以0.1dB到10dB
   for i=1:xxx                                                  %进行xxx次循环来计算MSE
       y=awgn(x,c,'measured');                                  %为信号添加一个信噪比为c的加性高斯白噪声，measured表示先计算信号功率再添加噪声，使得信噪比为c
       temp=fft(y,N);                                               
       pw1=abs(temp).^2/N;                                      %以上两行为直接法计算信号功率谱密度
      [pks,locs] = findpeaks(pw1(1:N/2),'SortStr','descend');   %按照降序，将所有峰值保存起来
      locs=locs(1:3);                                           %取信号峰值最大的三个点（因为我设定了三个f）
      locs=sort(locs);                                          %按照横坐标升序将三个点排序
      f1hat=locs(1)/N;                                          
      f2hat=locs(2)/N; 
      f3hat=locs(3)/N;                                          %把点数归一化为频率，定义为fhat
      M1=M1+(f1-f1hat)^2; 
      M2=M2+(f2-f2hat)^2; 
      M3=M2+(f3-f3hat)^2;                                       %计算三个估计的误差平方和
   end 
   t=t+1;                                                       %计数用于表示横坐标点数
   MSE1(1,t)=M1/xxx;                                            
   MSE2(1,t)=M2/xxx;
   MSE3(1,t)=M3/xxx;                                            %对误差平方和进行平均，求给定信噪比下的MSE
   M1=0; 
   M2=0;
   M3=0;                                                        %重置计数
   X(1,t)=c;                                                    %用信噪比来表示横坐标
end 

%% 绘制图像
plot(X,MSE1,'r'); 
hold on
plot(X,MSE2,'k'); 
hold on
plot(X,MSE3,'b'); 
legend('f1','f2','f3');
title('SNR与MSE的关系'); 
xlabel('SNR(dB)'); 
ylabel('MSE');
