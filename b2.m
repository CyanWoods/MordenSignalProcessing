%% 整体最小二乘法
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

%% 创建增广矩阵
Rx=xcorr(x);                                    %直接调用函数求自相关序列
M=100;                                          %满足M>=pe
pe=100;                                         %设定AR阶数
for i = 1:M
    for j = 1:pe+1
        B(i,j) = Rx(N+i-j+1);                   %构造增广矩阵B（100行 101列）
    end
end

%% 利用奇异值分解求增广矩阵的有效秩
[U,K,V]=svd (B);                                %执行矩阵 B 的奇异值分解，因此 B = U*K*V'。100x100 100x101 101x101
for k=1:pe
    if(K(k,k)/K(1,1)<0.025)     
           break;  
    end                                         %找到一个可以使对角线元素之比小于0.05的k作为阈值表示有效秩（p132 4.4.40）
end  
disp(k);                                        %打印有效秩

%% 最小二乘法
p=k;
S=zeros(p+1);                                   %定义一个空矩阵
for j=1:p
    for i=1:(pe+1-p)
        djj=K(j,j)*K(j,j);
        vij=V(i:i+p,j);                         %对应课本P134(4.4.53)
        S=S+djj*(vij*vij');                     %对应课本P135(4.4.56)
    end
end                                             %构造S^p
Sinv=inv(S);                                    %求逆矩阵S^(-p)
for i=1:p
    s(1,i)=Sinv(i+1,1)/Sinv(1,1);               
end                                             %对应课本P134(4.4.58)

%% 由参数计算功率谱估计并绘图
[h,f]=freqz(1,[1,s],N,Fs);                      %求频响 —— AR系统  a0 a1 ...a7 参数为1到共8个值 N个点频率响应 
PSD=10*log10(abs(h).^2);                        %将频响转换为分贝形式
plot(f,PSD);
title('采用SVD-TLS算法');
grid;

%% 寻找f1 f2 f3的估计值
[pks,locs]=findpeaks(PSD,'SortStr','descend');  %寻找峰值 并依次降序排列 
locs=locs(1:3);
locs=sort(locs);

for i=1:3                                       %显示最大值所在的坐标（即fhat）
    fprintf("F%d=%f\n",i,locs(i)/N/2)
end
