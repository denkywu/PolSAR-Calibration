function [T,R,S1,S2,S3,S_test,M1,M2,M3,M_test] = Gen_Data_PARC(delta_dB,f_dB,f_phase,SNR_dB)
% 函数：生成 PARC 算法（基于有源定标器的点目标极化定标算法）所需数据
% 
% 输入变量：
%   1）delta_dB      设置 R 和 T 失真矩阵的串扰大小,单位 dB（串扰相位取随机相位）
%   2）f_dB          设置 R 和 T 失真矩阵的幅度不平衡,单位 dB
%   3）f_phase       设置 R 和 T 失真矩阵的相位不平衡,单位 度（°）——可以指定，也可以随机生成
%   4）SNR_dB        设置系统模型中加性噪声水平,表示比绝对幅度A低 SNR_dB（信噪比）,单位 dB
% 输出变量：
%   1）T,R               分别表示失真矩阵的设置值（真值）
%   2）S1,S2,S3          分别表示选用的三个 PARC 的散射矩阵
%   3）M1,M2,M3          分别对应 S1,S2,S3 的观测矩阵
%   4）S_test 和 M_test  分别对应检测用 三面角反射器的散射矩阵和观测矩阵
%
%
% 生成极化定标所需要的几个矩阵；
%   1）3个 PARC 的理论散射矩阵 S，和1个三面角反射器的 S;
%   2）对应的目标观测矩阵 M;
% 采用的系统模型如下
%       M = A*exp(1j*phy)*RST + N
% 设置：
%   1）A 设置为 1，绝对幅度因子
%   2）phy 从生成的 [-π,π] 随机数中选取
%   3）R 和 T 的设置见下文
%   4) N 根据 SNR 设置或忽略
%
% 本程序截止至：2017.11.16. 19:23


%%
% 生成 8 个区间为 [-π,π] 的随机数，用来作为：
%   1）发射、接收矩阵中串扰的相位
%      发射、接收矩阵中各2个，总共需要 4 个随机数作为相位；
%   2）不同定标点目标所对应的绝对相位 phy（见系统模型中的 phy）
%      PARC算法需要 3 个不同的参考定标器，再加上 1 个检测用三面角反射器，
%      总共需要 4 个随机数作为绝对相位；

phy_rand_Num = 8; % 随机相位的个数

Phase_TR = rand(1,phy_rand_Num)*(2*pi) - pi;    % [-π,π] 中的随机数
% save('Phase_TR.mat','Phase_TR');
% load Phase_TR;


%%
% 设置系统失真矩阵

% 1） 接收失真矩阵 R
%     R = [ 1    , delta2 
%          delta1,   f1   ];
delta2_dB = delta_dB;    	% delta2，串扰，单位 dB
delta1_dB = delta_dB;      	% delta1，串扰，单位 dB
f1_dB = f_dB;           	% f1，幅度不平衡，单位 dB
% --------------------------------------------------------
% 这是随机生成相位不平衡数值；
% 如果需要指定，则注销该行即可。
f_phase = rand(1,1)*(360);  % [0, 360]（°）之间随机数
% --------------------------------------------------------
f1_phase = f_phase;       	% f1，相位不平衡，单位 度（°）

R = [           1,                              10^(delta2_dB/20)*exp(1j*Phase_TR(1,1));
    10^(delta1_dB/20)*exp(1j*Phase_TR(1,2)),    10^(f1_dB/20)*exp(1j*(f1_phase/180*pi))];

% 2） 发射失真矩阵 T
%     T = [ 1    , delta3
%          delta4,   f2   ];
%   2016.08.11. 假设发射和接收失真矩阵具有互易性（幅度）;相位依然随机
delta3_dB = delta1_dB;     	% delta3，串扰，单位 dB
delta4_dB = delta2_dB;  	% delta4，串扰，单位 dB
f2_dB = f1_dB;            	% f2，幅度不平衡，单位 dB
f2_phase = f1_phase;      	% f2，相位不平衡，单位 度（°）

T = [           1,                              10^(delta3_dB/20)*exp(1j*Phase_TR(1,3))
    10^(delta4_dB/20)*exp(1j*Phase_TR(1,4)),    10^(f2_dB/20)*exp(1j*(f2_phase/180*pi))];


%%
% 设置不同定标器的理论散射矩阵

% 1） PARC-1
S1 =    [ 0, 0
          1, 0 ];
% 2） PARC-2
S2 =    [ 0, 1
          0, 0 ];
% 3）PARC-3
S3 =    [ -1, -1
           1,  1 ];

% 4）检测用“三面角反射器”
S_test =    [ 1, 0
              0, 1 ];


%%
% 根据系统失真矩阵和设定的理论散射矩阵，计算得到对应的观测矩阵
A = 1;              	% 设置，绝对幅度因子

% 1） PARC-1
% M1 = A*exp(1j*Phase_TR(1,5))*R*S1*T;
M1 = A*exp(1j*Phase_TR(1,5))*R*S1*T + A*10^(-SNR_dB/20)*exp(1j.*(rand(2,2)*2*pi-pi));

% 2） PARC-2
% M2 = A*exp(1j*Phase_TR(1,6))*R*S2*T;
M2 = A*exp(1j*Phase_TR(1,6))*R*S2*T + A*10^(-SNR_dB/20)*exp(1j.*(rand(2,2)*2*pi-pi));

% 3）PARC-3
% M3 = A*exp(1j*Phase_TR(1,7))*R*S3*T;
M3 = A*exp(1j*Phase_TR(1,7))*R*S3*T + A*10^(-SNR_dB/20)*exp(1j.*(rand(2,2)*2*pi-pi));

% 4）检测用“三面角反射器”
% M_test = A*exp(1j*Phase_TR(1,8))*R*S_test*T;
M_test = A*exp(1j*Phase_TR(1,8))*R*S_test*T + A*10^(-SNR_dB/20)*exp(1j.*(rand(2,2)*2*pi-pi));


%%
% 至此,仿真数据完成


end





