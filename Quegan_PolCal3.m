function [u,v,w,z,alpha,r_Ohh_conjOvh,r_Ovv_conjOhv,r_Ohh_conjOvv] = Quegan_PolCal3(SLC_11,SLC_12,SLC_21,SLC_22)
% Polarization Calibration
%
% 输入：
%   1）SLC_11    选定区域的 HH 通道的 SLC 数据
%   2）SLC_12               HV
%   3）SLC_21               VH
%   4）SLC_22               VV
% 注：与国际惯例保持一致，代表“前收后发”;
%
% 输出：
%   1）定标参数 u,v,w,z,alpha
%   注意：
%       这里没有对接收通道不平衡 k 和 绝对系统增益 Y 定标（因为还额外需要角反射器数据）
%       k 和 Y 由另外的程序，利用三面角反射器数据进行求解。
%   2）同极化和交叉极化通道之间的归一化交叉极化相关系数
%       r_Ohh_conjOvh 和 r_Ovv_conjOhv
% 2016.12.05. 添加 r_Ohh_conjOvv 的计算
%   3）两个同极化通道 HH 和 VV 之间的归一化相关系数
%       r_Ohh_conjOvv
%
% 本程序在计算观测矩阵O的自相关矩阵C（4×4矩阵）时，只对非零数值进行求平均
%
% 更新至：2017.11.03. 15:59


%%
% -------------------------------------------------------------------------
%                                Quegan 算法
% -------------------------------------------------------------------------
%
% 输入：某一块区域（分布目标）的全极化（四通道） SLC 数据
%
% 输出：定标参数结果，包括 u,v,w,z 和 alpha；不包括 Y 和 k。

%%
% 计算 SLC_11 非零区域的数值个数
%   注：另外三个通道的非零区域个数应该是一致的
N_nonzero = numel(find(SLC_11~=0));

%%
% 计算观测矩阵O的自相关矩阵C（4×4矩阵）
% C11 = < O11 * conj(O11) >
C11 = sum(sum(SLC_11.*conj(SLC_11)))/N_nonzero;
% C12 = < O11 * conj(O21) >
C12 = sum(sum(SLC_11.*conj(SLC_21)))/N_nonzero;
% C13 = < O11 * conj(O12) >
C13 = sum(sum(SLC_11.*conj(SLC_12)))/N_nonzero;
% C14 = < O11 * conj(O22) >
C14 = sum(sum(SLC_11.*conj(SLC_22)))/N_nonzero;

% C21 = < O21 * conj(O11) >
C21 = sum(sum(SLC_21.*conj(SLC_11)))/N_nonzero;
% C22 = < O21 * conj(O21) >
C22 = sum(sum(SLC_21.*conj(SLC_21)))/N_nonzero;
% C23 = < O21 * conj(O12) >
C23 = sum(sum(SLC_21.*conj(SLC_12)))/N_nonzero;
% C24 = < O21 * conj(O22) >
C24 = sum(sum(SLC_21.*conj(SLC_22)))/N_nonzero;

% C31 = < O12 * conj(O11) >
C31 = sum(sum(SLC_12.*conj(SLC_11)))/N_nonzero;
% C32 = < O12 * conj(O21) >
C32 = sum(sum(SLC_12.*conj(SLC_21)))/N_nonzero;
% C33 = < O12 * conj(O12) >
C33 = sum(sum(SLC_12.*conj(SLC_12)))/N_nonzero;
% C34 = < O12 * conj(O22) >
C34 = sum(sum(SLC_12.*conj(SLC_22)))/N_nonzero;

% C41 = < O22 * conj(O11) >
C41 = sum(sum(SLC_22.*conj(SLC_11)))/N_nonzero;
% C42 = < O22 * conj(O21) >
C42 = sum(sum(SLC_22.*conj(SLC_21)))/N_nonzero;
% C43 = < O22 * conj(O12) >
C43 = sum(sum(SLC_22.*conj(SLC_12)))/N_nonzero;
% C44 = < O22 * conj(O22) >
C44 = sum(sum(SLC_22.*conj(SLC_22)))/N_nonzero;

% 至此，观测矩阵O的自相关矩阵C（4×4矩阵）已经计算完毕。

%%
% 计算定标参数，包括 u,v,w,z 和 alpha ；

% 1）The cross-talk ratios
C_delta = C11*C44 - abs(C14)^2;
%   a）u = r21 / r11
u = ( C44*C21 - C41*C24 )/C_delta;
%   b）v = t21 / t22
v = ( C11*C24 - C21*C14 )/C_delta;
%   c）w = r12 / r22
w = ( C11*C34 - C31*C14 )/C_delta;% 已修改（原来和 z 写反了）!
%   d）z = t12 / t11
z = ( C44*C31 - C41*C34 )/C_delta;% 已修改（原来和 w 写反了）!
% 至此，u,v,w,z 已计算得到。

% 2）The ratio of the receive and transmit channel imbalances
%      alpha = ( r22 / r11 )*( t11 / t22 )
alpha1 = ( C22 - u*C12 - v*C42 )/( C32 - z*C12 - w*C42 );
alpha2 = conj( C32 - z*C12 - w*C42 )/( C33 - conj(z)*C31 - conj(w)*C34 );
% 由 alpha1 和 alpha2 组合得到 alpha 的解
abs_alpha = (abs(alpha1*alpha2) - 1 + sqrt( (abs(alpha1*alpha2)-1)^2 + 4*abs(alpha2)^2 ))...
    /( 2*abs(alpha2) );
    % 理论上，alpha1 和 alpha2 的相位是相同的，任选一个作为 alpha 的相位即可；
    % 但本程序中，我取 alpha1 和 alpha2 的相位的平均值作为 alpha 的相位；
    % 2016.11.01.取平均我认为会出现一个问题，就是因为2π的混叠，取平均可能会出错;所以还是任选一个;
% phase_alpha = ( angle(alpha1) + angle(alpha2) )/2;  % 取平均值
phase_alpha = angle(alpha1);  % 取 alpha1 的相位
alpha = abs_alpha*exp(1j*phase_alpha);
% 至此，alpha 已计算得到。


%%
% 同极化和交叉极化散射信号之间的归一化交叉极化相关系数
% 表征同极化通道和交叉极化通道间的相关系数（Quegan算法假设其相关性为0）
%
% 参考文献：'Improvement of Polarimetric SAR Calibration based on the Quegan
% Algorithm'中式（5）；

% 1）Ohh 和 Ovh 计算
r_Ohh_conjOvh = C12/(sqrt(C11*C22));
% 2）Ovv 和 Ohv 计算
r_Ovv_conjOhv = C43/(sqrt(C44*C33));

% 2016.12.05. 添加下述输出
% 3）Ohh 和 Ovv 计算
r_Ohh_conjOvv = C14/(sqrt(C11*C44));


end