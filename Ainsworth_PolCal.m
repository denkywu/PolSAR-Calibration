function [u_op,v_op,w_op,z_op,alpha_op] = Ainsworth_PolCal(SLC_HH,SLC_HV,SLC_VH,SLC_VV,Flag)
% Polarization Calibration Based on Ainsworth Algorithm.
%
% 输入：
%   1）SLC_HH    选定区域的 HH 通道的 SLC 数据
%   2）SLC_HV               HV
%   3）SLC_VH               VH
%   4）SLC_VV               VV
% 注：与国际惯例保持一致，代表“前收后发”（同Ainsworth算法）;
%   5）Flag      标记：a）为0则仅输出迭代最终解；b）为1则输出迭代过程全部值（中间值和最终解）;
%                除非需要调试，则推荐选0;
%
% 输出：
%   1）定标参数 u_op,v_op,w_op,z_op,alpha_op
%   注1：
%       如果 Flag==1，则输出是向量：第1个元素表示初值，第2个至最后依次表示各迭代结果；
%       如果 Flag==0，则输出是一个值，表示迭代最终解；
%   注2：
%       这里没有对接收通道不平衡 k 和 绝对系统增益 Y 定标（因为还额外需要角反射器数据）
%       k 和 Y 由另外的程序，利用三面角反射器数据进行求解。
%
% 本程序截止至：2017.12.13. 21:19


%%
% -------------------------------------------------------------------------
%                           Ainsworth 算法
% -------------------------------------------------------------------------
%
% 输入：某一块区域（分布目标）的全极化（四通道） SLC 数据
%
% 输出：定标参数结果，包括 u,v,w,z 和 alpha；不包括 Y 和 k。


%% 计算观测矩阵O的自相关矩阵C（4×4矩阵）
C = zeros(4,4);
% C11 = < Ohh * conj(Ohh) >
C(1,1) = mean(mean(SLC_HH.*conj(SLC_HH)));
% C12 = < Ohh * conj(Ohv) >
C(1,2) = mean(mean(SLC_HH.*conj(SLC_HV)));
% C13 = < Ohh * conj(Ovh) >
C(1,3) = mean(mean(SLC_HH.*conj(SLC_VH)));
% C14 = < Ohh * conj(Ovv) >
C(1,4) = mean(mean(SLC_HH.*conj(SLC_VV)));

% C21 = < Ohv * conj(Ohh) >
C(2,1) = mean(mean(SLC_HV.*conj(SLC_HH)));
% C22 = < Ohv * conj(Ohv) >
C(2,2) = mean(mean(SLC_HV.*conj(SLC_HV)));
% C23 = < Ohv * conj(Ovh) >
C(2,3) = mean(mean(SLC_HV.*conj(SLC_VH)));
% C24 = < Ohv * conj(Ovv) >
C(2,4) = mean(mean(SLC_HV.*conj(SLC_VV)));

% C31 = < Ovh * conj(Ohh) >
C(3,1) = mean(mean(SLC_VH.*conj(SLC_HH)));
% C32 = < Ovh * conj(Ohv) >
C(3,2) = mean(mean(SLC_VH.*conj(SLC_HV)));
% C33 = < Ovh * conj(Ovh) >
C(3,3) = mean(mean(SLC_VH.*conj(SLC_VH)));
% C34 = < Ovh * conj(Ovv) >
C(3,4) = mean(mean(SLC_VH.*conj(SLC_VV)));

% C41 = < Ovv * conj(Ohh) >
C(4,1) = mean(mean(SLC_VV.*conj(SLC_HH)));
% C42 = < Ovv * conj(Ohv) >
C(4,2) = mean(mean(SLC_VV.*conj(SLC_HV)));
% C43 = < Ovv * conj(Ovh) >
C(4,3) = mean(mean(SLC_VV.*conj(SLC_VH)));
% C44 = < Ovv * conj(Ovv) >
C(4,4) = mean(mean(SLC_VV.*conj(SLC_VV)));
% 至此，观测矩阵O的自相关矩阵C（4×4矩阵）已经计算完毕。


%%
% 该脚本中，参数 u,v,w,z,alpha 表示为一个向量
% 第 1 个元素表示初值
% 第 2 个到最后依次表示各迭代结果


%% 初始化参数
% 1）串扰因子：u0 = v0 = w0 = z0 = 0
u(1) = 0;
v(1) = 0;
w(1) = 0;
z(1) = 0;
% 2）k 置为1：k0 = 1 —— 由于该脚本并不计算因子k，因此在该脚本中始终为k0.
k0 = 1;
% 3）计算 alpha 的初值 alpha0
alpha0_abs = abs(C(3,3)/C(2,2))^(1/4);
% -----------------------------------------------------------------
% alpha0_phase = atan(C(3,2)) / 2;% 原始文献中的表达式，我认为有误，改为用下式计算
alpha0_phase = angle(C(3,2)) / 2;
% -----------------------------------------------------------------
alpha(1) = alpha0_abs * exp(1j*alpha0_phase);
clear alpha0_abs;clear alpha0_phase;


%% 迭代计算定标参数 u,v,w,z,alpha
N_iter = 100;% 迭代次数
Epsilon = 1e-15;% 阈值

% 迭代结束条件：（1）达到迭代次数；或（2）alpha 的变化小于给定阈值
for i_iter = 1 : N_iter
    if 1 == i_iter  % 第1次计算时采用下式。此时的 D_Matrix(M) 由于不考虑串扰，因此简化为对角阵 G.
                    % 这是直接将 inv(G) * C * inv(G') 展开后的表达式.
        Sigma_1_Matrix = [
            C(1,1)/( abs(k0*alpha(i_iter))^2 ),                    C(1,2)*conj(alpha(i_iter))/(k0*alpha(i_iter)), C(1,3)/(k0*abs(alpha(i_iter))^2),              C(1,4)*conj(k0)*conj(alpha(i_iter))/(k0*alpha(i_iter));
            C(2,1)*alpha(i_iter)/(conj(k0)*conj(alpha(i_iter))),   C(2,2)*abs(alpha(i_iter))^2,                   C(2,3)*alpha(i_iter)/(conj(alpha(i_iter))),    C(2,4)*conj(k0)*abs(alpha(i_iter))^2;
            C(3,1)/(conj(k0)*abs(alpha(i_iter))^2),                C(3,2)*conj(alpha(i_iter))/alpha(i_iter),      C(3,3)/(abs(alpha(i_iter))^2),                 C(3,4)*conj(k0)*conj(alpha(i_iter))/alpha(i_iter);
            C(4,1)*k0*alpha(i_iter)/(conj(k0)*conj(alpha(i_iter))),C(4,2)*k0*abs(alpha(i_iter))^2,                C(4,3)*k0*alpha(i_iter)/(conj(alpha(i_iter))), C(4,4)*abs(k0*alpha(i_iter))^2;
        ];% 原始文献式(10).
    end
    
    A = ( Sigma_1_Matrix(2,1) + Sigma_1_Matrix(3,1) )/2;
    B = ( Sigma_1_Matrix(2,4) + Sigma_1_Matrix(3,4) )/2;
    
    X_Matrix = [
        Sigma_1_Matrix(2,1) - A;
        Sigma_1_Matrix(3,1) - A;
        Sigma_1_Matrix(2,4) - B;
        Sigma_1_Matrix(3,4) - B;
    ];
    Zeta_Matrix = [
        0,                      0,                      Sigma_1_Matrix(4,1),    Sigma_1_Matrix(1,1);
        Sigma_1_Matrix(1,1),    Sigma_1_Matrix(4,1),    0,                      0;
        0,                      0,                      Sigma_1_Matrix(4,4),    Sigma_1_Matrix(1,4);
        Sigma_1_Matrix(1,4),    Sigma_1_Matrix(4,4),    0,                      0;
    ];
    Tau_Matrix = [
        0,                      Sigma_1_Matrix(2,2),    Sigma_1_Matrix(2,3),    0;
        0,                      Sigma_1_Matrix(3,2),    Sigma_1_Matrix(3,3),    0;
        Sigma_1_Matrix(2,2),    0,                      0,                      Sigma_1_Matrix(2,3);
        Sigma_1_Matrix(3,2),    0,                      0,                      Sigma_1_Matrix(3,3);
    ];

    X_Matrix_Need = [
        real(X_Matrix);
        imag(X_Matrix);
    ];% 8 × 1
    Zeta_Tau_Matrix_Need = [
        real( Zeta_Matrix + Tau_Matrix ),   -1*imag( Zeta_Matrix - Tau_Matrix );
        imag( Zeta_Matrix + Tau_Matrix ),      real( Zeta_Matrix - Tau_Matrix );
    ];% 8 × 8
    
    % 求解线性方程组，得到串扰因子修正量（列向量）：
    %       线性方程组为： X_Matrix_Need = Zeta_Tau_Matrix_Need * delta_Solve;
    delta_Solve = linsolve( Zeta_Tau_Matrix_Need, X_Matrix_Need );% 利用MATLAB函数进行求解
    
    delta_u = delta_Solve(1) + 1j*delta_Solve(5);
    delta_v = delta_Solve(2) + 1j*delta_Solve(6);
    delta_w = delta_Solve(3) + 1j*delta_Solve(7);
    delta_z = delta_Solve(4) + 1j*delta_Solve(8);
    
    % 更新串扰因子 u,v,w,z
    u(i_iter+1) = u(i_iter) + delta_u;
    v(i_iter+1) = v(i_iter) + delta_v;
    w(i_iter+1) = w(i_iter) + delta_w;
    z(i_iter+1) = z(i_iter) + delta_z;
    
    % 计算 Sigma_2_Matrix
    Matrix_uvwz = [
        1,                          v(i_iter+1),                w(i_iter+1),                v(i_iter+1)*w(i_iter+1);
        z(i_iter+1),                1,                          w(i_iter+1)*z(i_iter+1),    w(i_iter+1);
        u(i_iter+1),                u(i_iter+1)*v(i_iter+1),    1,                          v(i_iter+1);
        u(i_iter+1)*z(i_iter+1),    u(i_iter+1),                z(i_iter+1),                1;
    ];
%     Sigma_2_Matrix = inv(Matrix_uvwz) * Sigma_1_Matrix * inv( Matrix_uvwz' ); % 注意最后那个是共轭转置。
    Sigma_2_Matrix = Matrix_uvwz \ Sigma_1_Matrix;
    Sigma_2_Matrix = Sigma_2_Matrix / ( Matrix_uvwz' );
    
    % 计算更新量 alpha_Updata 并更新 alpha
    alpha_Updata_abs = abs(Sigma_2_Matrix(3,3)/Sigma_2_Matrix(2,2))^(1/4);
% -----------------------------------------------------------------
%     alpha_Updata_phase = atan(Sigma_2_Matrix(3,2)) / 2;% 原始文献中的表达式，我认为有误，改为用下式计算
    alpha_Updata_phase = angle(Sigma_2_Matrix(3,2)) / 2;
% -----------------------------------------------------------------
    alpha_Updata = alpha_Updata_abs * exp(1j*alpha_Updata_phase);
    alpha(i_iter+1) = alpha(i_iter) * alpha_Updata;
    
    % 计算定标矩阵
    %   注：第一个矩阵的 alpha 是old而不是本次更新后的，详见原始文献。
    D_Matrix = diag( [alpha(i_iter), 1/alpha(i_iter), alpha(i_iter), 1/alpha(i_iter)] )...
        * Matrix_uvwz...
        * diag( [alpha_Updata, 1/alpha_Updata, alpha_Updata, 1/alpha_Updata] );
% --------------------------------------------------------------------
%     % （1） 更新观测协方差矩阵为 C = Sigma_2_Matrix;——根据张红王超的书，这里也有疑问！
%     C = Sigma_2_Matrix;——不对
    % （2）还是用假设 u=v=w=z=0 时的矩阵 G 来计算 Sigma_1_Matrix，
    % 故直接返回到本循环开头，啥也不干 —— 不对
    
    % （3）用包括 u,v,w,z 的矩阵 D_Matrix(M) 来更新 Sigma_1_Matrix.
%     Sigma_1_Matrix = inv(D_Matrix) * C * inv(D_Matrix');% 表达式
    Sigma_1_Matrix = D_Matrix \ C;
    Sigma_1_Matrix = Sigma_1_Matrix / (D_Matrix');  
% --------------------------------------------------------------------
    
    % 如果达到阈值则结束迭代
    if ( abs(alpha(i_iter+1) - alpha(i_iter)) ) < Epsilon
        break;
    end
end


%% 返回值
if Flag == 0
    % 输出最终解
    u_op = u(end);
    v_op = v(end);
    w_op = w(end);
    z_op = z(end);
    alpha_op = alpha(end);
else
    % 输出迭代过程全部解，包括初始值，中间值，最终解。
    u_op = u;
    v_op = v;
    w_op = w;
    z_op = z;
    alpha_op = alpha;
end


end