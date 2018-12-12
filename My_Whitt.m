function [T_solve_fin,R_solve_fin,A_solve_fin] = My_Whitt(S1,M1,S2,M2,S3,M3,M_test,Flag)
% 我的改进Whitt算法的 主函数
%
% 通过引入额外的反射器数据,完成完整的求解过程,筛选出所需的系统失真矩阵求解值;
%
% 输入变量：
%   1）S1 和 M1 是三面角反射器的散射矩阵和观测矩阵；
%   2）S2 和 M2 是0°二面角反射器的散射矩阵和观测矩阵；
%   3）S3 和 M3 是45°二面角反射器的散射矩阵和观测矩阵；
%   默认为以上顺序，最好不要改变；
%   4）M_test    是额外引入的角反射器数据（观测矩阵）,用作筛选参考；
%   5）Flag      是标记,表示 M_test 的类型：
%                Flag == 0，表示 22.5°二面角反射器;
%                Flag == 1，表示 三面角反射器;
%
% 输出变量：
%   1）T_solve_fin 是求解得到的发射失真矩阵T;
%   2）R_solve_fin 是求解得到的接收失真矩阵R;
%   3）A_solve_fin 是求解得到的绝对幅度因子A;
% 这是筛选后的结果：
%   a)若 Flag == 0,则可筛选出唯一解,
%           T_solve_fin 和 R_solve_fin 都是 2×2 的矩阵;
%           A_solve_fin 是1×1的数值.
%   b)若 Flag == 1,则只能筛选出两个解,
%           T_solve_fin 和 R_solve_fin 都是 4×2 的矩阵 —— 前两行表示一组解,后两行表示另一组解;
%           A_solve_fin 是 2×1 的矩阵 —— (1,1)表示第一组解对应的A, (2,1)表示另一组解对应的A;
%
% 该程序截止至：2016.10.10. 19:14


%%  Flag == 0,使用22.5°二面角反射器的筛选过程

% 1）首先考察校准结果的 VV 通道相位：应该是±180°,故如果是0°附近,则舍弃;
% 2）再考察 HV 和 VH 通道相位：应该是0°,故如果是±180°附近,则舍弃;
% 3）最终筛选出唯一解：HV 和 VH 通道相位应该是0°附近,VV 通道相位应该是±180°附近;
if Flag == 0
    for Flag1 = 0:1:1       % Flag1 分别置为0和1
        for Flag2 = 0:1:1   % Flag2 分别置为0和1
            [T_solve,R_solve,A_solve] = My_Whitt_sub(S1,M1,S2,M2,S3,M3,Flag1,Flag2);
            % 对"检测用 22.5°二面角反射器"进行极化校准
            M_After_PC = M_Polarization(reshape(M_test.',1,[]),T_solve,R_solve,A_solve);
            gyh_M_After_PC = M_After_PC./M_After_PC(1,1);   % 归一化
            
            tmp = angle(gyh_M_After_PC)/pi*180;             % 取相位（°）
            % 考察 tmp 的相位关系
            if abs(tmp(1,4)) < 90   % 如果0°附近则舍弃，这里用90°作为阈值
                continue;
            end
            if abs(tmp(1,2)) > 90 ...
                    || abs(tmp(1,3) > 90)  % 如果180°附近则舍弃，这里用90°作为阈值
                continue;
            end
            T_solve_fin = T_solve;
            R_solve_fin = R_solve;
            A_solve_fin = A_solve;
        end
    end
%     disp('-------------------------------------------------------------------')
%     disp('求解完成');
%     disp('利用22.5°二面角反射器数据，筛选得到系统失真矩阵唯一解');
%     disp('-------------------------------------------------------------------')
end


%%  Flag == 1,使用三面角反射器的筛选过程

% 1）考察校准结果的 VV 通道相位（相位不平衡）：应该是0°,故如果是180°附近,则舍弃;
% 2）最终只能筛选出两组解,对于不需要串扰相位的应用,任选一组皆可。
count = 0;
if Flag == 1
    for Flag1 = 0:1:1       % Flag1 分别置为0和1
        for Flag2 = 0:1:1   % Flag2 分别置为0和1
            [T_solve,R_solve,A_solve] = My_Whitt_sub(S1,M1,S2,M2,S3,M3,Flag1,Flag2);
            % 对"检测用三面角反射器"进行极化校准
            M_After_PC = M_Polarization(reshape(M_test.',1,[]),T_solve,R_solve,A_solve);
            gyh_M_After_PC = M_After_PC./M_After_PC(1,1);   % 归一化            
            
            tmp = angle(gyh_M_After_PC)/pi*180;             % 取相位（°）
            % 考察 tmp 的相位关系
            if abs(tmp(1,4)) > 90   % 如果180°附近则舍弃，这里用90°作为阈值
                continue;
            end
            
            if count == 0       % 第一组可能解
                T_solve_fin = T_solve;
                R_solve_fin = R_solve;
                A_solve_fin = A_solve;
                count = 1;
            end
            if count == 1      % 第二组可能解
                T_solve_fin(3:4,:) = T_solve;
                R_solve_fin(3:4,:) = R_solve;
                A_solve_fin(2,1) = A_solve;                
            end
        end
    end
    disp('-------------------------------------------------------------------')
    disp('求解完成');
    disp('利用三面角反射器数据，筛选得到系统失真矩阵两组解');
    disp('-------------------------------------------------------------------')
end


end