function [Y,k] = Quegan_PolCal_Y_k_FromTri(u,v,w,z,alpha,O)
% 基于三面角反射器数据，求解 Quegan参数 Y 和 k。
% 注：  Quegan原始文献中的下标是“前收 后发”；
%       而室里的下标都是“前发 后收”；
%       该函数已经统一为室里的表达；
%       因此注意下述公式和原始文献的对应情况。
%
% 输入：
%   1）已定标参数：u,v,w,z,alpha
%   2）三面角反射器的观测矩阵 O = [ O11,O12,O21,O22 ];% 这是“前发 后收”
%
% 输出：
%   定标参数 Y 和 k ;
%
% 本程序更新至：2017.12.06. 17:15


%%
% 注：由于我的输入是按照室里的“前发 后收”表示的，其对应关系刚好Quegan要求一致，无需更改 ！

A = [   alpha,      v+alpha*w,      v*w;
        alpha*u,    alpha,          v;
        alpha*z,    1,              w;
        alpha*u*z,  u+alpha*z,      1   ];

S_O = (A'*A)\(A')*(O.');    % 也即 S_O = inv(A'*A)*(A')*(O.');
                            % 这里根据 MATLAB 提示将 INV(A)*b 修改为 A\b;  
% S_O 理论上等于下式
%       S_O = Y.*[  k^2*S11;
%                   k*S21;
%                   S22     ];

%%
% 下面利用 S_O 求解 Y 和 k 的大小
%   1）求解 Y
Y = S_O(3);

%   2）求解 k
arg_k = angle( S_O(1)/Y )/2;                    % k 的相位（会有180°的模糊）
abs_k = ( abs(S_O(1))^2/(abs(Y)^2) )^(1/4);  	% k 的幅度
k = abs_k*exp(1j*arg_k);


%%
% 至此，Y 和 k 的数值已求解得到。


end