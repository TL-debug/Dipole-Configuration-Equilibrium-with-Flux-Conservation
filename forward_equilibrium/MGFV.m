function [E0D] = MGFV(MRI,res,e0,NIT,ome)%仅返回迭代后值
% MRI 第一列弛豫矩阵 第二列插值矩阵 第三和第四列表示网格
% ---- V Cycle 测试ok 
EV = {};%用于记录迭代后的结果作为初值
% - V cycle for FMG 返回迭代后值 
% 第一次迭代 的初值 是输入值（E0 RHS）
% 输出V cycle 后的迭代值
for k = 1:2*size(MRI,1) - 1
    if ( k <= size(MRI,1) )   % 网格向下（粗）% 保留迭代结果
        if ( k == 1 ) % 首次迭代初始由函数输入
            RHS = cell2mat(res);E0 = e0;
        else 
            tmpMR = cell2mat(MRI(k,1));
            RHS = tmpMR*RES(:);RHS = reshape(RHS,sqrt(size(tmpMR,1)),[]);
            E0  = 0*RHS;
        end
        MR  = cell2mat(MRI(k,3));  
        MZ  = cell2mat(MRI(k,4)); 
        [E1,RES] = GSSOR2M(E0,RHS,MR,MZ,NIT,ome);
        EV = [EV;{E1}];
        if ( k == size(MRI,1) )
            E0D = E1;
        end
    elseif ( k > size(MRI,1) )% 网格向上（细）%  迭代结果向上传递
        MI  = cell2mat(MRI(2*size(MRI,1) - k + 1,2));
        RHS = MI*RES(:); RHS = reshape(RHS,sqrt(size(MI,1)),[]);
        
        E0  = reshape(MI*E0D(:),sqrt(size(MI,1)),[]) +...
            cell2mat(EV(2*size(MRI,1) - k));
        
        MR  = cell2mat(MRI(2*size(MRI,1) - k,3));
        MZ  = cell2mat(MRI(2*size(MRI,1) - k,4));
        [E1,RES] = GSSOR2M(E0,RHS,MR,MZ,NIT,ome);
        E0D = E1;
    end 
end
end

