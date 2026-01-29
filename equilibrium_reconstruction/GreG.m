function [GrG] = GreG(R,Z)
% GREG 计算网格对网格的格林函数
%      网格总数对应的行；网格总数对应的列；暨使用一次矩阵运算便能得到所有网格电流在
%      线圈位置的磁通量采用周围四点的delta sigma 方法计算平均值（平均格林函数）
tmp_R = reshape(R,[],1);tmp_R = tmp_R';
tmp_Z = reshape(Z,[],1);tmp_Z = tmp_Z';
indr0 = find(tmp_R == 0);
dr    = R(2,1);dz = Z(1,2) - Z(1,1);% 步长
GrG   = zeros(length(tmp_R),length(tmp_Z));

for jg = 1:length(tmp_R)
    tmp_r = tmp_R( jg );
    tmp_z = tmp_Z( jg );

    if ( tmp_r == 0 )
        tmp_GrG = 0*tmp_r;
    else
        k2 = 4*tmp_R*tmp_r./( (tmp_r + tmp_R).^2 + (tmp_z - tmp_Z).^2 );
        [ MK,ME ] = ellipke(k2);

        tmp_GrG = sqrt( tmp_R.*tmp_r )./(2*pi*sqrt(k2)).*( (2 - k2).*MK - 2*ME );

        % 每次循环都需要对R=0的位置赋值为零 以及 无定义点（自己本身）
        % R=0的点
        tmp_GrG(indr0) = 0;
        % 无定义点
        Rrlud    = [ tmp_r - dr/16, tmp_r + dr/16, tmp_r       , tmp_r       ];
        Zrlud    = [ tmp_z       , tmp_z       , tmp_z + dz/16, tmp_z - dz/16];
        k2rlud   = 4*tmp_r*Rrlud./( ( tmp_r + Rrlud ).^2 + ( tmp_z - Zrlud ).^2 );
        [MK,ME]  = ellipke(k2rlud);
        tmpGrlud = sqrt( tmp_r.*Rrlud )./(2*pi*sqrt(k2rlud)).*( (2 - k2rlud).*MK -2*ME );
        tmp_GrG(jg) = sum( tmpGrlud )/length(tmpGrlud);
    end
    GrG(jg,:) = tmp_GrG;
end

end

