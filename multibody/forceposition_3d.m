function [qq,qq_T,xtip,ytip,ztip,jforca] = forceposition_3d(n,r,L0)
    syms q [2*n 1]
    xtip = 0;
    ytip = 0;
    ztip = 0;
    qq = 0;
    qq_T = 0;
    for i = 1:r
        qq = qq + q(i);
        qq_T = qq_T + q(n+i);
        xtip = xtip + L0*cos(qq);%-(d0-a)*sin(qq)*sin(qq_T);
        ytip = 0;%(d0-a)*cos(qq_T);
        ztip = ztip + L0*sin(qq);% +(d0-a)*sin(qq_T)*cos(qq);
    end
    jforca = transpose(jacobian([xtip;ytip;ztip],q));
end