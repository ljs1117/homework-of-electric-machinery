%% 1851754 李玖思 电机大作业

%清除工作空间
clc;
clear;
close all;

%定义变量
syms id iq te

%% 100kW PMSM电机参数

Rs= 8.3e-3;           %定子相电阻(Ω)
fluxLinkage_f=0.071;  %永磁磁链峰值（Wb）
pn=4;                 %电机极对数
J= 0.1;               %电机转动惯量（kg·m^2）
L_leakage=30e-6;      %定子绕组相漏电感（H）
Lmd=174e-6;            %d轴励磁电感（H）
Lmq=293e-6;            %q轴励磁电感（H）
n=4700;               %电机额定转速（rpm）
temax=256;            %峰值转矩（Nm）

Ld = Lmd+L_leakage;   %d轴同步电感(H)
Lq = Lmq+L_leakage;   %q轴同步电感(H)

%假设电机的最高转速
nmax=12000;

%% 基本公式

%恒转矩公式（5-61）
Equ.const_torque(id,iq,te) = iq == 2*te/3/pn/(fluxLinkage_f+(Ld-Lq)*id);
%MTPA曲线（5-64）
Equ.MTPA = iq == sqrt(id*id+id*fluxLinkage_f/(Ld-Lq));
%转矩公式（5-53）
Equ.te(id,iq) = 3/2*pn*(fluxLinkage_f * iq + (Ld-Lq)*id*iq);
% 电流极限圆（5-67）
Equ.IsCircle(id,iq) = id^2 + iq^2;
%电压极限椭圆(5-75)
Equ.Vellipse(id,iq)  = Ld^2 *(id+fluxLinkage_f/Ld)^2 + Lq^2 * iq^2;
% MTPV曲线（5-89）
Equ.MTPV = iq== sqrt(Ld*(Ld*id+fluxLinkage_f)*((Ld-Lq)*id+fluxLinkage_f)/(Ld-Lq)/Lq^2);

%% （1） 峰值转矩对应峰值相电流

%联立（5-61）和（5-64）求解峰值转矩对应的峰值相电流（设定子供电电流矢量峰值为500A，以便联立数值求解时设寻找最值）
[idmax,iqmax] = vpasolve(Equ.const_torque(id,iq,temax), Equ.MTPA, [id, iq],[-500,500]);
idmax = double(idmax);
iqmax = double(iqmax);
%极坐标表示
[beta, Ismax] = cart2pol(idmax,iqmax);
disp("峰值相电流="+num2str (Ismax)+"A")

%% （2）恒转矩曲线（渐近线）、等转矩线

%画出恒转矩线及其渐近线，如图5-9

figure(1);
% 恒转矩渐近线
id_asymptote = double(fluxLinkage_f /(Lq-Ld));
%限制id、iq轴范围
xyinterval = [-Ismax id_asymptote*2  -Ismax Ismax ];
%画出恒转矩线组（转矩间隔为30Nm）
fcontour(Equ.te, xyinterval, 'LevelList',linspace(10,temax,10));
hold on;
%渐近线
plot([id_asymptote,id_asymptote],[-Ismax-1,Ismax-1],'--')
grid on;
title('100kW PMSM电机恒转矩特性曲线');
xlabel('i_d/A');
ylabel('i_q/A');
%调整坐标轴至原点处
ax = gca;ax.XAxisLocation = 'origin';ax.YAxisLocation = 'origin';

%% （3）电压极限椭圆（簇）（4）MTPA曲线 （5）MTPV曲线

%画出电压极限椭圆簇，如图5-15
figure(2);

%设置求解区域
xyinterval = [-Ismax*3 Ismax*1.5  -Ismax Ismax ];
axis(xyinterval);
axis equal;
hold on;
grid on;
title('100kW PMSM电机电压极限椭圆簇'); 
xlabel('i_d/A');
ylabel('i_q/A');
%调整坐标轴至原点处
ax = gca;ax.XAxisLocation = 'origin';ax.YAxisLocation = 'origin';
%电流极限圆
equ = Equ.IsCircle==Ismax^2;
fimplicit(equ,xyinterval);

%电压极限椭圆

% 依据（5-77），计算Usm
wrb = n*2*pi/60 * pn; %基速
Usm = wrb * sqrt((fluxLinkage_f + Ld *idmax)^2 + (Lq *iqmax)^2);
% 在额定转速到最高转速的范围内选3个速度画椭圆
wrmax = nmax*2*pi/60*pn; 
wr = linspace(wrb,wrmax,3);
% 画出电压极限椭圆
for k=1:3
    fimplicit(Equ.Vellipse == Usm^2/wr(k)^2,xyinterval);
end

%画出电压极限椭圆、MTPA、MTPV，如图5-18
figure(3);

%设置画图区域
xyinterval = [-Ismax*3 Ismax*1.5  -Ismax Ismax ];
axis(xyinterval);
axis equal;
hold on;
grid on;
title('100kW PMSM电机电压极限椭圆簇+MTPA+MTPV'); 
xlabel('i_d/A');
ylabel('i_q/A');
%调整坐标轴至原点处
ax = gca;ax.XAxisLocation = 'origin';ax.YAxisLocation = 'origin';
%电流极限圆
equ = Equ.IsCircle==Ismax^2;
fimplicit(equ,xyinterval);

%电压极限椭圆

% 依据（5-77），计算Usm
wrb = n*2*pi/60 * pn; %基速
Usm = wrb * sqrt((fluxLinkage_f + Ld *idmax)^2 + (Lq *iqmax)^2);
% 在额定转速到最高转速的范围内选3个速度画椭圆
wrmax = nmax*2*pi/60*pn; 
wr = linspace(wrb,wrmax,3);
% 画出电压极限椭圆
for k=1:3
    fimplicit(Equ.Vellipse == Usm^2/wr(k)^2,xyinterval);
end
% MTPA 曲线（5-64）
xyinterval = [idmax 0  0 iqmax ];
fimplicit(Equ.MTPA,xyinterval,'r-','linewidth',2);

% MTPV曲线（5-89）
% 求解区域
id_mtpv_max = double(-fluxLinkage_f/Ld);
xyinterval = [id_mtpv_max*1.5 id_mtpv_max+100  0 iqmax ];
%画MTPV曲线
fimplicit(Equ.MTPV, xyinterval,'k-','linewidth',2);

%% （6）机械特性（外特性）

wrmax= nmax*2*pi/60 * pn; 

% 观察画出的电压极限椭圆图+MTPA+MTPV可知，Ismax>|id_mtpv_max|

% 区间II：wr1<wr<wr2

wr1=wrb;
% 联立（5-67）和（5-89）可得wr2对应的id和iq
equ1 = Equ.IsCircle == Ismax^2;
equ2 = Equ.MTPV;   
[id_wr2,iq_wr2]=vpasolve(equ1,equ2,[id iq],[idmax,iqmax]);
id_wr2 = double(id_wr2);
iq_wr2 = double(iq_wr2);
% 根据（5-75）求wr2
wr2 = double(Usm./sqrt( Equ.Vellipse(id_wr2,iq_wr2)));
% 根据（5-53）求wr2对应的te
te_wr2 = Equ.te(id_wr2,iq_wr2);
%极坐标表示电流
[beta_wr2,Is_wr2]=cart2pol(id_wr2,iq_wr2);
% 计算区间II内的id、iq、te、wr。（位于电流极限圆上）
numte = 20;%取20个点计算
id_II = linspace(id_wr2,idmax,numte);
id_II = flip(id_II);%id增大，wr2减小，故将id序列反向
iq_II = sqrt(Ismax^2-id_II.^2);
te_II = Equ.te(id_II,iq_II);
wr_II = Usm./sqrt( Equ.Vellipse(id_II,iq_II));

% 区间III：wr2<wr<wrmax

% 计算区间III内的id、iq、te、wr。（位于MTPV曲线上）
numte = 20;%取20个点计算
wr_III = linspace(wr2,wrmax,numte);
% 初始化
id_III = ones(1,numte);
iq_III = ones(1,numte);    
for k=1:numte
    equ1=Equ.Vellipse == (Usm/ wr_III(k))^2 ;
    equ2= Equ.MTPV ;
    [id_temp, iq_temp]=vpasolve(equ1,equ2, [id iq],[id_wr2 iq_wr2]);      
    id_III(k) = double(id_temp);
    iq_III(k) = double(iq_temp);
end
te_III = double(Equ.te(id_III,iq_III));

figure(4)
hold on;
grid on;
title('100kW PMSM电机外特性曲线');
xlabel('n/rpm');
ylabel('t_e,P_e,i_s,U_s');
wr = [0,wrb,wr_II,wr_III];
n  = wr*60/2/pi/pn;
te = [temax,temax,te_II,te_III];
Pe = [0,wrb*temax/1000/pn,wr_II.*te_II/1000/pn,wr_III.*te_III/1000/pn];
is = [Ismax,Ismax,sqrt(id_II.^2+iq_II.^2),sqrt(id_III.^2+iq_III.^2)];
Us = [0,sqrt(Equ.Vellipse(idmax,iqmax)).*wrb,sqrt(Equ.Vellipse(id_II,iq_II)).*wr_II,sqrt(Equ.Vellipse(id_III,iq_III)).*wr_III];

plot(n,te,'r-')
plot(n,Pe,'b-')
plot(n,is,'g-')
plot(n,Us,'k-')
legend('t_e(N·m)','P_e(kW)','i_s(A)','U_s(V)')

%% （7）额定转速、最大转矩工作点的空间矢量图

% 额定转速、最大转矩工作点的稳态矢量图，如图5-6
figure(5);
hold on;
grid on;
axis equal;
title('100kW PMSM电机额定转速、最大转矩工作点的稳态矢量图'); 
axis([-Usm Usm -Usm*0.6 Usm*0.6]);
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
xlabel('d');
ylabel('q');

%为了使电流、磁链与电压可画在同一张矢量图上，且长度适宜，需加电流和磁链比例尺
scale_psi = 800;
scale_current = 0.25;

% 永磁磁链
X=0; Y=0; U= fluxLinkage_f*scale_psi; V=0; 
V_Psi_f = quiver(X,Y,U,V,1,'k','LineWidth',2);
% Lq*iq
X=V_Psi_f.XData+V_Psi_f.UData; Y=V_Psi_f.YData+V_Psi_f.VData; U= 0; V=Lq*iqmax*scale_psi;
V_PsiLq = quiver(X,Y,U,V,1,'r','LineWidth',2);
% Ld*id
X=V_PsiLq.XData+V_PsiLq.UData; Y=V_PsiLq.YData+V_PsiLq.VData; U= Ld*idmax*scale_psi; V=0;
V_PsiLd = quiver(X,Y,U,V,1,'g','LineWidth',2);
% 定子磁链
X=0; Y=0; U=V_PsiLd.XData+V_PsiLd.UData; V=V_PsiLd.YData+V_PsiLd.VData;
V_Psi_s = quiver(X,Y,U,V,1,'b','LineWidth',2);
% iq
X=0; Y=0; U=0; V=iqmax*scale_current;
V_iq= quiver(X,Y,U,V,1,'y','LineWidth',2);
% id
X=0; Y=V_iq.VData; U= idmax*scale_current; V=0;
V_id= quiver(X,Y,U,V,1,'m','LineWidth',2);
% Is 电流矢量
X=0; Y=0; U=V_id.XData+V_id.UData; V=V_id.YData+V_id.VData;
V_Is = quiver(X,Y,U,V,1,'c','LineWidth',2);
% e0 反电动势
X=0; Y=0; U=0; V=wrb*fluxLinkage_f;
V_e0 = quiver(X,Y,U,V,1,'Color',[1/2 0 0],'LineWidth',2);
% Is*Rs 定子电阻压降
X=0; Y=V_e0.VData; U=Rs*idmax; V=Rs*iqmax;
V_RsIs = quiver(X,Y,U,V,1,'Color',[1/2 1/2 0],'LineWidth',2);
% Lq*iq*wr 
X=V_RsIs.XData+V_RsIs.UData; Y=V_RsIs.YData+V_RsIs.VData; U=-Lq*iqmax*wrb; V=0;
V_ELq = quiver(X,Y,U,V,1,'Color',[1/2 1/2 1/2],'LineWidth',2);
% Ld*id*wr 
X=V_ELq.XData+V_ELq.UData; Y=V_ELq.YData+V_ELq.VData; U=0; V=Ld * idmax * wrb;
V_ELd = quiver(X,Y,U,V,1,'Color',[2/3 0 1],'LineWidth',2);
% Us
X=0; Y=0; U=V_ELd.XData+V_ELd.UData; V=V_ELd.YData+V_ELd.VData;
V_Us = quiver(X,Y,U,V,1,'Color',[1 1/2 0],'LineWidth',2);
drawnow
legend('\psi_f','L_qi_q','L_di_d','\psi_s','i_q','i_d','i_s','e_0','R_si_s','\omega_rL_qi_q','\omega_rL_di_d','u_s')

%% （8）额定转速、最大制动转矩回馈制动空间矢量图

% 额定转速、最大制动转矩回馈制动空间矢量图，如图5-19（b）
figure(6);
hold on;
grid on;
axis equal;
title('100kW PMSM电机额定转速、最大制动转矩回馈制动空间矢量图'); 
axis([-Usm Usm -Usm*0.6 Usm*0.6]);
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
xlabel('d');
ylabel('q');

%为了使电流、磁链与电压可画在同一张矢量图上，且长度适宜，需加电流和磁链比例尺
scale_psi = 800;
scale_current = 0.25;

% 永磁磁链
X=0; Y=0; U= fluxLinkage_f*scale_psi; V=0; 
V_Psi_f_brake = quiver(X,Y,U,V,1,'k','LineWidth',2);
% Lq*iq
X=V_Psi_f_brake.XData+V_Psi_f_brake.UData; Y=V_Psi_f_brake.YData+V_Psi_f_brake.VData; U= 0; V=-Lq*iqmax*scale_psi;
V_PsiLq_brake = quiver(X,Y,U,V,1,'r','LineWidth',2);
% Ld*id
X=V_PsiLq_brake.XData+V_PsiLq_brake.UData; Y=V_PsiLq_brake.YData+V_PsiLq_brake.VData; U= Ld*idmax*scale_psi; V=0;
V_PsiLd_brake = quiver(X,Y,U,V,1,'g','LineWidth',2);
% 定子磁链
X=0; Y=0; U=V_PsiLd_brake.XData+V_PsiLd_brake.UData; V=V_PsiLd_brake.YData+V_PsiLd_brake.VData;
V_Psi_s_brake = quiver(X,Y,U,V,1,'b','LineWidth',2);
% iq
X=0; Y=0; U=0; V=-iqmax*scale_current;
V_iq_brake= quiver(X,Y,U,V,1,'y','LineWidth',2);
% id
X=0; Y=V_iq_brake.VData; U= idmax*scale_current; V=0;
V_id_brake= quiver(X,Y,U,V,1,'m','LineWidth',2);
% Is 电流矢量
X=0; Y=0; U= V_id_brake.XData+V_id_brake.UData; V=V_id_brake.YData+V_id_brake.VData;
V_Is_brake = quiver(X,Y,U,V,1,'c','LineWidth',2);
% e0 反电动势
X=0; Y=0; U=0; V=wrb*fluxLinkage_f;
V_e0_brake = quiver(X,Y,U,V,1,'Color',[1/2 0 0],'LineWidth',2);
% Is*Rs 定子电阻压降
X=0; Y=V_e0_brake.VData; U=Rs*idmax; V=-Rs*iqmax;
V_RsIs_brake = quiver(X,Y,U,V,1,'Color',[1/2 1/2 0],'LineWidth',2);
% Lq*iq*wr 
X=V_RsIs_brake.XData+V_RsIs_brake.UData; Y=V_RsIs_brake.YData+V_RsIs_brake.VData; U=Lq*iqmax*wrb; V=0;
V_ELq_brake = quiver(X,Y,U,V,1,'Color',[1/2 1/2 1/2],'LineWidth',2);
% Ld*id*wr 
X=V_ELq_brake.XData+V_ELq_brake.UData; Y=V_ELq_brake.YData+V_ELq_brake.VData; U=0; V=Ld * idmax * wrb;
V_ELd_brake = quiver(X,Y,U,V,1,'Color',[2/3 0 1],'LineWidth',2);
% Us
X=0; Y=0; U=V_ELd_brake.XData+V_ELd_brake.UData; V=V_ELd_brake.YData+V_ELd_brake.VData;
V_Us_brake = quiver(X,Y,U,V,1,'Color',[1 1/2 0],'LineWidth',2);
drawnow
legend('\psi_f','L_qi_q','L_di_d','\psi_s','i_q','i_d','i_s','e_0','R_si_s','\omega_rL_qi_q','\omega_rL_di_d','u_s')

%% （9）MTPA对应的最佳电流分配曲线

%给定20个转矩求对应的id，iq，is
numte = 20;
te =linspace(0,temax,numte);
%初始化id，iq，is，beta序列
id_mtpa=zeros(1,numte);
iq_mtpa=zeros(1,numte);
is_mtpa=zeros(1,numte);

for k=1:numte
    [id_k,iq_k] = vpasolve(Equ.const_torque(id,iq,te(k)), Equ.MTPA, [id, iq],[idmax*k/numte iqmax/2*k/numte]);
    id_mtpa(k) = double(id_k);
    iq_mtpa(k) = double(iq_k); 
    is_mtpa(k) = sqrt(id_mtpa(k)^2+iq_mtpa(k)^2);
end

% MTPA对应的最佳电流分配曲线
figure(7);
plot(is_mtpa,id_mtpa);
hold on;grid on;
plot(is_mtpa,iq_mtpa);
title('100kW PMSM电机MTPA对应的最佳电流分配曲线');
xlabel('i_s/A');
ylabel('i_q、i_d/A');
legend('i_d','i_q')
ax = gca;ax.XAxisLocation = 'origin';ax.YAxisLocation = 'origin';

%画出恒转矩线+MTPA，如图5-11
figure(8);
%限制id、iq轴范围
xyinterval = [-Ismax id_asymptote+20  0 Ismax ];
%画出恒转矩线组（转矩间隔为30Nm）
fcontour(Equ.te, xyinterval, 'LevelStep',30);
grid on;
hold on;
plot(id_mtpa,iq_mtpa,'r-','linewidth',2);
title('100kW PMSM电机恒转矩特性曲线+MTPA');
xlabel('i_d/A');
ylabel('i_q/A');
%调整坐标轴至原点处
ax = gca;ax.XAxisLocation = 'origin';ax.YAxisLocation = 'origin';

%% （10）不同电流幅值下的“矩角特性”曲线 （11）磁阻转矩/电磁转矩曲线（第一、第二象限） 

figure(9);
hold on;
title('100kW PMSM电机不同电流幅值下的矩角特性');
xlabel('\beta/(°)'),ylabel('t_e/N·m'); 
grid on;

%取10个电流幅值画图，如图5-12
numIs = 10;
Is = linspace(0,Ismax,numIs);
beta =(0:5:180)*pi/180;%5°取一个点画（平滑）      

%初始化电磁转矩、磁阻转矩、合成转矩
torque_mag = ones(numIs,length(beta)); %电磁转矩
torque_reluctance = ones(numIs,length(beta)); %磁阻转矩
te = ones(numIs,length(beta));%合成转矩

for i=1:numIs   
    id_i = Is(i)*cos(beta); 
    iq_i = Is(i)*sin(beta);
    te_i = Equ.te(id_i,iq_i);
    te(i,:) = te_i;%合成转矩
    torque_mag(i,:) = 3/2*pn*fluxLinkage_f*iq_i;%电磁转矩
    for j = 1:length(id_i)
            id_ij = id_i(j);
            iq_ij = iq_i(j);
            torque_reluctance(i,j) = double(3/2*pn*(Ld-Lq)*iq_ij.*id_ij);%磁阻转矩
    end
    %不同电流幅值下的“矩角特性”曲线
    plot(beta*180/pi,te_i);
end
% MTPA 曲线（5-64）
beta_MTPA = atan(iq_mtpa./id_mtpa);
beta_MTPA(1) = -pi/2;
te_MTPA = double(Equ.te(id_mtpa,iq_mtpa));
plot(beta_MTPA*180/pi+180,te_MTPA,'r-','linewidth',2)

%磁阻转矩/电磁转矩/合成转矩曲线，如图5-8
figure(10);
%电磁转矩
subplot(3,1,1);
hold on;
title('100kW PMSM电机不同电流幅值下的电磁转矩'); 
xlabel('\beta/(°)');
ylabel('t_e/N·m'); 
grid on;
for k=1:numIs   
    plot(beta*180/pi,torque_mag(k,:));
end
%磁阻转矩
subplot(3,1,2);
hold on;
title('100kW PMSM电机不同电流幅值下的磁阻转矩');
xlabel('\beta/(°)');
ylabel('t_e/N·m'); 
grid on;
for k=1:numIs   
    plot(beta*180/pi,torque_reluctance(k,:));
end
%合成转矩
subplot(3,1,3);
hold on;
title('100kW PMSM电机不同电流幅值下的合成转矩');
xlabel('\beta/(°)');
ylabel('t_e/N·m'); 
grid on;
for k=1:numIs   
    plot(beta*180/pi,te(k,:));
end


