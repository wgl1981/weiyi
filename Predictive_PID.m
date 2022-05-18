
simu_time   = 10;
ts          = 0.002;
N           = simu_time/ts;
k0          = 200;
DataSave    = 0;
dis_OP      = 1;
% LPH
ro          = 4;
t1          = 0.5;
M           = M1_M2_PID(ro,t1);
M1          = M(4:4+ro,1:3);
M2          = M(4:4+ro,4:4+ro);
g0          = [1,zeros(1,ro)];

c1          = 10;
c2          = 0.1;
c3          = 0.05;

% Geometrical Parameters
P0    = 46601.6;
S     = 0.040877;
mass  = 204.01;
vs    = 315.89;
dr    = 0.2286;
Iy    = 247.42;
Ca    = 0.3;
wa    = 150;
ita_a = 0.7;
an    = 0.000103;
bn    = -0.00945;
cn    = -0.1696;
dn    = -0.034;
am    = 0.000215;
bm    = -0.0195;
cm    = 0.051;
dm    = -0.206;
c0    = 180/pi;

Ka    = 0.7*P0*S/(mass*vs);
Kq    = 0.7*P0*S*dr/Iy;
Kz    = 0.7*P0*S/mass;
ax    = 0.7*P0*S*Ca/mass;

time        = zeros(N,1);
u           = zeros(N,1);

yd          = zeros(N,1);
x1          = zeros(N+1,1);
x1(1)       = 0;
x1_k        = x1(1);

x2          = zeros(N+1,1);
x2(1)       = 0;
x2_k        = x2(1);

x3          = zeros(N+1,1);
x3(1)       = 0;
x3_k        = x3(1);

x4          = zeros(N+1,1);
x4(1)       = 0;
x4_k        = x4(1);

e_DSC       = zeros(N+1,1);
e_DSC(1)    = 0;

Mach        = 3;
e_sum       = 0;
e_0         = 0;

for   k     = 1:N
time(k)     = k*ts;

if  time(k) <= 2
yd(k)       = 0.3;
elseif (time(k) > 2) & (time(k) <= 5)
yd(k)       = -0.2;
elseif (time(k) > 5) & (time(k) <= 8)
yd(k)       = 0.2;
else
yd(k)       = -0.1;
end

d1_yd = 0;
d2_yd = 0;
d3_yd = 0;
d4_yd = 0;

Cn_a  = an*(x1_k*c0)^3+bn*x1_k*abs(x1_k)*c0^2+cn*(2-Mach/3)*x1_k*c0;
Cm_a  = am*(x1_k*c0)^3+bm*x1_k*abs(x1_k)*c0^2+cm*(-7+8*Mach/3)*x1_k*c0+dm*c0*x3_k;

f1    = Ka*Mach*Cn_a+x2_k;
f2    = Kq*Mach^2*Cm_a;
f3    = x4_k;
f4    = -wa^2*x3_k-2*ita_a*wa*x4_k;

dz1_dx1  = 2*k0/(pi*(1+k0^2*x1_k^2));
dz2_dx1  = -4*k0^3*x1_k/(pi*(1+k0^2*x1_k^2)^2);

p1f1_px1 = Ka*Mach*(an*c0^3*3*x1_k^2+bn*c0^2*(abs(x1_k)+x1_k*sign(x1_k))+cn*(2-Mach/3)*c0);
p2f1_px1 = Ka*Mach*(an*c0^3*6*x1_k+bn*c0^2*(2*sign(x1_k)+x1_k*dz1_dx1));
p3f1_px1 = Ka*Mach*(an*c0^3*6+bn*c0^2*(3*dz1_dx1+x1_k*dz2_dx1));

p1f2_px1 = Kq*Mach^2*(am*c0^3*3*x1_k^2+bm*c0^2*(abs(x1_k)+x1_k*sign(x1_k))+cm*(-7+8*Mach/3)*c0);
p2f2_px1 = Kq*Mach^2*(am*c0^3*6*x1_k+bm*c0^2*(2*sign(x1_k)+x1_k*dz1_dx1));
p1f2_px3 = Kq*Mach^2*dm*c0;

Lf1h     = f1;
Lf2h     = p1f1_px1*f1+f2;
Lf3h     = (p1f1_px1^2+p2f1_px1*f1+p1f2_px1)*f1+p1f1_px1*f2+p1f2_px3*f3;

pLf3h_px1= p3f1_px1*f1^2+p2f1_px1*4*f1*p1f1_px1+p1f1_px1^3+p2f2_px1*f1+2*p1f1_px1*p1f2_px1+p2f1_px1*f2;
pLf3h_px2= 2*f1*p2f1_px1+p1f1_px1^2+p1f2_px1;
pLf3h_px3= p1f1_px1*p1f2_px3;
pLf3h_px4= p1f2_px3;

Lf4h     = pLf3h_px1*f1+pLf3h_px2*f2+pLf3h_px3*f3+pLf3h_px4*f4;
LgLf3h   = p1f2_px3*wa^2;

e1_k     = x1_k-yd(k);
e_sum    = e_sum+e1_k;

%[e_DSC(k+1),de_k]  = DSC(ts,0.02,e_DSC(k),e1_k);

s_k      = c1*e1_k+c2*e_sum*ts+c3*(Lf1h-d1_yd);
d1s_k    = c1*(Lf1h-d1_yd)+c2*(x1_k-yd(k))+c3*(Lf2h-d2_yd);
d2s_k    = c1*(Lf2h-d2_yd)+c2*(Lf1h-d1_yd)+c3*(Lf3h-d3_yd);

H0       = [s_k;d1s_k;d2s_k];

u(k)     = -1/(c3*LgLf3h)*(c1*(Lf3h-d3_yd)+c2*(Lf2h-d2_yd)+c3*(Lf4h-d4_yd)+g0*inv(M2)*M1*H0);

%f1_x        = Ka*Mach*(Cn_a+dn*c0*x3_k)*cos(x1_k)+x2_k;
f1_x        = f1;
f2_x        = f2;
f3_x        = f3;
f4_x        = f4;

% dis_k       = 5*sin(pi*time(k));
dis_k       = 0;
%dis_k       = 5;

x1(k+1)     = x1_k + ts*f1_x;
x2(k+1)     = x2_k + ts*f2_x;
x3(k+1)     = x3_k + ts*f3_x;
x4(k+1)     = x4_k + ts*(f4_x+wa^2*u(k)+dis_OP*dis_k);

x1_k        = x1(k+1);
x2_k        = x2(k+1);
x3_k        = x3(k+1);
x4_k        = x4(k+1);

e_0         = e1_k;
end

if DataSave  == 1
alpha = x1(1:N,1);
fin_p = u*57.29;
save time.txt            -ascii time
save yd.txt              -ascii yd
save alfa_PID.txt        -ascii alpha
save fin_PID.txt         -ascii fin_p
end

figure(1)
plot(time,yd,'r',time,x1(1:N,1),':b','linewidth',2)
xlabel('Time: s')
ylabel('\alpha: rad')

figure(2)
subplot(211)
plot(time,u*57.29,'r','linewidth',2)
xlabel('Time: s')
ylabel('u: deg')
subplot(212)
plot(time,x2(1:N,1),'r','linewidth',2)
xlabel('Time: s')
ylabel('q: rad')