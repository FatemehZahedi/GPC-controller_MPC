clear
clc
[n1,d1,n2,d2]=Inputsys(1);
Gs1 = tf(n1,d1);
Ts=0.1;
Gd1 = c2d(Gs1,Ts,'impulse');
[num1,den1]=tfdata(Gd1,'v');
Gs2 = tf(n2,d2);
Gd2 = c2d(Gs2,Ts,'impulse');
[num2,den2]=tfdata(Gd2,'v');
sys_info = stepinfo(Gd1);
ts1 = sys_info.SettlingTime;
tr1=sys_info.RiseTime; 
sys_info = stepinfo(Gd2);
ts2 = sys_info.SettlingTime;
tr2=sys_info.RiseTime;
t=1:Ts:10;
[g1,t1] = impulse(Gd1,t);
[g2,t2] = impulse(Gd2,t);
P1=floor(tr1/Ts);
P2=floor(tr2/Ts);
N1=floor( ts1/Ts);
N2=floor( ts2/Ts);
P=max(P1,P2);
N=max(N1,N2);
M=P;
%.....................Toeplitz Matrix.................................
b1 = zeros(1,P); b1(1,1)= g1(2);
a1 = g1(2:P+1);
G1 = toeplitz(a1,b1);
G1(:,M) = G1(:,M:P)*ones(P-M+1,1);
G1 = G1(:,1:M);
%........................................................
b2 = zeros(1,P); b2(1,1)= g2(2);
a2 = g2(2:P+1);
G2 = toeplitz(a2,b2);
G2(:,M) = G2(:,M:P)*ones(P-M+1,1);
G2 = G2(:,1:M);
G=[G1 G2];
%........................................................................................
%A~=1-2.564z^-1+2.2365z^-2-0.6725z^-3
%According to the discrete transfer function define below parameters
na=3;
nb1=1; nb2=1;
nb=nb1;
d=0;
N1=d+1;
N2=d+P;
a_=[1 -2.564 2.2365 -0.6725];
b1_=num1(2:end);
b2_=num2(2:end);
C=1;  % because of using white noise
f=zeros(P+d,na+1);
f(1,1:3)=-1*a_(2:4);
for j=1:P+d-1
    for i=1:na
        f(j+1,i)=f(j,i+1)-f(j,1)*a_(i+1);
    end
end
F=f(N1:N2,1:na);
%.......................................
E1=zeros(P);
E1(:,1)=1;
for j=1:P-1
    E1(j+1:P,j+1)=f(j,1);
end
B1=zeros(P,P+nb);
for k=1:P
        B1(k,k:k+1)=b1_;
end
m1_=E1*B1;
M1_=zeros(P,nb+d);
for k=1:P
    M1_(k,:)=m1_(k,k+1);
end
%............................
E2=zeros(P);
E2(:,1)=1;
for j=1:P-1
    E2(j+1:P,j+1)=f(j,1);
end
B2=zeros(P,P+nb);
for k=1:P
        B2(k,k:k+1)=b2_;
end
m2_=E2*B2;
M2_=zeros(P,nb+d);
for k=1:P
    M2_(k,:)=m2_(k,k+1);
end
M_=[M1_ M2_];
%...............................................................................
gamma =0.001;
gain_DC=(num1(1)+num1(2)+num1(3))/(den1(1)+den1(2)+den1(3));
gain_DC2=(num2(1)+num2(2)+num2(3))/(den2(1)+den2(2)+den2(3));
Q = eye(P);
R1 =((1.4)^2)*gamma*gain_DC^2*eye(M);
R2=gamma*gain_DC2^2*eye(M);
R=[R1 zeros(M); zeros(M) R2];
alpha=0.5;
Kgpc=(G'*Q*G+R)\(G'*Q);
x01=0.0882;
x02=441.2;
%...................................................................................
U1_ = zeros(nb+d,length(t));
U2_ = zeros(nb+d,length(t));
U_=[U1_ ; U2_];
d=zeros(1,length(t));
%y1=0; %linear
y1=441.2;
u_1=[];
u_2=[];
ym=[];
y=0;
Y_d=zeros(P,length(t));
Y_past=zeros(P,length(t));
Y_m=zeros(P,length(t));
D=zeros(P,length(t));
E_bar=zeros(P,length(t));
U1=zeros(M,length(t));
U2=zeros(M,length(t));
U=[U1;U2];
Y_=zeros(na,length(t));
%..................step...........................
r =ones(length(t),1);
%...................sine..............................
%[r,t1]= gensig('sine',length(t)*Ts/2,length(t)*Ts,Ts);

for i=1:length(t)-1 
    
for j=1:P
  Y_d(j,i)=(alpha^j)*y+(1-(alpha)^j)*r(i); % Programmed
end 
    
Y_past(:,i)=M_*U_(:,i)+F*Y_(:,i);
D(:,i)=d(i)*ones(P,1);
    
E_bar(:,i)=Y_d(:,i)-Y_past(:,i)-D(:,i);
U(:,i)=Kgpc*E_bar(:,i);
U1(:,i)=U(1:M,i);
U2(:,i)=U(M+1:2*M,i);
U(:,i)=[U1(:,i);U2(:,i)];
    
Y_m(:,i)=G*U(:,i)+Y_past(:,i);


U1_(2:nb+d,i+1) = U1_(1:nb+d-1,i);
U1_(1,i+1)=U1(1,i);
U2_(2:nb+d,i+1) = U2_(1:nb+d-1,i);
U2_(1,i+1)=U2(1,i);
U_(:,i+1)=[U1_(:,i+1);U2_(:,i+1)];

Y_(2:na,i+1)=Y_(1:na-1,i);
Y_(1,i+1)=Y_m(1,i);

u1=U1(1,i);
u2=U2(1,i);
sim('Model')
%d(i+1)=yl(end)-Y_m(1,i); %linear
d(i+1)=y(end)-Y_m(1,i);
%y1=[y1;yl(end)];  % linear
y1=[y1; y(end)+441.2];
x01=x1(end);
x02=x2(end);
%y=yl(end); % linear
y=y(end);    % nonlinear
ym=[ym; Y_m(1,i)];
u_1=[u_1; u1];
u_2=[u_2; u2];
end
figure(3);
plot(y1,'b');
hold on
plot(r+441.2,'r');
 grid on
%axis([0 45 439 447]);
legend('y','r');
title('Response of the nonlinear system');
xlabel('sample');
figure(4);
plot(y1-441.2,'b');
hold on
plot(ym,'r');
grid on
xlabel('sample');
title('Ym and Yp without bias');
legend('YPlant','YModel');
figure(5);
plot(u_1,'b');
grid on
xlabel('sample');
title('Control law for input 1 without bias');
figure(6);
plot(u_2,'b');
grid on
xlabel('sample');
title('Control law for input 2 without bias');