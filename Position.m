clc; clear; close all;
x=-40:1:40;y=zeros(length(x));
rn=[];xnt=[];ynt=[];znt=[];
for j=-40:1:40
    xt=[0; j; -1]; %Measuring point location
%     rn=[];
u0=4*pi*1e-7;
mx=10;    my=1851;   mz=29;    M=[mx;my;mz];   %Magnetic moment of magnetic source
d1=0.1;   d2=0.1;    %Baseline distance
xm=0; ym=0; zm=0;    %Magnetic source location

x1=xt(1); y1=xt(2); z1=xt(3);        
x2=x1; y2=y1; z2=z1-d2;      x3=x1; y3=y1; z3=z1+d2;      %The sensor position at the measuring point
x4=x1; y4=y1-d1; z4=z1;     x5=x1; y5=y1+d1; z5=z1;      X=[x1 x2 x3 x4 x5]; Y=[y1 y2 y3 y4 y5];  Z=[z1 z2 z3 z4 z5];
xc=zeros(1,5);yc=zeros(1,5);zc=zeros(1,5);rc=zeros(1,5);

for i=1:5
    xc(i)=X(i)-xm;  yc(i)=Y(i)-ym;  zc(i)=Z(i)-zm;   %The vector from the magnetic source to the measuring point
    rc(i)=(xc(i)^2+yc(i)^2+zc(i)^2)^0.5;
end

B=zeros(3,5);
for i=1:5
    B(:,i)=u0/(4*pi*rc(i)^5)*[2*xc(i)^2-yc(i)^2-zc(i)^2 3*xc(i)*yc(i)             3*xc(i)*zc(i);
                              3*xc(i)*yc(i)             2*yc(i)^2-xc(i)^2-zc(i)^2 3*yc(i)*zc(i);
                              3*xc(i)*zc(i)             3*yc(i)*zc(i)             2*zc(i)^2-xc(i)^2-yc(i)^2]*M;
end
Bxy=(B(1,5)-B(1,4))/(2*d1);  Bxz=(B(1,3)-B(1,2))/(2*d2);  Byy=(B(2,5)-B(2,4))/(2*d1);
Byz=(B(2,3)-B(2,2))/(2*d2);  Bzz=(B(3,3)-B(3,2))/(2*d2);
Byx=Bxy;Bzx=Bxz;Bzy=Byz;Bxx=-Bzz-Byy;
% G=[Bxx Bxy Bxz;Byx Byy Byz;Bzx Bzy Bzz]; %MGT

%0.01nT=0.1pT                 y = ROUNDN(x,n) %To the exact number n after the decimal point  
% Bxy1=roundn(Bxy,-11);Bxz1=roundn(Bxz,-11);Byy1=roundn(Byy,-11);Byz1=roundn(Byz,-11);Bzz1=roundn(Bzz,-11);   
% Byx1=roundn(Byx,-11);Bzx1=roundn(Bzx,-11);Bzy1=roundn(Bzy,-11);Bxx1=roundn(Bxx,-11);
% B(:,1)=roundn(B(:,1),-11);
% G=[Bxx1 Bxy1 Bxz1;Byx1 Byy1 Byz1;Bzx1 Bzy1 Bzz1];

%1pT=10^-12T=10^-3nT
Bxy1=roundn(Bxy,-12);Bxz1=roundn(Bxz,-12);Byy1=roundn(Byy,-12);Byz1=roundn(Byz,-12);Bzz1=roundn(Bzz,-12);   
Byx1=roundn(Byx,-12);Bzx1=roundn(Bzx,-12);Bzy1=roundn(Bzy,-12);Bxx1=roundn(Bxx,-12);
B(:,1)=roundn(B(:,1),-12);
G=[Bxx1 Bxy1 Bxz1;Byx1 Byy1 Byz1;Bzx1 Bzy1 Bzz1];

r=-3*(G\B(:,1)); 
rn=[rn r];
xnt=[xnt x1];ynt=[ynt y1];znt=[znt z1];


end

dnx=rn(1,:)-xnt;  dnx1=dnx';
dny=rn(2,:)-ynt;  dny1=dny';
dnz=rn(3,:)-znt;  dnz1=dnz';

qq=ones(81,1);dxmax=abs(max(dnx1)*qq);dymax=abs(max(dny1)*qq);dzmax=abs(max(dnz1)*qq);%xyz error£®max£©


xx=x';dn=[xx dnx1 dny1 dnz1 dxmax dymax dzmax];%x,y,z ERROR
yy=xx;zz=ones(81,1)*z1;
dnr=[(dnx1.^2+dny1.^2+dnz1.^2).^0.5 xx.*0 yy zz];%Total error


% xlswrite('D:\¡ı÷“œÈ\matlab data\dn_position',dn)
% xlswrite('D:\¡ı÷“œÈ\matlab data\dnr_position',dnr)

plot(x,dnx,'red');hold on
plot(x,dny,'k');hold on
plot(x,dnz);hold on