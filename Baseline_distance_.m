clc; clear; close all;
rn=[];xnt=[];ynt=[];znt=[];rnn=[];sizexy=20;
c1=0.05:0.005:0.1;c2=0.05:0.05:0.55;dd=length(c1)^2;
for k=1:1:4
x1=5*k; y1=0; z1=-1; %Measuring point location
for d1=0.05:0.005:0.1        %Baseline distance
%    rnn=[rn;rnn];  
for d2=0.05:0.05:0.55     %Baseline distance

    
u0=4*pi*1e-7;
mx=10;    my=1851;   mz=29;    M=[mx;my;mz];   %Magnetic moment of magnetic source
xm=0; ym=0; zm=0;    %Magnetic source location

          
x2=x1; y2=y1; z2=z1-d2;      x3=x1; y3=y1; z3=z1+d2;      %The sensor position at the measuring point
x4=x1; y4=y1-d1; z4=z1;     x5=x1; y5=y1+d1; z5=z1;      X=[x1 x2 x3 x4 x5]; Y=[y1 y2 y3 y4 y5];  Z=[z1 z2 z3 z4 z5];
xc=zeros(1,5);yc=zeros(1,5);zc=zeros(1,5);rc=zeros(1,5);

for i=1:5                             
    xc(i)=X(i)-xm;  yc(i)=Y(i)-ym;  zc(i)=Z(i)-zm;   %The vector from the magnetic source to the measuring point
    rc(i)=(xc(i)^2+yc(i)^2+zc(i)^2)^0.5;  %scalar
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

% 0.1pT                 y = ROUNDN(x,n) %To the exact number n after the decimal point
% Bxy1=roundn(Bxy,-13);Bxz1=roundn(Bxz,-13);Byy1=roundn(Byy,-13);Byz1=roundn(Byz,-13);Bzz1=roundn(Bzz,-13);   
% Byx1=roundn(Byx,-13);Bzx1=roundn(Bzx,-13);Bzy1=roundn(Bzy,-13);Bxx1=roundn(Bxx,-13);
% B(:,1)=roundn(B(:,1),-13);
% G=[Bxx1 Bxy1 Bxz1;Byx1 Byy1 Byz1;Bzx1 Bzy1 Bzz1];

% 0.01pT                 y = ROUNDN(x,n) 
% Bxy1=roundn(Bxy,-14);Bxz1=roundn(Bxz,-14);Byy1=roundn(Byy,-14);Byz1=roundn(Byz,-14);Bzz1=roundn(Bzz,-14);   
% Byx1=roundn(Byx,-14);Bzx1=roundn(Bzx,-14);Bzy1=roundn(Bzy,-14);Bxx1=roundn(Bxx,-14);
% B(:,1)=roundn(B(:,1),-14);
% G=[Bxx1 Bxy1 Bxz1;Byx1 Byy1 Byz1;Bzx1 Bzy1 Bzz1];

% 0.001pT                 y = ROUNDN(x,n) 
Bxy1=roundn(Bxy,-15);Bxz1=roundn(Bxz,-15);Byy1=roundn(Byy,-15);Byz1=roundn(Byz,-15);Bzz1=roundn(Bzz,-15);   
Byx1=roundn(Byx,-15);Bzx1=roundn(Bzx,-15);Bzy1=roundn(Bzy,-15);Bxx1=roundn(Bxx,-15);
B(:,1)=roundn(B(:,1),-15);
G=[Bxx1 Bxy1 Bxz1;Byx1 Byy1 Byz1;Bzx1 Bzy1 Bzz1];

%1pT=10^-12T=10^-3nT
% Bxy1=roundn(Bxy,-12);Bxz1=roundn(Bxz,-12);Byy1=roundn(Byy,-12);Byz1=roundn(Byz,-12);Bzz1=roundn(Bzz,-12);   
% Byx1=roundn(Byx,-12);Bzx1=roundn(Bzx,-12);Bzy1=roundn(Bzy,-12);Bxx1=roundn(Bxx,-12);
% B(:,1)=roundn(B(:,1),-12);
% G=[Bxx1 Bxy1 Bxz1;Byx1 Byy1 Byz1;Bzx1 Bzy1 Bzz1];

r=-3*(G\B(:,1)); 
rr=sqrt((r(1)-x1)^2+(r(2)-y1)^2+(r(3)-z1)^2);
rr1=[rr x1 d1 d2];

rn=[rn;rr1];


end    
  

end

end

rn1=rn(1:dd,:);              rn2=rn(dd+1:2*dd,:);         rn3=rn(2*dd+1:3*dd,:);       rn4=rn(3*dd+1:4*dd,:);
rn12=reshape(rn1(:,1),11,11);rn22=reshape(rn2(:,1),11,11);rn32=reshape(rn3(:,1),11,11);rn42=reshape(rn4(:,1),11,11);

figure
contourf(c1,c2,rn12);xlabel('d_1(m)','FontSize',sizexy),ylabel('d_2(m)','FontSize',sizexy),title('(a) x=5m','FontSize',sizexy);

set(gca,'xtick',[0.05 0.055 0.06 0.065 0.07 0.075 0.08 0.085 0.09 0.095 0.1]);
set(gca,'xtickLabel',{'0.05',' ','0.06',' ','0.07',' ','0.08',' ','0.09',' ','0.10'});set(gca,'LineWidth',1.5)
set(gca,'ytick',[0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55]);
set(gca,'ytickLabel',{'0.05',' ','0.15',' ','0.25',' ','0.35',' ','0.45',' ','0.55'});set(gca,'LineWidth',1.5)
set(gca,'FontName','Times New Roman','FontSize',sizexy)

axis square
axis([0.05 0.1 0.05 0.55])
set (gcf,'units','centimeters','Position',[35,12,13.5,13]);
colorbar('position',[0.85 0.215 0.04 0.65],'FontSize',sizexy);

set(gca,'Position',[0.185 0.215 0.65 0.65]);
% 
figure
contourf(c1,c2,rn22);xlabel('d_1(m)','FontSize',sizexy),ylabel('d_2(m)','FontSize',sizexy),title('(b) x=10m','FontSize',sizexy);

set(gca,'xtick',[0.05 0.055 0.06 0.065 0.07 0.075 0.08 0.085 0.09 0.095 0.1]);
set(gca,'xtickLabel',{'0.05',' ','0.06',' ','0.07',' ','0.08',' ','0.09',' ','0.10'});set(gca,'LineWidth',1.5)
set(gca,'ytick',[0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55]);
set(gca,'ytickLabel',{'0.05',' ','0.15',' ','0.25',' ','0.35',' ','0.45',' ','0.55'});set(gca,'LineWidth',1.5)
set(gca,'FontName','Times New Roman','FontSize',sizexy)

axis square
axis([0.05 0.1 0.05 0.55])
set (gcf,'units','centimeters','Position',[25,12,13.5,13]);
colorbar('position',[0.85 0.215 0.04 0.65],'FontSize',sizexy);

set(gca,'Position',[0.185 0.215 0.65 0.65]);

figure
contourf(c1,c2,rn32);xlabel('d_1(m)','FontSize',sizexy),ylabel('d_2(m)','FontSize',sizexy),title('(c) x=15m','FontSize',sizexy);

set(gca,'xtick',[0.05 0.055 0.06 0.065 0.07 0.075 0.08 0.085 0.09 0.095 0.1]);
set(gca,'xtickLabel',{'0.05',' ','0.06',' ','0.07',' ','0.08',' ','0.09',' ','0.10'});set(gca,'LineWidth',1.5)
set(gca,'ytick',[0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55]);
set(gca,'ytickLabel',{'0.05',' ','0.15',' ','0.25',' ','0.35',' ','0.45',' ','0.55'});set(gca,'LineWidth',1.5)
set(gca,'FontName','Times New Roman','FontSize',sizexy)

axis square
axis([0.05 0.1 0.05 0.55])
set (gcf,'units','centimeters','Position',[18,12,13.5,13]);
colorbar('position',[0.85 0.215 0.04 0.65],'FontSize',sizexy);

set(gca,'Position',[0.185 0.215 0.65 0.65]);

figure
contourf(c1,c2,rn42);xlabel('d_1(m)','FontSize',sizexy),ylabel('d_2(m)','FontSize',sizexy),title('(d) x=20m','FontSize',sizexy);

set(gca,'xtick',[0.05 0.055 0.06 0.065 0.07 0.075 0.08 0.085 0.09 0.095 0.1]);
set(gca,'xtickLabel',{'0.05',' ','0.06',' ','0.07',' ','0.08',' ','0.09',' ','0.10'});set(gca,'LineWidth',1.5)
set(gca,'ytick',[0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55]);
set(gca,'ytickLabel',{'0.05',' ','0.15',' ','0.25',' ','0.35',' ','0.45',' ','0.55'});set(gca,'LineWidth',1.5)
set(gca,'FontName','Times New Roman','FontSize',sizexy)

axis square
axis([0.05 0.1 0.05 0.55])
set (gcf,'units','centimeters','Position',[10,12,13.5,13]);
colorbar('position',[0.85 0.215 0.04 0.65],'FontSize',sizexy);

set(gca,'Position',[0.185 0.215 0.65 0.65]);

% xlswrite('D:\¡ı÷“œÈ\matlab data\baseline_rn',rn)
% xlswrite('D:\¡ı÷“œÈ\matlab data\baseline_x5',rn12)
% xlswrite('D:\¡ı÷“œÈ\matlab data\baseline_x10',rn22)
% xlswrite('D:\¡ı÷“œÈ\matlab data\baseline_x15',rn32)
% xlswrite('D:\¡ı÷“œÈ\matlab data\baseline_x20',rn42)


