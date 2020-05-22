%NAMA : MUHAMMAD REFRI ANSYARI
%NIM  : 12916028

%MODUL IV HIDRODINAMIKA 2D
%TUGAS 1 (DANAU)

clear all;
%%  INISIALISASI VARIABEL
%   Domain Model
L  =1000;            % Ukuran Domain Model (meter)
dx =50;              % Ukuran Grid ruang x (meter)
dy =50;              % Ukuran Grid ruang y (meter)
T  =7200;            % Waktu Simulasi (detik)
d  =10;              % Kedalaman Danau (meter)
dt =2;               % Langkah Waktu

%  Parameter Model
g    =9.81;          % Percepatan Gravitasi (m/s^2)
A    =1.5;           % Amplitudo Gelombang (meter)
pg   =12;            % Periode Gelombang (detik)
r    =0.005;         % Gesekan Dasar
lamda=0.001;         % Gesekan pemukaan
Sigma=((2*pi)/pg);   % Frekuensi sudut
Wx   =-1;            % Kecepatan angin arah x konstan (m/s)
Wy   =1;             % Kecepatan angin arah y konstan (m/s)

%Jumlah Iterasi Program
imax=L/dx;
jmax=L/dy;
nmax=T/dt;
nout=300;


%% INPUT
%  Input Domain Utama
danau= load ('danau.txt');

%  Nilai Awal
for j=1:jmax
    for i=1:imax
        Ut(i,j,1)=0;
        Vt(i,j,1)=0;
        elev(i,j,1)=0;
    end
end
Ut=zeros(imax,jmax,nmax);
Vt=zeros(imax,jmax,nmax);
u=zeros(imax,jmax,nmax);
v=zeros(imax,jmax,nmax);

%% PROSES
%  Langkah Iterasi
for n=1:nmax
    for j=2:jmax-1
        for i=2:imax-1
    %Variabel ustar dan vstar
            Vs=(Vt(i,j)+Vt(i+1,j)+Vt(i,j-1)+Vt(i+1,j-1))/4;
            Us=(Ut(i,j)+Ut(i-1,j)+Ut(i-1,j+1)+Ut(i,j+1))/4;
            
    %Variabel HU dan HV
            Hu=(elev(i+1,j,n)+elev(i,j,n)+(2*d))/2;
            Hv=(elev(i,j+1,n)+elev(i,j,n)+(2*d))/2;
            
    %Variabel RX dan RY
            Rx=1/(1+((r*dt*sqrt((Ut(i,j,n))^2+Us^2))/(Hu^2)));
            Ry=1/(1+((r*dt*sqrt((Vt(i,j,n))^2+Vs^2))/(Hv^2)));
            
    %Variabel U dan V
            Ut(i,j,n+1)=(Ut(i,j,n)-((g*dt*Hu/dx)*(elev(i+1,j,n)-elev(i,j,n)))+(dt*lamda*Wx*(sqrt(Wx^2+Wy^2))))*Rx;
            Vt(i,j,n+1)=(Vt(i,j,n)-((g*dt*Hv/dy)*(elev(i,j+1,n)-elev(i,j,n)))+(dt*lamda*Wy*(sqrt(Wx^2+Wy^2))))*Ry;
            
    %Menghitung Elevasi
            elev(i,j,n+1)=elev(i,j,n)-(dt*(((Ut(i,j,n+1)-Ut(i-1,j,n+1))/dx)+((Vt(i,j,n+1)-Vt(i,j-1,n+1))/dy)));
            Hmax=max(d+elev(i,j,n));
            end
    end
    
    %Interpretasi Hasil Pada Domain
    for j=1:jmax
        for i=1:imax
            u(i,j,n+1)=(Ut(i,j,n+1)/Hu)*danau(i,j);
            v(i,j,n+1)=(Vt(i,j,n+1)/Hv)*danau(i,j);
            elevasi(i,j,n+1)=elev(i,j,n+1)*danau(i,j);
            kec(i,j,n+1)=sqrt(u(i,j,n+1)^2+v(i,j,n+1)^2);
        end
    end
    
    % Syarat Kestabilan Model
CFL=dx/sqrt(2*g*Hmax);
if dt <= CFL
    disp('Model Stabil');
else
    disp ('Model Tidak Stabil');
end
    
%% OUTPUT
%  Plot Elevasi Danau Sintetik
figure (1)    
if mod(n,nout)==0
    elevasif=flipud(elevasi(:,:,n));
    pcolor(elevasif);
    h=colorbar; 
    colormap(cool)
    shading faceted;
    xlabel('Bujur (m)');
    ylabel('Lintang (m)');
    set (gca,'XTick',0:5:20); set(gca,'XTickLabel',{' ','250','500','750','1000'});
    set (gca,'YTick',0:5:20); set(gca,'YTickLabel',{' ','250','500','750','1000'});
    grid on; box on;
    ylabel(h,'Elevasi Danau(m)');
    if n==nmax
        t='120 Menit';
    else
        t=[num2str((n/nout)*10) ' Menit'];
    end
    title(['Variasi Nilai Elevasi di Danau Sintetik Saat t=' t])
    filename=['evlake',num2str(n/nout),'.jpg'];
    saveas(gcf,filename)
    clf;
end

%  Plot Magnitudo Kecepatan Arus Danau Sintetik
figure(2)
if mod(n,nout)==0
    kecf=flipud(kec(:,:,n));
    uf=flipud(u(:,:,n));
    vf=flipud(v(:,:,n));
    pcolor(kecf);
    h=colorbar; 
    colormap(winter)
    ylabel(h,'Kecepatan Arus (m/s)');
    shading faceted;
    hold on;
    quiver(uf,vf,'red');
    xlabel('Bujur (m)');
    ylabel('Lintang (m)');
    set (gca,'XTick',0:5:20); set(gca,'XTickLabel',{' ','250','500','750','1000'});
    set (gca,'YTick',0:5:20); set(gca,'YTickLabel',{' ','250','500','750','1000'});
    grid on; box on;
    if n==nmax
        t='120 Menit (2 Jam)';
    else
        t=[num2str((n/nout)*10) ' Menit'];
    end
    title(['Variasi Kecepatan Arus di Danau Sintetik Saat t=' t])
    filename=['vlake',num2str(n/nout),'.jpg'];
    saveas(gcf,filename)
    clf;
end
end


