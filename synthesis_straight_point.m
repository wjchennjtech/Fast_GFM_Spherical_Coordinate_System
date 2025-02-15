function[gravity_potential,gravity_abnomal,gravity_gradient]=synthesis_straight_point(cnm_snm,nmax,Lat,Lon,height)

R=6371*1e+3;
GM=3986005*1e8;
row=size(Lat,1);

gravity_potential=zeros(row,1);
gravity_abnomal=zeros(row,1);
gravity_gradient=zeros(row,1);
rp=height+R;%计算高度

Cnm=zeros(nmax+1,nmax+1);
Snm=zeros(nmax+1,nmax+1);

nmnumber=0;
for i=1:nmax+1
    nmnumber=i+nmnumber;
end

for i=1: nmnumber

    Cnm(cnm_snm(i,1)+1,cnm_snm(i,2)+1)=cnm_snm(i,3);
    Snm(cnm_snm(i,1)+1,cnm_snm(i,2)+1)=cnm_snm(i,4);

end%系数文件拆分
for m=0:nmax

    n=0:nmax;
    coef=GM/R;

    potential=(R/rp).^(n+1);
    abnomal=(n+1).*potential./rp;
    gradient=(n+1).*(n+2).*potential./(rp*rp);
    potential=potential(ones(1,row),:);
    abnomal=abnomal(ones(1,row),:);
    gradient=gradient(ones(1,row),:);

    gravity_potential=coef.*potential.*plm(n,m,90-Lat)*Cnm(:,m+1).*cosd(m*Lon)...
        +coef.*potential.*plm(n,m,90-Lat)*Snm(:,m+1).*sind(m*Lon)...
        +gravity_potential;
    gravity_abnomal=coef.*abnomal.*plm(n,m,90-Lat)*Cnm(:,m+1).*cosd(m*Lon)...
        +coef.*abnomal.*plm(n,m,90-Lat)*Snm(:,m+1).*sind(m*Lon)...
        +gravity_abnomal;
    gravity_gradient=coef.*gradient.*plm(n,m,90-Lat)*Cnm(:,m+1).*cosd(m*Lon)...
        +coef.*gradient.*plm(n,m,90-Lat)*Snm(:,m+1).*sind(m*Lon)...
        +gravity_gradient;

end

end

