function[gravity_potential,gravity_abnomal,gravity_gradient]=synthesis_slice_point(Cnm,Snm,nmax,Lat,Lon,height,nz)

R=6371*1e+3;
GM=3986005*1e+8;

row=size(Lat,1);

gravity_potential=zeros(row,1);
gravity_abnomal=zeros(row,1);
gravity_gradient=zeros(row,1);
rp=height+R;%观测高度

CNM=zeros(nmax+1,nmax+1);
SNM=zeros(nmax+1,nmax+1);

for i=1:nz
CNM=CNM+Cnm(:,:,i);
SNM=SNM+Snm(:,:,i);
end

for m=0:nmax

    n=0:nmax;
    coef=GM/R;

    potential=(R/rp).^(n+1);
    abnomal=(n+1).*potential./rp;
    gradient=(n+1).*(n+2).*potential./(rp*rp);
    potential=potential(ones(1,row),:);
    abnomal=abnomal(ones(1,row),:);
    gradient=gradient(ones(1,row),:);

    gravity_potential=coef.*potential.*plm(n,m,90-Lat)*CNM(:,m+1).*cosd(m*Lon)...
        +coef.*potential.*plm(n,m,90-Lat)*SNM(:,m+1).*sind(m*Lon)...
        +gravity_potential;
    gravity_abnomal=coef.*abnomal.*plm(n,m,90-Lat)*CNM(:,m+1).*cosd(m*Lon)...
        +coef.*abnomal.*plm(n,m,90-Lat)*SNM(:,m+1).*sind(m*Lon)...
        +gravity_abnomal;
    gravity_gradient=coef.*gradient.*plm(n,m,90-Lat)*CNM(:,m+1).*cosd(m*Lon)...
        +coef.*gradient.*plm(n,m,90-Lat)*SNM(:,m+1).*sind(m*Lon)...
        +gravity_gradient;

end

end



