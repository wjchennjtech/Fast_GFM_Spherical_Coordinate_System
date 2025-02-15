function[gravity_potential,gravity_abnomal,gravity_gradient]=synthesis_slice(Cnm,Snm,nmax,Lat,Lon,height,nz)

R=6371*1e+3;
GM=3986005*1e+8;
nlat=size(Lat,2);
nlon=size(Lon,2);

gravity_potential=zeros(nlat,nlon);
gravity_abnomal=zeros(nlat,nlon);
gravity_gradient=zeros(nlat,nlon);
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
    potential=potential(ones(1,nlat),:);
    abnomal=abnomal(ones(1,nlat),:);
    gradient=gradient(ones(1,nlat),:);


    gravity_potential=coef.*potential.*plm(n,m,90-Lat)*CNM(:,m+1)*cosd(m*Lon)...
        +coef.*potential.*plm(n,m,90-Lat)*SNM(:,m+1)*sind(m*Lon)...
        +gravity_potential;
    gravity_abnomal=coef.*abnomal.*plm(n,m,90-Lat)*CNM(:,m+1)*cosd(m*Lon)...
        +coef.*abnomal.*plm(n,m,90-Lat)*SNM(:,m+1)*sind(m*Lon)...
        +gravity_abnomal;
    gravity_gradient=coef.*gradient.*plm(n,m,90-Lat)*CNM(:,m+1)*cosd(m*Lon)...
        +coef.*gradient.*plm(n,m,90-Lat)*SNM(:,m+1)*sind(m*Lon)...
        +gravity_gradient;

end

end



