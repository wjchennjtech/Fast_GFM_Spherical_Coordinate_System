function   [Cnm,Snm]=SHA(data,nmax,dlat,dlon,nlat,nlon)
nmax1=nmax+1;
mmax=nmax;
Cnm=zeros(nmax1,nmax1);
Snm=zeros(nmax1,nmax1);

latmax=max(data(:,2));
latmin=min(data(:,2));
lonmax=max(data(:,1));
lonmin=min(data(:,1));

Lon = linspace(lonmin,lonmax,nlon);
Lat = linspace(latmax,latmin,nlat);

data_t=reshape(data(:,3),nlon,nlat);
data_t=data_t';

d=cosd(Lat);
cosLAT=d(ones(1,nmax+1),:);

dlat=dlat*pi/180;
dlon=dlon*pi/180;

for m=0:mmax

         n=0:nmax;

    Cnm(:,m+1)=1/(4*pi)*((dlat*dlon).*(plm(n,m,90-Lat)'.*cosLAT)*(data_t*cosd(m*Lon)'));
    Snm(:,m+1)=1/(4*pi)*((dlat*dlon).*(plm(n,m,90-Lat)'.*cosLAT)*(data_t*sind(m*Lon)'));


end

end