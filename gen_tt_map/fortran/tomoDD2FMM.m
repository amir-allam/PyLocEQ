%Convert tomoDD velocity model for use by FMM code
%  This guy is a nice and easy ASCII format:
%    ngrids ntypes
%    dradius dlat dlon     !grid spacings
%    oradius olat olon     !origin of the grid
%    V(1,1,1)
%    V(2,1,1)
%    V(1,2,1)
%    V(1,1,2)
%    etc...

% load combined_vel_models %Velocity model smoothed CVM and ABZ2012
%% Convert UTM to lat/lon
clear lon
clear lat
for iy = 1:length(C.qy)
    for ix = 1:length(C.qx)
        [lat(iy,ix) lon(iy,ix)]=utm2deg(C.qx(ix),C.qy(iy),'11 S');
    end
end
%% Find lon/lat grid spacings
dx = distance(lat(1,1),lon(1,1),lat(1,2),lon(1,2));
dy = distance(lat(1,1),lon(1,1),lat(2,1),lon(2,1));
dz = C.newqz(2)-C.newqz(1);
%%Find origin: SW corner, shallowest depth !!!I THINK
clat=min(min(lat));
clon=min(min(lon));
cz = min(C.newqz);
%Convert depth to radius
earthrad=6376.0;
oradius=earthrad-max(C.newqz);
%% Make a vector out of all velocities
ic=0;
clear output
for iz = 1:length(C.newqz)
    for iy = 1:C.ny
        for ix =1:C.nx
            ic=ic+1;
            output(ic)=C.smoothVp(iy,ix,iz);
        end
    end
end
%% Write to an ASCII file
outfnam='poop.txt';
% The first line contains the number of grids (1) and the number of velocities (1 for just Vp)
dlmwrite(outfnam,[1 1],'delimiter',' ','precision','%u');
dlmwrite(outfnam,[length(C.newqz) C.ny C.nx ],'delimiter',' ','precision','%u','-append');
dlmwrite(outfnam,[dz dy dx],'delimiter',' ','precision','%8.5f','-append');
dlmwrite(outfnam,[cz clat clon],'delimiter',' ','precision','%8.5f','-append');
dlmwrite(outfnam,output','delimiter',' ','precision','%5.3f','-append');

%% Now make the interfaces file
% This is just an adhoc version with 2 interfaces top and bottom of
% constant value
outfnam='tmp_interfaces';
r1=6371;
r2=r1-max(C.newqz); %Just a little lower than the actual minimum
nlines_out=C.nx*C.ny; %The number of output lines per interfaces
ninterfaces=2;
dlmwrite(outfnam,ninterfaces,'delimiter',' ','precision','%u');
dlmwrite(outfnam,[C.nx C.ny],'delimiter',' ','precision','%u','-append');
dlmwrite(outfnam,[dy dx],'delimiter',' ','precision','%8.5f','-append');
dlmwrite(outfnam,[clat clon],'delimiter',' ','precision','%8.5f','-append');
% dlmwrite(outfnam,zeros(nlines_out,1)+r1,'delimiter',' ','precision','%5.2f','-append');
%dlmwrite(outfnam,zeros(nlines_out,1)+r2,'delimiter',' ','precision','%5.2f','-append');
%% Read topo file and output topography as the top interface
fn='/Users/aallam/gmt/100m_poop.grd'
x=ncread(fn,'x'); %lon
y=ncread(fn,'y'); %lat
z=ncread(fn,'z'); %elevation
z=z';
if 0 %plot
    figure(4);clf;hold on;
    pcolor(x,y,z); shading flat; axis equal
    plot(faultlon,faultlat,'k');
    xlim([x(1) x(end)])
    ylim([y(1) y(end)])
end
%find the topo at the same points as the output velocity
clear output
ic=0;
for iy = 1:C.ny
    fprintf('%u\n',iy)
    for ix =1:C.nx
        ic=ic+1;
        near_ix=findnearest(lon(iy,ix),x);
        near_iy=findnearest(lat(iy,ix),y);
        topo_out(iy,ix)=z(near_iy,near_ix);
        output(ic)=topo_out(iy,ix);
    end
end
if 0
    figure(5);clf;hold on;
    pcolor(lon,lat,topo_out); shading flat; axis equal
    plot(faultlon,faultlat,'k');
    xlim([lon(1) lon(end)])
    ylim([lat(1) lat(end)])
end
%% Write the topography to the interfaces file
dlmwrite(outfnam,output','delimiter',' ','precision','%5.3f','-append');
dlmwrite(outfnam,zeros(nlines_out,1)+r2,'delimiter',' ','precision','%5.2f','-append');
%% Testing topo
% tst=ginput(1)
% tmpx=findnearest(tst(1),x);
% tmpy=findnearest(tst(2),y);
% z(tmpy,tmpx)
