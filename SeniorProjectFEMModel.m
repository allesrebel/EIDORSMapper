%%  Senior Project FEM Model
%   Based off Actual Electrode Positions in the Model

z_contact= 0.01;
n_elec= 8;
nodes_per_elec= 5;
elec_width= 0.1;
elec_spacing= .25;

xllim=-12.72; xrlim= 12; ydepth=-12.72;
[x,y] = meshgrid( linspace(xllim,xrlim,49), linspace(ydepth,0,31) );
vtx= [x(:),y(:)];

% Refine points close to electrodes - don't worry if points overlap
[x,y] = meshgrid( -5:.25:5, -7:.25:-3 );
vtx= [vtx; x(:),y(:)];

%   First Strip of Electrodes
xgrid=  linspace(-elec_width/2, +elec_width/2, nodes_per_elec)';
x2grid= elec_width* [-5,-4,-3,3,4,5]'/4;
for i=1:n_elec
% Electrode center
  x0= 0;
  y0= (i-1-(n_elec-1)/2)*elec_spacing - 5.5;
  elec_nodes{i}= [x0+ xgrid, y0+0*xgrid];
  vtx= [ vtx; ...
        [x0 + x2grid   ,  y0             + 0*x2grid];
        [x0 + xgrid*1.5,  y0-elec_width/2+ 0*xgrid];
        [x0 + x2grid*1.5, y0-elec_width/2+ 0*x2grid];
        [x0 + xgrid*2   , y0-elec_width  + 0*xgrid];
        [x0 + xgrid*2   , y0-elec_width*2+ 0*xgrid]];
end

%   Second Strip of Electrodes
xgrid=  linspace(-elec_width/2, +elec_width/2, nodes_per_elec)';
x2grid= elec_width* [-5,-4,-3,3,4,5]'/4;
for i=1:n_elec
% Electrode center
  x0= -4.2;
  y0= (i-1-(n_elec-1)/2)*elec_spacing - 4.5;
  elec_nodes{i+n_elec}= [x0+ xgrid, y0+0*xgrid];
  vtx= [ vtx; ...
        [x0 + x2grid   ,  y0             + 0*x2grid];
        [x0 + xgrid*1.5,  y0-elec_width/2+ 0*xgrid];
        [x0 + x2grid*1.5, y0-elec_width/2+ 0*x2grid];
        [x0 + xgrid*2   , y0-elec_width  + 0*xgrid];
        [x0 + xgrid*2   , y0-elec_width*2+ 0*xgrid]];
end

%   Thrid Strip of Electrodes
xgrid=  linspace(-elec_width/2, +elec_width/2, nodes_per_elec)';
x2grid= elec_width* [-5,-4,-3,3,4,5]'/4;
for i=1:n_elec
% Electrode center
  x0= 4.2;
  y0= (i-1-(n_elec-1)/2)*elec_spacing - 4.5;
  elec_nodes{i+n_elec*2}= [x0+ xgrid, y0+0*xgrid];
  vtx= [ vtx; ...
        [x0 + x2grid   ,  y0             + 0*x2grid];
        [x0 + xgrid*1.5,  y0-elec_width/2+ 0*xgrid];
        [x0 + x2grid*1.5, y0-elec_width/2+ 0*x2grid];
        [x0 + xgrid*2   , y0-elec_width  + 0*xgrid];
        [x0 + xgrid*2   , y0-elec_width*2+ 0*xgrid]];
end

fmdl= mk_fmdl_from_nodes( vtx, elec_nodes, z_contact, 'sq_m1');

show_fem(fmdl); axis image

print_convert FEM_ElectodeStrips.png '-density 175'
