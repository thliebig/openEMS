function [CSX mesh] = CreateCRLH(CSX, mesh, CRLH, resolution, translate)
% function [CSX mesh] = CreateCRLH(CSX, mesh, CRLH, resolution, translate)
% 
% support function to create a CRLH unit cell
% 
% currently used by Tutorials/CRLH_Extraction
% 
% Tested with
%  - Matlab 2009b
%  - openEMS v0.0.23
%
% (C) 2011 Thorsten Liebig <thorsten.liebig@gmx.de>

if (nargin<5)
    translate = [0 0 0];
end

CSX = AddMetal(CSX, 'metal_top');
one_two_third = [-resolution/3 2*resolution/3];

start = [-CRLH.LL/2 -CRLH.LW/2 CRLH.TopSig]+translate;
stop  = [-CRLH.GLT/2  CRLH.LW/2 CRLH.TopSig]+translate;
CSX = AddBox(CSX, 'metal_top', 10, start, stop);
mesh.x = [mesh.x start(1) stop(1)+one_two_third];
mesh.y = [mesh.y start(2)-one_two_third stop(2)+one_two_third];

start = [+CRLH.LL/2  -CRLH.LW/2 CRLH.TopSig]+translate;
stop  = [+CRLH.GLT/2  CRLH.LW/2 CRLH.TopSig]+translate;
CSX = AddBox(CSX, 'metal_top', 10, start, stop);
mesh.x = [mesh.x start(1) stop(1)-one_two_third];

CSX = AddMetal(CSX, 'metal_bot');
start = [-(CRLH.LL-CRLH.GLB)/2 -CRLH.LW/2 CRLH.BottomSig]+translate;
stop  = [+(CRLH.LL-CRLH.GLB)/2  CRLH.LW/2 CRLH.BottomSig]+translate;
CSX = AddBox(CSX, 'metal_bot', 10, start, stop);
mesh.x = [mesh.x start(1)-one_two_third stop(1)+one_two_third];

start = [-CRLH.SW/2 -CRLH.LW/2-CRLH.SL CRLH.BottomSig]+translate;
stop  = [+CRLH.SW/2  CRLH.LW/2+CRLH.SL CRLH.BottomSig]+translate;
CSX = AddBox(CSX, 'metal_bot', 10, start, stop);
mesh.x = [mesh.x start(1)-one_two_third stop(1)+one_two_third];
mesh.y = [mesh.y start(2) stop(2)];

CSX = AddMetal(CSX, 'via');
start = [0 -CRLH.LW/2-CRLH.SL+CRLH.SW/2 0]+translate;
stop  = [0 -CRLH.LW/2-CRLH.SL+CRLH.SW/2 CRLH.BottomSig]+translate;
CSX = AddCylinder(CSX, 'via', 10, start, stop, CRLH.VR);
mesh.x = [mesh.x start(1)+[-1 0 1]*CRLH.VR];
mesh.y = [mesh.y start(2)+[-1 0 1]*CRLH.VR];

start(2) = -start(2);
stop(2)  = -stop(2);
CSX = AddCylinder(CSX, 'via', 10, start, stop, CRLH.VR);
mesh.y = [mesh.y start(2)+[-1 0 1]*CRLH.VR];
end