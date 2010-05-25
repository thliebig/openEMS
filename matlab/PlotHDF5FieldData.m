function PlotHDF5FieldData(file, PlotArgs)
% function PlotHDF5FieldData(file, PlotArgs)
%
% e.g.
% PlotArgs.slice = {0 [10 20] 0};
% PlotArgs.pauseTime=0.01;
% PlotArgs.component=2;
% PlotArgs.Limit = 'auto';
% 
% PlotHDF5FieldData('tmp/Et.h5',PlotArgs)
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

component = PlotArgs.component;

if (isfield(PlotArgs,'pauseTime'))
    pauseT = PlotArgs.pauseTime;
else
    pauseT = 0.01;
end

mesh = ReadHDF5Mesh(file);
fields = ReadHDF5FieldData(file);

if (mesh.type==0)
    % matlab is building a meshgrid incorrect, adressing it like X(y,x,z)
    % ??it's not a bug, it's a feature??
    % check if this workaround is correct  (TL)
    [X Y Z] = meshgrid(mesh.lines{2},mesh.lines{1},mesh.lines{3});
else
    disp(['PlotHDF5FieldData:: Error: unknown mesh type ' num2str(mesh.type)]);
end

max_amp = 0;

if (component>0)
    for n=1:numel(fields.values)
        Field{n} = fields.values{n}(:,:,:,component);
    end
else 
    for n=1:numel(fields.values)
        fx = fields.values{n}(:,:,:,1);
        fy = fields.values{n}(:,:,:,2);
        fz = fields.values{n}(:,:,:,3);
        Field{n} = sqrt(fx.^2 + fy.^2 + fz.^2);
    end        
end

for n=1:numel(Field)
    amp = max(max(max(abs(Field{n}))));
    if (amp>max_amp)
        max_amp = amp;
    end
end

if (max_amp==0)
    disp('max found amplitude was 0 --> nothing to plot');
    return
end

for n=1:numel(Field)
    hsurfaces = slice(X,Y,Z, Field{n} , PlotArgs.slice{:});
    set(hsurfaces,'FaceColor','interp','EdgeColor','none');
    title(fields.names{n});
    %view(3)
    axis equal
    if (isfield(PlotArgs,'Limit'))
        if ~ischar(PlotArgs.Limit)
            caxis(PlotArgs.Limit);
        elseif strcmp(PlotArgs.Limit,'auto')    
            if (component>0)
                caxis([-max_amp,max_amp]);
            else
                caxis([0,max_amp]);
            end
        end
    end
    
    drawnow
    pause(pauseT)
end