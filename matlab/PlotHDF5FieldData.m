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
    % Cartesian mesh
    [X Y Z] = meshgrid(mesh.lines{1},mesh.lines{2},mesh.lines{3});
    for n=1:numel(fields.TD.values)
        % since Matlab 7.1SP3 the field needs to be reordered
        fields.TD.values{n} = permute(fields.TD.values{n},[2 1 3 4]); % reorder: y,x,z (or y,x)
    end
else
    disp(['PlotHDF5FieldData:: Error: unknown mesh type ' num2str(mesh.type)]);
end

max_amp = 0;

if (component>0)
    for n=1:numel(fields.TD.values)
        Field{n} = fields.TD.values{n}(:,:,:,component);
    end
else 
    for n=1:numel(fields.TD.values)
        fx = fields.TD.values{n}(:,:,:,1);
        fy = fields.TD.values{n}(:,:,:,2);
        fz = fields.TD.values{n}(:,:,:,3);
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
    if size(Field{n},3) > 1
        % Field is a volume
        hsurfaces = slice(X,Y,Z, Field{n} , PlotArgs.slice{:});
        set(hsurfaces,'FaceColor','interp','EdgeColor','none');
    else
        % Field is already a 2D cut
        pcolor(X,Y,Field{n});
        shading( 'interp' );
        xlabel( 'x' );
        ylabel( 'y' );
    end
    title(fields.TD.names{n});
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
