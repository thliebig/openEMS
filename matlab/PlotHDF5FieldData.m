function PlotHDF5FieldData(file, PlotArgs)

plane = PlotArgs.plane;
component = PlotArgs.component;

if (isfield(PlotArgs,'pauseTime'))
    pauseT = PlotArgs.pauseTime;
else
    pauseT = 0.01;
end

mesh = ReadHDF5Mesh(file);
fields = ReadHDF5FieldData(file);

max_amp = 0;

if (strcmp(plane,'zx'))
    [X1 X2] = meshgrid(mesh.lines{3}, mesh.lines{1});
    
    if (component>0)
        for n=1:numel(fields.values)
            Z{n} = squeeze(double(fields.values{n}(:,1,:,component)));
        end
    else 
        for n=1:numel(fields.values)
            fx = squeeze(double(fields.values{n}(:,1,:,1)));
            fy = squeeze(double(fields.values{n}(:,1,:,2)));
            fz = squeeze(double(fields.values{n}(:,1,:,3)));
            Z{n} = sqrt(fx.^2 + fy.^2 + fz.^2);
        end        
    end
end



for n=1:numel(Z)
    amp = max(max(abs(Z{n})));
    if (amp>max_amp)
        max_amp = amp;
    end
end
if (max_amp==0)
    disp('max found amplitude was 0 --> nothing to plot');
    return
end

for n=1:numel(Z)
    surf(X1,X2,Z{n})
    title(fields.names{n});
    
    if (isfield(PlotArgs,'zlim'))
        if ~ischar(PlotArgs.zlim)
            zlim(PlotArgs.zlim);
        elseif strcmp(PlotArgs.zlim,'auto')
            zlim([-max_amp*(component>0) max_amp]);
        end
    end
      
    pause(pauseT)
end