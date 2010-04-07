function PlotHDF5FieldData(file, PlotArgs)

component = PlotArgs.component;

if (isfield(PlotArgs,'pauseTime'))
    pauseT = PlotArgs.pauseTime;
else
    pauseT = 0.01;
end

mesh = ReadHDF5Mesh(file);
fields = ReadHDF5FieldData(file);

[X Y Z] = meshgrid(double(mesh.lines{1}),double(mesh.lines{2}),double(mesh.lines{3}));

max_amp = 0;

if (component>0)
    for n=1:numel(fields.values)
        Field{n} = double(fields.values{n}(:,:,:,component));
    end
else 
    for n=1:numel(fields.values)
        fx = double(fields.values{n}(:,:,:,1));
        fy = double(fields.values{n}(:,:,:,2));
        fz = double(fields.values{n}(:,:,:,3));
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
%     if (isfield(PlotArgs,'zlim'))
%         if ~ischar(PlotArgs.zlim)
%             zlim(PlotArgs.zlim);
%         elseif strcmp(PlotArgs.zlim,'auto')
%             zlim([-max_amp*(component>0) max_amp]);
%         end
%     end
%       
    drawnow
    pause(pauseT)
end