function [] = plot_shape(shape, X, Y, xlimits, ylimits, E)
%plot_shape Plots a shape structure with the electric field overlayed on
% top.

maxv = max(E(:)) + 0.0001; % Offset required to produce correct color output
minv = min(E(:));

pshape = shape;
[M, N] = size(shape);

for x = 1:N
    for y = 1:M
        if (shape(x, y) == 0)
            pshape(x, y) = minv;
        else
            pshape(x, y) = maxv;
        end
    end
end

figure;
surf(X, Y, E);
hold on
surf(X, Y, pshape);
hold off

title('');
xlabel('Distance (m)');
ylabel('Distance (m)');

xlim(xlimits);
ylim(ylimits);

shading interp;
view(2);

h = colorbar;
ylabel(h, 'Total Electric Field (dB)');

end