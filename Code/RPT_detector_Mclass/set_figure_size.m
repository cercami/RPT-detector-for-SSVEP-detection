function set_figure_size( x, y )

    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) x y]);

end

