function save_figure_compat(fig,filename,resolution)
%SAVE_FIGURE_COMPAT Save PNG using exportgraphics when available, else print.
if nargin<3, resolution=180; end
if exist('exportgraphics','file')==2
    exportgraphics(fig,filename,'Resolution',resolution);
else
    set(fig,'PaperPositionMode','auto');
    print(fig,filename,'-dpng',sprintf('-r%d',resolution));
end
end
