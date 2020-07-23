function plot_peakinessSubset(cells,handle)

handle;
category = 'peakiness';


plot_lineWithSEM(cells.(category),[])
title('Peakiness');

end