% https://www.mathworks.com/help/matlab/matlab_external/passing-data-to-python.html
import reesaurora.rees_model.loadaltenergrid 

minalt=30; 
nalt=286;
isotropic=false;
glat=65; 
glon=-148;

t = '2013-03-31T12:00:00Z'

disp('not working, perhaps a limitation of Matlab right now')
[z,E] = loadaltenergrid(minalt,nalt);
