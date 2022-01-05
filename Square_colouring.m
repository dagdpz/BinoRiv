function Square_coloring(PX, colour, PY, Bace )
hold on;

xlimit = get(gca, 'XLim');                                                  
ylimit = get(gca, 'YLim');                                                  
colourdflt = [1.0 1.0 0.9];                                                  
if nargin<4, Bace = ylimit(1); end                                          
if nargin<3, PY = [ylimit(2), ylimit(2)]; end                               
if nargin<2, colour = colourdflt; end                                          
if nargin<1, PX = [xlimit(1), xlimit(2)]; end                                

Area_handle = area(PX , PY, Bace);                                         

%hold off;                                                                    
set(Area_handle,'FaceColor', colour);                                        
set(Area_handle,'LineStyle','none');                                        
set(Area_handle,'ShowBaseline','off');                                      
set(gca,'layer','top');                                                     
set(Area_handle.Annotation.LegendInformation, 'IconDisplayStyle','off');    
children_handle = get(gca, 'Children');                                     
set(gca, 'Children', circshift(children_handle,[-1 0]));                    

end