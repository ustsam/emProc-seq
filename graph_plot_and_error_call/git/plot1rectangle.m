function plot1rectangle(x1,x2,y1,y2,cc)

x = [x1, x2, x2, x1];
y = [y1, y1, y2, y2];
patch(x,y,cc);