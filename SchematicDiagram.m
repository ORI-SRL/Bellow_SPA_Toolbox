function out = SchematicDiagram(EffectiveLength,AverageRadius,BellowNum,InnerRadius,WallThickness)


OuterRadius = AverageRadius * 2 - InnerRadius;
UpperArcRadius = EffectiveLength/BellowNum/2/2;
LowerArcRadius = EffectiveLength/BellowNum/2/2;
flank = OuterRadius - InnerRadius - UpperArcRadius - LowerArcRadius;
ContactWidth = AverageRadius * 2;

r1 = UpperArcRadius;
r2 = LowerArcRadius;
Phi1 = pi/2;
Phi2 = pi/2; 
f = flank; 
Beta = 0; 
ro = OuterRadius;
ri = InnerRadius;
rm = AverageRadius;
t = WallThickness;
L0 = EffectiveLength;
D = ContactWidth;

% frame connector 
fx_1 = (ro+t-ri-2*t)*[1 1]; % left vertical line
fy_1 = [0 6]; 
fx_2 = [(ro+t-ri-2*t) (t+ro+ro+2*t)]; % upper horizontal line
fy_2 = [6 6];
fx_3 = (t+ro+ro+2*t)*[1 1]; % right vertical line
fy_3 = [6 0];
fx_4 = [(t+ro+ro+2*t) (ro+t-ri-2*t)]; % lower horizontal line
fy_4 = [0 0];

SchematicDiagram = figure('visible','off');
plot(fx_1, fy_1,'-b');
hold on
plot(fx_2, fy_2,'-b');
plot(fx_3, fy_3,'-b');
plot(fx_4, fy_4,'-b');

fX = [fx_1 fx_2 fx_3 fx_4];
fY = [fy_1 fy_2 fy_3 fy_4];
fill(fX, fY, [0 0.4470 0.7410],'LineStyle','none');


% Bellow Layer
n = 1;
bx_1 = (ro+t-ri-t)*[1 1]; % left vertical line to the frame connector
by_1 = [0 -t];
bx_2 = (ro+t-ri)*[1 1]; % middle vertical line to the frame connector
by_2 = [-t 0];
bX = [bx_1 bx_2];
bY = [by_1 by_2];

plot(bx_1,by_1, '-k');
plot(bx_2,by_2, '-k');

fill(bX, bY, [0.5 0.5 0.5],'LineStyle','none');

b1X = [];
b1Y = [];
b2X = [];
b2Y = [];

for n = 1:1:BellowNum
    %Inner Bellow Layer
    LeftLowerArcXCenter = t + ro - ri - r2;
    LeftLowerArcYCenter = -t-(n-1)*2*(r1+r2);
    UpperArcXCenter = t + r1 ;
    UpperArcYCenter = -t-(n-1)*2*(r1+r2)-(r1+r2);
    RightLowerArcXCenter = t + ro - ri - r2;
    RightLowerArcYCenter = -t-n*2*(r1+r2);
    ThetaLeftLower = pi*3/2+Phi2: -0.001: pi*3/2;
    ThetaUpper = pi/2: 0.01 : pi/2+Phi1*2;
    ThetaRightLower = Phi2: -0.001: 0;
    bx_1 = LeftLowerArcXCenter + r2*cos(ThetaLeftLower);
    by_1 = LeftLowerArcYCenter + r2*sin(ThetaLeftLower);
    bx_3 = UpperArcXCenter + r1*cos(ThetaUpper);
    by_3 = UpperArcYCenter + r1*sin(ThetaUpper);
    bx_2 = [bx_1(length(bx_1)) bx_3(1)];
    by_2 = [by_1(length(by_1)) by_3(1)];
    bx_5 = RightLowerArcXCenter + r2*cos(ThetaRightLower);
    by_5 = RightLowerArcYCenter + r2*sin(ThetaRightLower);
    bx_4 = [bx_3(length(bx_3)) bx_5(1)];
    by_4 = [by_3(length(by_3)) by_5(1)];
    
    plot(bx_1,by_1, '-k');
    plot(bx_2,by_2, '-k');
    plot(bx_3,by_3, '-k');
    plot(bx_4,by_4, '-k');
    plot(bx_5,by_5, '-k');
    
    b1X = [b1X bx_1 bx_2 bx_3 bx_4 bx_5];
    b1Y = [b1Y by_1 by_2 by_3 by_4 by_5];
    
    %Outer Bellow Layer
    LeftLowerArcXCenter = t + ro - ri - r2;
    LeftLowerArcYCenter = -t-(n-1)*2*(r1+r2);
    UpperArcXCenter = t + r1 ;
    UpperArcYCenter = -t-(n-1)*2*(r1+r2)-(r1+r2);
    RightLowerArcXCenter = t + ro - ri - r2;
    RightLowerArcYCenter = -t-n*2*(r1+r2);
    ThetaLeftLower = pi*3/2: 0.001: pi*3/2+Phi2;
    ThetaUpper = pi/2+Phi1*2: -0.001 : pi/2;
    ThetaRightLower = 0: 0.001: Phi2;
    bx_1 = LeftLowerArcXCenter + (r2-t)*cos(ThetaLeftLower);
    by_1 = LeftLowerArcYCenter + (r2-t)*sin(ThetaLeftLower);
    bx_3 = UpperArcXCenter + (r1+t)*cos(ThetaUpper);
    by_3 = UpperArcYCenter + (r1+t)*sin(ThetaUpper);
    bx_2 = [bx_3(length(bx_3)) bx_1(1)];
    by_2 = [by_3(length(by_3)) by_1(1)];
    bx_5 = RightLowerArcXCenter + (r2-t)*cos(ThetaRightLower);
    by_5 = RightLowerArcYCenter + (r2-t)*sin(ThetaRightLower);
    bx_4 = [bx_5(length(bx_5)) bx_3(1)];
    by_4 = [by_5(length(by_5)) by_3(1)];
    
    plot(bx_1,by_1, '-k');
    plot(bx_2,by_2, '-k');
    plot(bx_3,by_3, '-k');
    plot(bx_4,by_4, '-k');
    plot(bx_5,by_5, '-k');
    
    b2X = [bx_5 bx_4 bx_3 bx_2 bx_1 b2X];
    b2Y = [by_5 by_4 by_3 by_2 by_1 b2Y];
end
fill([b1X b2X], [b1Y b2Y], [0.9 0.85 0],'LineStyle','none');


%Bottom Layer
box_1 = (t+ro+ro+2*t)*[1 1]; % Outer bottom line
boy_1 = [-(L0+2*t) 0];
box_2 = (t+ro+ri)*[1 1];  
boy_2 = [0 -t]; 
box_3 = [(ro+t+ri) (ro+t+ro)]; 
boy_3 = -t*[1 1];
box_4 = (t+ro+ro)*[1 1]; % inner bottom line
boy_4 = [-t -(L0+t)];
box_5 = [(t+ro+ro) (ro+t-ri)]; 
boy_5 = -(L0+t)*[1 1];
box_6 = (ro+t-ri-t)*[1 1]; 
boy_6 = [-(L0+t) -(L0+2*t)];
box_7 = [(ro+t-ri-t) (t+ro+ro+2*t)]; 
boy_7 = -(L0+2*t)*[1 1];

boX = [box_1 box_2 box_3 box_4 box_5 box_6 box_7];
boY = [boy_1 boy_2 boy_3 boy_4 boy_5 boy_6 boy_7];

plot(box_1, boy_1, '-k');
plot(box_2, boy_2, '-k');
plot(box_3, boy_3, '-k');
plot(box_4, boy_4, '-k');
plot(box_5, boy_5, '-k');
plot(box_6, boy_6, '-k');
plot(box_7, boy_7, '-b');

fill(boX,boY,[0.5 0.5 0.5],'LineStyle','none')


% End Cap
cx_1 = [(t+ro+ro+2*t) (ro+t-ri-t)]; % upper horizontal line
cy_1 = -(L0+2*t)*[1 1]; 
cx_2 = (ro+t-ri-t)*[1 1]; % left vertical line
cy_2 = [-(L0+2*t) -(L0+3*t)];
cx_3 = [(ro+t-ri-t) (ro+t+ro+2*t)]; % lower horizontal line
cy_3 = -(L0+3*t)*[1 1];
cx_4 = (ro+t+ro+2*t)*[1 1]; % right vertical line
cy_4 = [-(L0+3*t) -(L0+2*t)];

plot(cx_1, cy_1, '-b');
plot(cx_2, cy_2, '-b');
plot(cx_3, cy_3, '-b');
plot(cx_4, cy_4, '-b');

cX = [cx_1 cx_2 cx_3 cx_4];
cY = [cy_1 cy_2 cy_3 cy_4];

fill(cX,cY,[0 0.4470 0.7410],'LineStyle','none')






axis equal
grid on
box on

saveas(SchematicDiagram, 'SchematicDiagram.jpg')
set(SchematicDiagram, 'visible', 'on'); 
hold off
close(SchematicDiagram)
end