function out = DeformedSchematicDiagram(EffectiveLength,AverageRadius,BellowNum, InnerRadius, WallThickness, ri_new, rm_new, ro_new, r1_new, r2_new, Phi1_new, Phi2_new, AngularDeflection)




OuterRadius = AverageRadius * 2 - InnerRadius;
UpperArcRadius = EffectiveLength/BellowNum/2/2;
LowerArcRadius = EffectiveLength/BellowNum/2/2;
flank = OuterRadius - InnerRadius - UpperArcRadius - LowerArcRadius;
ContactWidth = AverageRadius * 2;


r1 = UpperArcRadius;
r2 = LowerArcRadius;
Phi1 = 0.5*pi;
Phi2 = 0.5*pi; 
Beta = 0;
f = flank; 
Beta = 0; 
ro = OuterRadius;
ri = InnerRadius;
rm = AverageRadius;
t = WallThickness;
L0 = EffectiveLength;
D = ContactWidth;



% Original Plot
% frame connector 
fx_1 = (ro+t-ri-2*t)*[1 1]; % left vertical line
fy_1 = [0 6]; 
fx_2 = [(ro+t-ri-2*t) (t+ro+ro+2*t)]; % upper horizontal line
fy_2 = [6 6];
fx_3 = (t+ro+ro+2*t)*[1 1]; % right vertical line
fy_3 = [6 0];
fx_4 = [(t+ro+ro+2*t) (ro+t-ri-2*t)]; % lower horizontal line
fy_4 = [0 0];

DeformedSchematicDiagram = figure('visible','off');
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
boy_1 = [-(EffectiveLength+2*t) 0];
box_2 = (t+ro+ri)*[1 1];  
boy_2 = [0 -t]; 
box_3 = [(ro+t+ri) (ro+t+ro)]; 
boy_3 = -t*[1 1];
box_4 = (t+ro+ro)*[1 1]; % inner bottom line
boy_4 = [-t -(EffectiveLength+t)];
box_5 = [(t+ro+ro) (ro+t-ri)]; 
boy_5 = -(EffectiveLength+t)*[1 1];
box_6 = (ro+t-ri-t)*[1 1]; 
boy_6 = [-(EffectiveLength+t) -(EffectiveLength+2*t)];
box_7 = [(ro+t-ri-t) (t+ro+ro+2*t)]; 
boy_7 = -(EffectiveLength+2*t)*[1 1];

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
cy_1 = -(EffectiveLength+2*t)*[1 1]; 
cx_2 = (ro+t-ri-t)*[1 1]; % left vertical line
cy_2 = [-(EffectiveLength+2*t) -(EffectiveLength+3*t)];
cx_3 = [(ro+t-ri-t) (ro+t+ro+2*t)]; % lower horizontal line
cy_3 = -(EffectiveLength+3*t)*[1 1];
cx_4 = (ro+t+ro+2*t)*[1 1]; % right vertical line
cy_4 = [-(EffectiveLength+3*t) -(EffectiveLength+2*t)];

plot(cx_1, cy_1, '-b');
plot(cx_2, cy_2, '-b');
plot(cx_3, cy_3, '-b');
plot(cx_4, cy_4, '-b');

cX = [cx_1 cx_2 cx_3 cx_4];
cY = [cy_1 cy_2 cy_3 cy_4];

fill(cX,cY,[0 0.4470 0.7410],'LineStyle','none')


%Deformed Plot
r1 = r1_new;
r2 = r2_new;
Phi1 = Phi1_new;
Phi2 = Phi2_new; 
Beta = 0;
f = flank; 
Beta = 0; 
ro = ro_new;
ri = ri_new;
rm = rm_new;
AngularDeflection = AngularDeflection;

%bottom layer
DeformedCenterX = (WallThickness+OuterRadius*2)+(EffectiveLength)/AngularDeflection;
DeformedCenterY = -t;

theta_i = pi+AngularDeflection:-0.001:pi;
theta_o = pi:0.001:pi+AngularDeflection;

box_1 = DeformedCenterX + (EffectiveLength)/AngularDeflection*cos(theta_i); % Inner bottom layer
boy_1 = DeformedCenterY + (EffectiveLength)/AngularDeflection*sin(theta_i);
box_2 = DeformedCenterX + ((EffectiveLength)/AngularDeflection-2*t)*cos(theta_o); % outer bottom layer
boy_2 = DeformedCenterY + ((EffectiveLength)/AngularDeflection-2*t)*sin(theta_o);
box_3 = (t+ro+ri)*[1 1];
boy_3 = [0 -t];
box_4 = [(ro+t+ri) (ro+t+ro)]; 
boy_4 = -t*[1 1];
box_5 = [box_1(1) (box_1(1) - (ri+OuterRadius)*cos(AngularDeflection))]; % from inner bottom layer to bellow layer
boy_5 = [boy_1(1) (boy_1(1) - (ri+OuterRadius)*sin(AngularDeflection))];
box_6 = [box_2(length(box_2)) (box_2(length(box_2)) + t*sin(AngularDeflection))]; % from outer bottom layer to the end cap
boy_6 = [boy_2(length(boy_2)) (boy_2(length(boy_2)) - t*cos(AngularDeflection))];
box_7 = [box_6(2) (box_6(2) - (ri+OuterRadius+3*t)*cos(AngularDeflection))]; % along the end cap
boy_7 = [boy_6(2) (boy_6(2) - (ri+OuterRadius+3*t)*sin(AngularDeflection))];
box_8 = fx_2(2)*[1 1]; % frame connector to the outer bottom layer
boy_8 = [0 -t];

boX = [box_8 box_2 box_6 box_7 box_5(2) box_5(1) box_1 box_4(2) box_4(1) box_3(2) box_3(1)];
boY = [boy_8 boy_2 boy_6 boy_7 boy_5(2) boy_5(1) boy_1 boy_4(2) boy_4(1) boy_3(2) boy_3(1)];

plot(box_1, boy_1, '-k');
plot(box_2, boy_2, '-k');
plot(box_3, boy_3, '-k');
plot(box_4, boy_4, '-k');
plot(box_5, boy_5, '-k');
plot(box_6, boy_6, '-k');
plot(box_7, boy_7, '-b');
plot(box_8, boy_8, '-k');

fill(boX,boY,[0.5 0.5 0.5],'LineStyle','none', 'FaceAlpha', 0.5)

% Bellow part
n = 1;
bx_1 = (OuterRadius+WallThickness-InnerRadius-WallThickness)*[1 1];
by_1 = [0 -WallThickness];
bx_2 = (OuterRadius+WallThickness-InnerRadius)*[1 1];
by_2 = [-WallThickness 0 ];

bX = [bx_1 bx_1(2) bx_2(1) bx_2 bx_2(2) bx_1(1)];
bY = [by_1 by_1(2) by_2(1) by_2 by_2(2) by_1(1)];

plot(bx_1,by_1, '-k');
plot(bx_2,by_2, '-k');

fill(bX, bY, [0.5 0.5 0.5],'LineStyle','none', 'FaceAlpha', 0.5);

n = 1;
% Inner Bellow
LeftLowerArcXCenter = DeformedCenterX + ((EffectiveLength)/AngularDeflection  + OuterRadius + ri + r2) * cos(pi + (n-1)*AngularDeflection/BellowNum);
LeftLowerArcYCenter = DeformedCenterY + ((EffectiveLength)/AngularDeflection  + OuterRadius + ri + r2) * sin(pi + (n-1)*AngularDeflection/BellowNum);
UpperArcXCenter = DeformedCenterX + ((EffectiveLength)/AngularDeflection + OuterRadius + ro - r1) * cos(pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2);
UpperArcYCenter = DeformedCenterY + ((EffectiveLength)/AngularDeflection + OuterRadius + ro - r1) * sin(pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2);
RightLowerArcXCenter = DeformedCenterX + ((EffectiveLength)/AngularDeflection + OuterRadius + ri + r2) * cos(pi + (n)*AngularDeflection/BellowNum);
RightLowerArcYCenter = DeformedCenterY + ((EffectiveLength)/AngularDeflection + OuterRadius + ri + r2) * sin(pi + (n)*AngularDeflection/BellowNum);
ThetaLeftLower = (n-1)*AngularDeflection/BellowNum: -0.001: (n-1)*AngularDeflection/BellowNum - Phi2;
ThetaUpper = pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2 - Phi1: 0.001 : pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2 + Phi1;
ThetaRightLower = (n)*AngularDeflection/BellowNum + Phi2: -0.001: (n)*AngularDeflection/BellowNum;


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

bX = [bx_1 bx_2 bx_3 bx_4 bx_5];
bY = [by_1 by_2 by_3 by_4 by_5];

plot(bx_1,by_1, '-k');
plot(bx_2,by_2, '-k');
plot(bx_3,by_3, '-k');
plot(bx_4,by_4, '-k');
plot(bx_5,by_5, '-k');


% Outer Bellow
LeftLowerArcXCenter = DeformedCenterX + ((EffectiveLength)/AngularDeflection  + OuterRadius + ri + r2) * cos(pi + (n-1)*AngularDeflection/BellowNum);
LeftLowerArcYCenter = DeformedCenterY + ((EffectiveLength)/AngularDeflection  + OuterRadius + ri + r2) * sin(pi + (n-1)*AngularDeflection/BellowNum);
UpperArcXCenter = DeformedCenterX + ((EffectiveLength)/AngularDeflection + OuterRadius + ro - r1) * cos(pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2);
UpperArcYCenter = DeformedCenterY + ((EffectiveLength)/AngularDeflection + OuterRadius + ro - r1) * sin(pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2);
RightLowerArcXCenter = DeformedCenterX + ((EffectiveLength)/AngularDeflection + OuterRadius + ri + r2) * cos(pi + (n)*AngularDeflection/BellowNum);
RightLowerArcYCenter = DeformedCenterY + ((EffectiveLength)/AngularDeflection + OuterRadius + ri + r2) * sin(pi + (n)*AngularDeflection/BellowNum);
ThetaLeftLower = (n-1)*AngularDeflection/BellowNum - Phi2: 0.001: (n-1)*AngularDeflection/BellowNum;
ThetaUpper = pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2 + Phi1: -0.001 : pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2 - Phi1;
ThetaRightLower = (n)*AngularDeflection/BellowNum: 0.001: (n)*AngularDeflection/BellowNum + Phi2;


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

bX = [bX bx_5 bx_4 bx_3 bx_2 bx_1];
bY = [bY by_5 by_4 by_3 by_2 by_1];

plot(bx_1,by_1, '-k');
plot(bx_2,by_2, '-k');
plot(bx_3,by_3, '-k');
plot(bx_4,by_4, '-k');
plot(bx_5,by_5, '-k');
    
fill(bX, bY, [0.9 0.85 0],'LineStyle','none', 'FaceAlpha', 0.5);

for n = 2:1:BellowNum-1
    % inner bellow
    LeftLowerArcXCenter = DeformedCenterX + ((EffectiveLength)/AngularDeflection  + OuterRadius + ri + r2) * cos(pi + (n-1)*AngularDeflection/BellowNum);
    LeftLowerArcYCenter = DeformedCenterY + ((EffectiveLength)/AngularDeflection  + OuterRadius + ri + r2) * sin(pi + (n-1)*AngularDeflection/BellowNum);
    UpperArcXCenter = DeformedCenterX + ((EffectiveLength)/AngularDeflection + OuterRadius + ro - r1) * cos(pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2);
    UpperArcYCenter = DeformedCenterY + ((EffectiveLength)/AngularDeflection + OuterRadius + ro - r1) * sin(pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2);
    RightLowerArcXCenter = DeformedCenterX + ((EffectiveLength)/AngularDeflection + OuterRadius + ri + r2) * cos(pi + (n)*AngularDeflection/BellowNum);
    RightLowerArcYCenter = DeformedCenterY + ((EffectiveLength)/AngularDeflection + OuterRadius + ri + r2) * sin(pi + (n)*AngularDeflection/BellowNum);
    ThetaLeftLower = (n-1)*AngularDeflection/BellowNum: -0.001: (n-1)*AngularDeflection/BellowNum - Phi2;
    ThetaUpper = pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2 - Phi1: 0.001 : pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2 + Phi1;
    ThetaRightLower = (n)*AngularDeflection/BellowNum + Phi2: -0.001: (n)*AngularDeflection/BellowNum;
    
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
    
    bX = [bx_1 bx_2 bx_3 bx_4 bx_5];
    bY = [by_1 by_2 by_3 by_4 by_5];
    
    plot(bx_1,by_1, '-k');
    plot(bx_2,by_2, '-k');
    plot(bx_3,by_3, '-k');
    plot(bx_4,by_4, '-k');
    plot(bx_5,by_5, '-k');
    
    %outer bellow
    LeftLowerArcXCenter = DeformedCenterX + ((EffectiveLength)/AngularDeflection  + OuterRadius + ri + r2) * cos(pi + (n-1)*AngularDeflection/BellowNum);
    LeftLowerArcYCenter = DeformedCenterY + ((EffectiveLength)/AngularDeflection  + OuterRadius + ri + r2) * sin(pi + (n-1)*AngularDeflection/BellowNum);
    UpperArcXCenter = DeformedCenterX + ((EffectiveLength)/AngularDeflection + OuterRadius + ro - r1) * cos(pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2);
    UpperArcYCenter = DeformedCenterY + ((EffectiveLength)/AngularDeflection + OuterRadius + ro - r1) * sin(pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2);
    RightLowerArcXCenter = DeformedCenterX + ((EffectiveLength)/AngularDeflection + OuterRadius + ri + r2) * cos(pi + (n)*AngularDeflection/BellowNum);
    RightLowerArcYCenter = DeformedCenterY + ((EffectiveLength)/AngularDeflection + OuterRadius + ri + r2) * sin(pi + (n)*AngularDeflection/BellowNum);
    ThetaLeftLower = (n-1)*AngularDeflection/BellowNum - Phi2: 0.001: (n-1)*AngularDeflection/BellowNum;
    ThetaUpper = pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2 + Phi1: -0.001 : pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2 - Phi1;
    ThetaRightLower = (n)*AngularDeflection/BellowNum: 0.001: (n)*AngularDeflection/BellowNum + Phi2;
    
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
    
    bX = [bX bx_5 bx_4 bx_3 bx_2 bx_1];
    bY = [bY by_5 by_4 by_3 by_2 by_1];

    plot(bx_1,by_1, '-k');
    plot(bx_2,by_2, '-k');
    plot(bx_3,by_3, '-k');
    plot(bx_4,by_4, '-k');
    plot(bx_5,by_5, '-k');

    fill(bX, bY, [0.9 0.85 0],'LineStyle','none', 'FaceAlpha', 0.5);
end

n = BellowNum;
% inner bellow
LeftLowerArcXCenter = DeformedCenterX + ((EffectiveLength)/AngularDeflection  + OuterRadius + ri + r2) * cos(pi + (n-1)*AngularDeflection/BellowNum);
LeftLowerArcYCenter = DeformedCenterY + ((EffectiveLength)/AngularDeflection  + OuterRadius + ri + r2) * sin(pi + (n-1)*AngularDeflection/BellowNum);
UpperArcXCenter = DeformedCenterX + ((EffectiveLength)/AngularDeflection + OuterRadius + ro - r1) * cos(pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2);
UpperArcYCenter = DeformedCenterY + ((EffectiveLength)/AngularDeflection + OuterRadius + ro - r1) * sin(pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2);
RightLowerArcXCenter = DeformedCenterX + ((EffectiveLength)/AngularDeflection + OuterRadius + ri + r2) * cos(pi + (n)*AngularDeflection/BellowNum);
RightLowerArcYCenter = DeformedCenterY + ((EffectiveLength)/AngularDeflection + OuterRadius + ri + r2) * sin(pi + (n)*AngularDeflection/BellowNum);
ThetaLeftLower = (n-1)*AngularDeflection/BellowNum: -0.001: (n-1)*AngularDeflection/BellowNum - Phi2;
ThetaUpper = pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2 - Phi1: 0.001 : pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2 + Phi1;
ThetaRightLower = (n)*AngularDeflection/BellowNum + Phi2: -0.001: (n)*AngularDeflection/BellowNum;


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
bX = [bx_1 bx_2 bx_3 bx_4 bx_5];
bY = [by_1 by_2 by_3 by_4 by_5];

plot(bx_1,by_1, '-k');
plot(bx_2,by_2, '-k');
plot(bx_3,by_3, '-k');
plot(bx_4,by_4, '-k');
plot(bx_5,by_5, '-k');

% outer bellow
LeftLowerArcXCenter = DeformedCenterX + ((EffectiveLength)/AngularDeflection  + OuterRadius + ri + r2) * cos(pi + (n-1)*AngularDeflection/BellowNum);
LeftLowerArcYCenter = DeformedCenterY + ((EffectiveLength)/AngularDeflection  + OuterRadius + ri + r2) * sin(pi + (n-1)*AngularDeflection/BellowNum);
UpperArcXCenter = DeformedCenterX + ((EffectiveLength)/AngularDeflection + OuterRadius + ro - r1) * cos(pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2);
UpperArcYCenter = DeformedCenterY + ((EffectiveLength)/AngularDeflection + OuterRadius + ro - r1) * sin(pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2);
RightLowerArcXCenter = DeformedCenterX + ((EffectiveLength)/AngularDeflection + OuterRadius + ri + r2) * cos(pi + (n)*AngularDeflection/BellowNum);
RightLowerArcYCenter = DeformedCenterY + ((EffectiveLength)/AngularDeflection + OuterRadius + ri + r2) * sin(pi + (n)*AngularDeflection/BellowNum);
ThetaLeftLower = (n-1)*AngularDeflection/BellowNum - Phi2: 0.001: (n-1)*AngularDeflection/BellowNum;
ThetaUpper = pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2 + Phi1: -0.001 : pi + (2*(n-1)+1)*AngularDeflection/BellowNum/2 - Phi1;
ThetaRightLower = (n)*AngularDeflection/BellowNum: 0.001: (n)*AngularDeflection/BellowNum + Phi2;


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

bx_6 = [(RightLowerArcXCenter + (r2-t)*cos(BellowNum*AngularDeflection/BellowNum)) (RightLowerArcXCenter + (r2-t)*cos(BellowNum*AngularDeflection/BellowNum) + t*sin(AngularDeflection))];
by_6 = [(RightLowerArcYCenter + (r2-t)*sin(BellowNum*AngularDeflection/BellowNum)) (RightLowerArcYCenter + (r2-t)*sin(BellowNum*AngularDeflection/BellowNum) - t*cos(AngularDeflection))];

bX = [bX bx_6(2) bx_6(1) bx_5 bx_4 bx_3 bx_2 bx_1];
bY = [bY by_6(2) by_6(1) by_5 by_4 by_3 by_2 by_1];

plot(bx_1,by_1, '-k');
plot(bx_2,by_2, '-k');
plot(bx_3,by_3, '-k');
plot(bx_4,by_4, '-k');
plot(bx_5,by_5, '-k');
fill(bX, bY, [0.9 0.85 0],'LineStyle','none', 'FaceAlpha', 0.5);





%end cap
cx_1 = [box_7(1) (box_7(1)+t*sin(AngularDeflection))]; % bottom side short t
cy_1 = [boy_7(1) (boy_7(1)-t*cos(AngularDeflection))];
cx_2 = [(box_7(2)+t*sin(AngularDeflection)) box_7(2)]; % bellow side short t
cy_2 = [(boy_7(2)-t*cos(AngularDeflection)) boy_7(2)];
cx_3 = [cx_1(2) cx_2(1)]; % 
cy_3 = [cy_1(2) cy_2(1)];

cX = [cx_1 cx_3 cx_2];
cY = [cy_1 cy_3 cy_2];

plot(cx_1, cy_1, '-b');
plot(cx_2, cy_2, '-b');
plot(cx_3, cy_3, '-b');

fill(cX,cY,[0 0.4470 0.7410],'LineStyle','none', 'FaceAlpha', 0.5)


% angle display

ax_1 = [((cx_3(1)+cx_3(2))/2 - 10*sin(AngularDeflection)) (cx_3(1)+cx_3(2))/2 ((cx_3(1)+cx_3(2))/2 + 10*sin(AngularDeflection))];
ay_1 = [((cy_3(1)+cy_3(2))/2 + 10*cos(AngularDeflection)) (cy_3(1)+cy_3(2))/2 ((cy_3(1)+cy_3(2))/2 - 10*cos(AngularDeflection))];
ax_2 = [(cx_3(1)+cx_3(2))/2 ((cx_3(1)+cx_3(2))/2 + 10)];
ay_2 = [(cy_3(1)+cy_3(2))/2 (cy_3(1)+cy_3(2))/2];

theta = 0:0.01:AngularDeflection;
ax_3 = (cx_3(1)+cx_3(2))/2 + (2*t)*cos(theta);
ay_3 = (cy_3(1)+cy_3(2))/2 + (2*t)*sin(theta);

plot(ax_1,ay_1, '--k');
plot(ax_2,ay_2, '--k');
plot(ax_3,ay_3, '-k');

axis equal
grid on
box on



saveas(DeformedSchematicDiagram, 'DeformedSchematicDiagram.jpg')
set(DeformedSchematicDiagram, 'visible', 'on'); 
hold off
close(DeformedSchematicDiagram)
end