function [PressureData, TheoreticalDeflection, InnerRadiusData, AverageRadiusData, OuterRadiusData, R1Data, R2Data, PHI1Data, PHI2Data] = TheoreticalModel(EffectiveLength, AverageRadius, BellowNum, InnerRadius, WallThickness, YoungModulus)


PressureData = [];
TheoreticalDeflection = [];
InnerRadiusData = [];
AverageRadiusData = [];
OuterRadiusData = [];
R1Data = [];
R2Data = [];
PHI1Data = [];
PHI2Data = [];


OuterRadius = AverageRadius * 2 - InnerRadius;
UpperArcRadius = EffectiveLength/BellowNum/2/2;
LowerArcRadius = EffectiveLength/BellowNum/2/2;
flank = OuterRadius - InnerRadius - UpperArcRadius - LowerArcRadius;


R1 = (UpperArcRadius + WallThickness/2)/1000; % initial value of upper arc radius
R2 = (LowerArcRadius - WallThickness/2)/1000; % initial value of lower arc radius

f = flank/1000; % initial value of flank

ro = OuterRadius/1000;
ri = InnerRadius/1000;
rm = AverageRadius/1000;
t = WallThickness/1000;
N = BellowNum;
E = YoungModulus*1000;



AngularDeflection = 0;
C = 1.6082e+04/(ri/rm)^(-0.3)/N^(0.9)/E^(0.6);
P = 100;


i = 1;
PressureData(i) = 0;
TheoreticalDeflection(i) = AngularDeflection;
InnerRadiusData(i) = ri;
AverageRadiusData(i) = rm;
OuterRadiusData(i) = ro;
R1Data(i) = R1;
R2Data(i) = R2;
PHI1Data(i) = pi/2;
PHI2Data(i) = pi/2;

while 1
    I = t^3 / 12; 
    p = P*i;    
        
    F = p*(ro-ri)/2 + p*ri/2; 
    MA = (1/12).*(2.*f+pi.*R1+pi.*R2).^(-1).*(12.*f.^2.*F+(-4).*f.^3.*p+24.*f.*F.*R1+(-12).*f.^2.*p.*R1+(-24).*F.*R1.^2+(-12).*f.*p.*R1.^2+12.*F.*pi.*R1.^2+24.*p.*R1.^3+(-9).*p.*pi.*R1.^3+12.*f.*F.*pi.*R2+(-6).*f.^2.*p.*pi.*R2+12.*F.*pi.*R1.*R2+(-12).*f.*p.*pi.*R1.*R2+(-6).*p.*pi.*R1.^2.*R2+24.*F.*R2.^2+(-24).*f.*p.*R2.^2+(-24).*p.*R1.*R2.^2+(-3).*p.*pi.*R2.^3);

    dx = (-1/3).*f.^3.*F+(1/2).*f.^2.*MA+(1/8).*f.^4.*p+(-1).*f.^2.*F.*R1+f.*MA.*R1+(1/2).*f.^3.*p.*R1+(-1).*f.*F.*R1.^2+(3/4).*f.^2.*p.*R1.^2+(1/2).*f.*p.*R1.^3+(1/24).*R1.^2.*(12.*MA.*((-2)+pi)+R1.*((-6).*F.*((-8)+3.*pi)+p.*((-44)+15.*pi).*R1))+(1/24).*R2.*(6.*pi.*(f+R1).*(2.*MA+(f+R1).*((-2).*F+p.*(f+R1)))+12.*(2.*MA+(f+R1).*((-4).*F+3.*p.*(f+R1))).*R2+3.*pi.*((-2).*F+3.*p.*(f+R1)).*R2.^2+8.*p.*R2.^3);
    DX = -dx/E/I/C;
    dy = (1/6).*R1.^2.*(6.*MA+R1.*((-3).*F+p.*R1))+(-1/2).*f.^2.*F.*R2+f.*MA.*R2+(1/6).*f.^3.*p.*R2+(-1).*f.*F.*R1.*R2+(1/2).*f.^2.*p.*R1.*R2+(1/2).*f.*p.*R1.^2.*R2+(1/24).*R2.*(6.*pi.*R1.*(2.*MA+(f+R1).*((-2).*F+p.*(f+R1)))+6.*(((-2).*f.*F+2.*MA+f.^2.*p).*((-2)+pi)+(-2).*(F+(-1).*f.*p).*pi.*R1+p.*(2+pi).*R1.^2).*R2+3.*((-4).*F+4.*f.*p+p.*(4+pi).*R1).*R2.^2+p.*((-4)+3.*pi).*R2.^3);
    DY = -dy/E/I/C;

    AngularDeflection = DX*N/(2*ro+2*t);  % in rad unit

    syms R1_new Phi1_new R2_new Phi2_new f_new Beta_new
    eqns = [R1_new*sin(Phi1_new)+R2_new*sin(Phi2_new)+f_new*sin(Beta_new) == R1 + R2 + DX/2 , R1_new*(1-cos(Phi1_new))+R2_new*(1-cos(Phi2_new))+f_new*cos(Beta_new) == R1 + R2 + f - DY/2, R1_new*Phi1_new/f_new == R1*pi/2/f, R2_new*Phi2_new/f_new == R2*pi/2/f, R1_new/R2_new == R1/R2, Phi1_new+Beta_new == pi/2];
    S = vpasolve(eqns, [R1_new Phi1_new R2_new Phi2_new f_new Beta_new], [R1, pi/2, R2, pi/2, f, 0]);
    

                             
    if (double(S.R1_new) <= 0)||(double(S.R2_new) <= 0)||(double(S.Phi1_new) <= 0.1)||(double(S.Phi2_new) <= 0.1)||(AngularDeflection >= pi)||(double(S.Beta_new)>=pi)
        break;
    end
    i = i + 1;

    PressureData(i) = p/1000;
    TheoreticalDeflection(i) = AngularDeflection;
    InnerRadiusData(i) = ri*1000;
    AverageRadiusData(i) = rm*1000;
    OuterRadiusData(i) = (ro-DY/2)*1000;
    R1Data(i) = (double(S.R1_new)-t/2)*1000;
    R2Data(i) = (double(S.R2_new)+t/2)*1000;
    PHI1Data(i) = double(S.Phi1_new);
    PHI2Data(i) = double(S.Phi2_new);
end

DeformedSchematicDiagram(EffectiveLength, AverageRadius, BellowNum, InnerRadius, WallThickness, InnerRadius, AverageRadius, (ro-DY/2)*1000, (double(S.R1_new)-t/2)*1000, (double(S.R2_new)+t/2)*1000, double(S.Phi1_new), double(S.Phi2_new), AngularDeflection);
end