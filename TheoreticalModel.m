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
ContactWidth = AverageRadius * 2;


r1 = UpperArcRadius + WallThickness/2; % initial value of upper arc radius
r2 = LowerArcRadius - WallThickness/2; % initial value of lower arc radius
Phi1 = pi/2; % initial value of upper arc angle
Phi2 = pi/2; % initial value of lower arc angle
f = flank; % initial value of flank
Beta = 0; % initial value of flank titled angle
ro = OuterRadius;
ri = InnerRadius;
rm = AverageRadius;
t = WallThickness;
L0 = EffectiveLength;
D = ContactWidth;
E = YoungModulus/1000; % N/mm^2 = MP

AngularDeflection = 0;
P = 10;
Pressure = P / 1000; % N/mm^2 = 1000 KPa
p = 100/10^6;
step = Pressure/p;
b = 3.8314 * 1/BellowNum/(InnerRadius/AverageRadius)/(E)^(1/3);

i = 1;
PressureData(i) = 0;
TheoreticalDeflection(i) = AngularDeflection;
InnerRadiusData(i) = ri;
AverageRadiusData(i) = rm;
OuterRadiusData(i) = ro;
R1Data(i) = r1;
R2Data(i) = r2;
PHI1Data(i) = Phi1;
PHI2Data(i) = Phi2;

while 1
        I = sqrt(pi) * rm * t^3 / 12; 
        
        
            % Beam Theory          
              
            dx_BC = ((f^3*p*r2)/6 + (f^4*p*cos(Beta))/8 + (Phi2*f^2*p*r2^2)/2 + (Phi2^2*f*p*r2^3)/2 - (f^3*p*r2*cos(Phi2))/6 - (Phi2*f^2*p*r2^2*cos(Phi2))/2 - (Phi2^2*f*p*r2^3*cos(Phi2))/2 + (Phi2^2*f^2*p*r2^2*cos(Beta))/4 + (Phi2*f^3*p*r2*cos(Beta))/3)/(E*I);
            dy_BC = (f*p*(3*f^3*sin(Beta) + 4*f^2*r2*sin(Phi2) + 12*Phi2^2*r2^3*sin(Phi2) + 8*Phi2*f^2*r2*sin(Beta) + 12*Phi2*f*r2^2*sin(Phi2) + 6*Phi2^2*f*r2^2*sin(Beta)))/(24*E*I);
            dx_AB = (r1*((p*r1*(2*Phi1*r1^2 + f^2*sin(Phi1) - 2*r1^2*sin(Phi1) + Phi2^2*r2^2*sin(Phi1) - 4*f*r1*(cos(Phi1)/2 - 1/2) - 4*Phi2*r1*r2*(cos(Phi1)/2 - 1/2) + 2*Phi2*f*r2*sin(Phi1)))/2 - (Phi1^3*p*r1^3*cos(Phi1))/6 + (Phi1^3*p*r1^2*r2)/6 + (Phi1*p*r2*(f + Phi1*r1 + Phi2*r2)^2)/2 + (Phi1^2*p*r1^2*cos(Phi1)*(f + Phi1*r1 + Phi2*r2))/2 + (Phi1^3*f*p*r1^2*cos(Beta))/6 - (Phi1^2*p*r1*r2*(f + Phi1*r1 + Phi2*r2))/2 - (Phi1^3*p*r1^2*r2*cos(Phi2))/6 + (Phi1*f*p*cos(Beta)*(f + Phi1*r1 + Phi2*r2)^2)/2 - (Phi1*p*r1*cos(Phi1)*(f + Phi1*r1 + Phi2*r2)^2)/2 - (Phi1*p*r2*cos(Phi2)*(f + Phi1*r1 + Phi2*r2)^2)/2 - (Phi1^2*f*p*r1*cos(Beta)*(f + Phi1*r1 + Phi2*r2))/2 + (Phi1^2*p*r1*r2*cos(Phi2)*(f + Phi1*r1 + Phi2*r2))/2))/(E*I);
            dy_AB = (r1*((Phi1^3*p*r1^3*sin(Phi1))/6 - (p*r1*(2*f^2*sin(Phi1/2)^2 - 4*r1^2*sin(Phi1/2)^2 + Phi1^2*r1^2 + 2*Phi2^2*r2^2*sin(Phi1/2)^2 + 2*Phi1*f*r1 - 2*f*r1*sin(Phi1) + 4*Phi2*f*r2*sin(Phi1/2)^2 - 2*Phi2*r1*r2*sin(Phi1) + 2*Phi1*Phi2*r1*r2))/2 - (Phi1^2*p*r1^2*sin(Phi1)*(f + Phi1*r1 + Phi2*r2))/2 + (Phi1^3*f*p*r1^2*sin(Beta))/6 + (Phi1^3*p*r1^2*r2*sin(Phi2))/6 + (Phi1*f*p*sin(Beta)*(f + Phi1*r1 + Phi2*r2)^2)/2 + (Phi1*p*r1*sin(Phi1)*(f + Phi1*r1 + Phi2*r2)^2)/2 + (Phi1*p*r2*sin(Phi2)*(f + Phi1*r1 + Phi2*r2)^2)/2 - (Phi1^2*f*p*r1*sin(Beta)*(f + Phi1*r1 + Phi2*r2))/2 - (Phi1^2*p*r1*r2*sin(Phi2)*(f + Phi1*r1 + Phi2*r2))/2))/(E*I);
            dx_CD = (r2*((Phi2^3*p*r2^3)/6 - (p*r2*(r2^2*sin(Phi2)*(Phi2^2 - 2) + 2*Phi2*r2^2*cos(Phi2)))/2))/(E*I);
            dy_CD = (p*r2^4*(2*cos(Phi2) + 2*Phi2*sin(Phi2) - Phi2^2*cos(Phi2) - 2))/(2*E*I);
            dx = (dx_AB + dx_BC + dx_CD)/b;
            dy = (dy_AB + dy_BC + dy_CD)/b;
            
            syms r1_new Phi1_new r2_new Phi2_new f_new Beta_new
            eqns = [r1_new*sin(Phi1_new)+r2_new*sin(Phi2_new)+f_new*sin(Beta_new) == r1*sin(Phi1) + r2*sin(Phi2) + f*sin(Beta) + dx , r1_new*(1-cos(Phi1_new))+r2_new*(1-cos(Phi2_new))+f_new*cos(Beta_new) == r1*(1-cos(Phi1)) + r2*(1-cos(Phi2)) + f*cos(Beta) - dy, r1_new*Phi1_new/f_new == r1*Phi1/f, r2_new*Phi2_new/f_new == r2*Phi2/f, r1_new/r2_new == r1/r2, Phi1_new+Beta_new == pi/2];
            S = vpasolve(eqns, [r1_new Phi1_new r2_new Phi2_new f_new Beta_new], [r1, Phi1, r2, Phi2, f, Beta]);
            r1 = double(S.r1_new);
            Phi1 = double(S.Phi1_new);
            r2 = double(S.r2_new);
            Phi2 = double(S.Phi2_new);
            f = double(S.f_new);
            Beta = double(S.Beta_new);
            
            ri = ro - ((r1-t/2)*(1-cos(Phi1))+(r2+t/2)*(1-cos(Phi2))+f*cos(Beta));
            rm = ri + (r2+t/2)*(1-cos(Phi2)) + f/2*cos(Beta);
%             rm = (ri+ro)/2;
            


        AngularDeflection = (((r1-t/2)*sin(Phi1)+(r2+t/2)*sin(Phi2)+f*sin(Beta))*BellowNum*2 - L0)/(sqrt(pi)*rm/2+ro+2*t);       % in rad unit
        if (r1 <= 0)||(r2 <= 0)||(Phi1 <= 0.1)||(Phi2 <= 0.1)||(AngularDeflection >= pi)
            break;
        end
        i = i + 1;
        PressureData(i) = p*(i-1)*1000;
        TheoreticalDeflection(i) = AngularDeflection;
        InnerRadiusData(i) = ri;
        AverageRadiusData(i) = rm;
        OuterRadiusData(i) = ro;
        R1Data(i) = r1;
        R2Data(i) = r2;
        PHI1Data(i) = Phi1;
        PHI2Data(i) = Phi2;
end

DeformedSchematicDiagram(EffectiveLength,AverageRadius,BellowNum,InnerRadius, WallThickness, ri, rm, ro, r1, r2, Phi1, Phi2, AngularDeflection);
end