using PyCall    
ENV["PYTHON"]="/usr/local/Caskroom/miniforge/base/bin/python3"
cb = pyimport("cbsyst")

function K1(S, T)
    T = T + 273.15  

    
    lnK1 = 2.83655 - (2307.1266 / T) - (1.5529413 * log(T)) - 
           ((0.207608410 + (4.0484 / T)) * sqrt(S)) + 
           (0.0846834 * S) - (0.00654208 * S^(3/2)) + 
           log(1 - (0.001005 * S))

    return exp(lnK1)
end

function K2(S, T)
    T = T + 273.15  
    
    lnK2 = -9.226508 - (3351.6106 / T) - (0.2005743 * log(T)) - ((0.106901773 + 23.9722 / T) * sqrt(S)) + (0.1130822 * S) - (0.00846934 * S^(3/2)) + log(1 - (0.001005 * S))
    
    return exp(lnK2)
end

function Kw(S,T)
    T=T + 273.15
    lnKw= 148.96502 - 13847.26/T - 23.6521*log(T)+
    (118.67/T- 5.977 + 1.0495*log(T))*sqrt(S) - 0.01615*S

    return exp(lnKw)
end

function KB(S, T)
    T = T + 273.15  
    
    
    lnKB = (
        (-8966.90 - 2890.53 * sqrt(S) - (77.942 * S) + (1.728 * S^(3/2)) - (0.0996 * (S^2))) / T
    ) + 148.0248 + (137.1942 * sqrt(S)) + 1.62142 * S - (
        (24.4344 + (25.085 * sqrt(S)) + (0.2474 * S)) * log(T)
    ) + (0.053105 * sqrt(S) * T)
    
    return exp(lnKB)
end

function I_val(S)
    return 19.924 * S / (1000 - 1.005 * S)
end

function KS(S,T)
    T=T+273.15
    I=I_val(S)
    lnKS= -4276.1/T + 141.328 - 23.093*log(T) + (-13856/T + 324.57 - 47.986*log(T))*sqrt(I)+(35474/T - 771.54 + 114.723*log(T))*I  - 2698/T*I^(3/2) + 1776/T*I^2 + log(1-0.001005 *S)
    return exp(lnKS)
end

function KF(S, T, Ks)
    T = T + 273.15
    I = I_val(S)
    ST = 0.02824 * S / 35

    lnKF = (1590.2 / T) - 12.641 + (1.525 * sqrt(I)) +
           log(1 - 0.001005 * S) + log(1 + (ST / Ks))

    return exp(lnKF)
end

function K1P(S,T)
    T=T+273.15
    lnK1P= -4576.752/T +115.525 - 18.453*log(T) + (-106.736/T +0.69171)*sqrt(S) +(-0.65643/T -0.01844)*S
    return exp(lnK1P)
end

function K2P(S,T)
    T=T+273.15
    lnK2P= -8814.715/T + 172.0883- 27.927*log(T) +(-160.34/T + 1.3566)*sqrt(S) + (0.37335/T -0.05778)*S
    return exp(lnK2P)
end

function K3P(S,T)
    T=T+273.15
    lnK3P= -3070.75/T - 18.141 +(17.27039/T + 2.81197)*sqrt(S) +(-44.99486/T - 0.09984)*S
    return exp(lnK3P)
end

function KSi(S,T)
    T=T+273.15
    I=I_val(S)
    lnKSi= -8904.2 /T +  117.385- 19.334*log(T) +(3.5913- 458.79/T)*sqrt(I) + (188.74/T - 1.5998)*I+(0.07871- 12.1652/T)*I^2 +log(1-0.001005*S)
    return exp(lnKSi)
end

function KspC(S, T)
    T = T + 273.15  # Convert Celsius to Kelvin

    logKspC = -171.9065 - 0.077993*T + 2839.319/T +
              71.595 * log10(T) +
              (-0.77712 + 0.0028426*T + 178.34/T) * sqrt(S) -
              0.07711*S + 0.0041249*S^(3/2)
    
    return 10^logKspC  # Return Ksp (not its negative log)
end

function KspA(S, T)
    T = T + 273.15  # Convert Celsius to Kelvin

    logKspA = -171.945 - 0.077993*T + 2903.293/T +
              71.595 * log10(T) +
              (-0.068393 + 0.0017276*T + 88.135/T) * sqrt(S) -
              0.10018*S + 0.0059415*S^(3/2)

    return 10^logKspA  # Return Ksp (not log)
end

function K0(S,T)
    T=T+273.15
    lnK0=9345.17/T- 60.2409 + 23.3585*log(T/100) +S*(0.023517 - 0.00023656*T + 0.0047036*(T/100)^2)
    return exp(lnK0)
end

function calculate_params(pH,pCO2,T=25.0,S=35)
    K1_val=K1(S,T)
    K2_val=K2(S,T)
    Kw_val=Kw(S,T)
    KB_val=KB(S,T)
    KS_val=KS(S,T)
    KF_val=KF(S,T,KS_val)   
    K1P_val=K1P(S,T)
    K2P_val=K2P(S,T)
    K3P_val=K3P(S,T)
    KSi_val=KSi(S,T)
    KspC_val=KspC(S,T)
    KspA_val=KspA(S,T)
    K0_val=K0(S,T)

    K_values = Dict(
    "K1" => K1_val,
    "K2" => K2_val,
    "KB" => KB_val,
    "KW" => Kw_val,
    "KS" => KS_val,
    "KF" => KF_val,
    "KP1" => K1P_val,
    "KP2" => K2P_val,
    "KP3" => K3P_val,
    "KSi" => KSi_val,
    "KspC" => KspC_val,
    "KspA" => KspA_val,
    "K0" => K0_val)

    system=cb.CBsys(pHtot=pH,pCO2=pCO2,T_in=T,S_in=S,Ks=K_values)

    CO2=system["CO2"][1]
    HCO3=system["HCO3"][1]
    CO3=system["CO3"][1]
    OH=system["OH"][1]
    BO3=system["BO3"][1]
    BO4=system["BO4"][1]

    return CO2,HCO3,CO3,OH,BO3,BO4
end




