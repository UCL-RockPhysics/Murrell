# Murrell functions, a set of functions used to process data from the 'Murrell'
# gas medium triaxial, RIPL, UCL
# Structure to store raw data

# Murell file directory interaction
function M_file(P,col)
    if Sys.iswindows()
        xf = XLSX.readxlsx("C:\\Users\\cwaha\\Dropbox\\My PC (DESKTOP-8JF2H49)\\Documents\\UCL\\raw_lab data\\log.xlsx")
    else
        xf = XLSX.readxlsx("/Users/christopherharbord/Dropbox/My PC (DESKTOP-8JF2H49)/Documents/UCL/raw_lab data/log.xlsx")
    end
    sh = xf[XLSX.sheetnames(xf)[1]]
    info = sh[:]
    P[:I]=[info[col,4],info[col,5],info[col,6]] # Indice at HP
    P[:T] = info[col,7] # Temperature of experiment
    P[:εr] = info[col,8] # Strain rate of test
    P[:K] = info[col,9]*1e-3 # Machine compliance [m kN^-1]
    P[:L] = info[col,10] # Sample length
    P[:d] = info[col,11]*1e-3 # Diameter of sample [m]
    fil = info[col,1]
    if Sys.iswindows()
        P[:fid] = "C:\\Users\\cwaha\\Dropbox\\My PC (DESKTOP-8JF2H49)\\Documents\\UCL\\raw_lab data\\"*fil*".tdms"
    else
        P[:fid] = "/Users/christopherharbord/Dropbox/My PC (DESKTOP-8JF2H49)/Documents/UCL/raw_lab data/"*fil*".tdms"
    end
end
# Murrell file reader
function M_read(P)
    #Read data from .tdms file
    tdmsIN = readtdms(P[:fid])
    P[:t_s] = tdmsIN.groups["Numeric"]["TimeStamp"].data
    P[:F_kN] = tdmsIN.groups["Numeric"]["Load"].data
    P[:Ua_mm] = tdmsIN.groups["Numeric"]["Displacement"].data
    P[:U1_mm] = tdmsIN.groups["Numeric"]["LVDT 1"].data
    P[:U2_mm] = tdmsIN.groups["Numeric"]["LVDT 2"].data
    P[:Pc1_MPa] = tdmsIN.groups["Numeric"]["Pc 700MPa"].data
    P[:Pc2_MPa] = tdmsIN.groups["Numeric"]["PC 1400 MPa"].data
    P[:Pf_MPa] = tdmsIN.groups["Numeric"]["Pore Pressure"].data
    P[:PpVol_mm3] = tdmsIN.groups["Numeric"]["PP vol"].data

    # Perform some data reduction
    P[:t_s_c] = P[:t_s].-P[:t_s][P[:I][1]]
    P[:F_kN_c] = P[:F_kN] .-P[:F_kN][P[:I][1]]
    P[:U_mm_c] = ((P[:U1_mm].+P[:U2_mm]).-(P[:U1_mm][P[:I][1]]+P[:U2_mm][P[:I][1]]))./2
    P[:U_mm_fc] = P[:U_mm_c] .-(P[:F_kN_c]*P[:K])
    P[:ε] = P[:U_mm_fc]./P[:L]
    P[:F_kN_j] = P[:F_kN_c] .-JR(P)*1e-3
    P[:σ_MPa] = P[:F_kN_c]./(0.25π*P[:d]^2) .*1e-3
    P[:σ_MPa_j] = P[:F_kN_j]./(0.25π*P[:d]^2) .*1e-3
end
# Hardening modulus calculations
function H_mod(P, ds)
    ε1 = Array{Float64,1}(undef,1)
    σ1 = Array{Float64,1}(undef,1)
    t1 = Array{Float64,1}(undef,1)
    ε = P[:ε][P[:I][1]:P[:I][2]]
    σ = P[:σ_MPa][P[:I][1]:P[:I][2]]
    t = P[:t_s][P[:I][1]:P[:I][2]]
    for i in 1:Int(floor(length(P[:ε][P[:I][1]:P[:I][2]])/ds))-1
        push!(ε1,sum(ε[Int(i*ds-ds+1):Int(i*ds)])/ds)
        push!(σ1,sum(σ[Int(i*ds-ds+1):Int(i*ds)])/ds)
        push!(t1,t[Int(i*ds)+Int(ceil(ds/2))])
    end
    deleteat!(ε1,1)
    deleteat!(σ1,1)
    deleteat!(t1,1)
    dσ = differentiate(t1,σ1,TotalVariation(),0.2,5e-3, maxit=2000)
    dε = differentiate(t1,ε1,TotalVariation(),0.2,5e-3, maxit=2000)
    P[:h] = dσ./dε*1e-3 # Calculate the tangent modulus in GPa
    P[:ε_h] = ε1
    h_max = maximum(P[:h])
    II = findfirst(P[:h] .== h_max)
    P[:E_GPa] = sum(P[:h][P[:h].>0.95h_max])/sum(P[:h].>0.95h_max)
    P[:H_GPa] = sum(P[:h][P[:h].<0.4h_max])/sum(P[:h].<0.4h_max)
    σ2 = σ1[II:end]
    h = P[:h][II:end]
    P[:YSa_MPa] = σ2[findfirst(h.<0.9h_max)]
    P[:YSb_MPa] = σ2[findfirst(h.<0.7h_max)]

end
# SFP log reader
function SFP_read(fil)
    SFP = Dict()
    path="/Users/christopherharbord/Dropbox/My PC (DESKTOP-8JF2H49)/Documents/UCL/Furnace_calibration/SFP_logging/"*fil
    dat = readdlm(path)
    headers = dat[:,1]
    I = findall(x->x==headers[1],dat[:,1])
    SFP[:t] = dat[:,1]
    SFP[:PF] = dat[:,3]
    SFP[:S1] = dat[:,4]
    SFP[:S2] = dat[:,5]
    SFP[:S3] = dat[:,6]
    SFP[:TC1] = dat[:,7]
    SFP[:TC2] = dat[:,8]
    close(dat)
    deleteat!(SFP[:t],I)
    deleteat!(SFP[:PF],I)
    deleteat!(SFP[:S1],I)
    deleteat!(SFP[:S2],I)
    deleteat!(SFP[:S3],I)
    deletat!(SFP[:TC1],I)
    deletat!(SFP[:TC2],I)
end
# Jacket Rheology function
function JR(P)
    #Inputs: passed as a dictionary
    #T.......Temperature (°C)
    #ε.......Sample strain (-)
    #ε̇.......Strain rate during test (s⁻¹)
    #d.......Sample diameter (m)
    #L.......Sample length (m)
    EaJ = 197000 # Copper activation enthalpy
    n1J = 4.8 # Empirical factor 1
    n2J = 0.22 # Empirical factor 2
    ε̇j = P[:εr]*(P[:d]/P[:L]) # Jacket strain rate
    Jr = 2*10^((log10(ε̇j*exp(EaJ/(8.3145*(P[:T]+278)))))/n1J-n2J) # Copper flow stress at experiment conditions
    Ja = π*P[:d]*P[:L] # Jacket area
    JR = Jr.*Ja.*(1 .+(P[:ε].*P[:d]/P[:L])) # Force due to jacket assuming linear increase due to incremental strain
end
