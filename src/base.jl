"""
    M_file(P, col::int64)

    Generate filepath and extract indices from spreadsheet log.xlsx

    Arguments
    - 'col::Int64' : column of file of interest in accompanying spread sheet
"""
function M_file(col::Int64)
    if Sys.iswindows()
        xf = XLSX.readxlsx("C:\\Users\\cwaha\\Dropbox\\My PC (DESKTOP-8JF2H49)\\Documents\\UCL\\raw_lab data\\log.xlsx")
    else
        xf = XLSX.readxlsx("/Users/christopherharbord/Dropbox/My PC (DESKTOP-8JF2H49)/Documents/UCL/raw_lab data/log.xlsx")
    end
    sh = xf[XLSX.sheetnames(xf)[1]]
    info = sh[:]
    exp_info = Dict()
    exp_info[:I]=[info[col,4],info[col,5],info[col,6]] # Indice at HP
    exp_info[:P_exp] = info[col,2]
    exp_info[:T] = info[col,7] # Temperature of experiment
    exp_info[:εr] = info[col,8] # Strain rate of test
    exp_info[:K_mm_kN] = info[col,9] # Machine compliance [mm kN^-1]
    exp_info[:L_mm] = info[col,10] # Sample length
    exp_info[:d_mm] = info[col,11] # Diameter of sample [m]
    fil = info[col,1]
    if Sys.iswindows()
        fid = "C:\\Users\\cwaha\\Dropbox\\My PC (DESKTOP-8JF2H49)\\Documents\\UCL\\raw_lab data\\"*fil*".tdms"
    else
        fid = "/Users/christopherharbord/Dropbox/My PC (DESKTOP-8JF2H49)/Documents/UCL/raw_lab data/"*fil*".tdms"
    end
    return fid, exp_info
end
"""
    M_read(P)

    Open and read contents of .tdms file
    #Arguments
    - fid : file path string
    return P
"""
function M_read(fid::String)
    tdmsIN = readtdms(fid)
    P = Dict()
    P[:t_s] = tdmsIN.groups["Numeric"]["TimeStamp"].data
    P[:F_kN] = tdmsIN.groups["Numeric"]["Load"].data
    P[:Ua_mm] = tdmsIN.groups["Numeric"]["Displacement"].data
    P[:U1_mm] = tdmsIN.groups["Numeric"]["LVDT 1"].data
    P[:U2_mm] = tdmsIN.groups["Numeric"]["LVDT 2"].data
    P[:Pc1_MPa] = tdmsIN.groups["Numeric"]["Pc 700MPa"].data
    P[:Pc2_MPa] = tdmsIN.groups["Numeric"]["PC 1400 MPa"].data
    P[:Pf_MPa] = tdmsIN.groups["Numeric"]["Pore Pressure"].data
    P[:PpVol_mm3] = tdmsIN.groups["Numeric"]["PP vol"].data
    return P
end

function t_conv!(P)
    P[:t_s] = unix2datetime(datetime2unix(DateTime(1904,1,1,0,0,0)) .+(P[:t_s]))
end

"""
    M_reduce!(P, exp_info)

    Reduce Murrell data
    #Arguments
    * P : dictionary containing raw data from the Murrell, must include indices of hit point and sample information
    * exp_info : dictionary containing experimental parameters
"""

function M_reduce!(P,exp_info)
    I1 = exp_info[:I][1]
    P[:t_s_c] = P[:t_s].-P[:t_s][I1]
    P[:F_kN_c] = P[:F_kN] .-P[:F_kN][I1]
    P[:U_mm_c] = ((P[:U1_mm].+P[:U2_mm]).-(P[:U1_mm][I1]+P[:U2_mm][I1]))./2
    P[:U_mm_fc] = P[:U_mm_c] .-(P[:F_kN_c]*exp_info[:K_mm_kN])
    P[:ε] = P[:U_mm_fc]./exp_info[:L_mm]
    P[:Jr] = JR!(P, exp_info)
    P[:F_kN_j] = P[:F_kN_c] .-P[:Jr]
    P[:σ_MPa] = P[:F_kN_c]./(0.25e-6π*exp_info[:d_mm]^2) .*1e-3
    P[:σ_MPa_j] = P[:F_kN_j]./(0.25e-6π*exp_info[:d_mm]^2) .*1e-3
end


"""
    SFP_read()

    Read and return data from program SFP_furnace_control_V2.vi

    #Arguments
    - 'SFP' : empty dictionary to store extracted data
    - 'fil::string' : path to file
"""
function SFP_read(fid)
    path="/Users/christopherharbord/Dropbox/My PC (DESKTOP-8JF2H49)/Documents/UCL/Furnace_calibration/SFP_logging/"*fil
    dat = readdlm(path)
    headers = dat[:,1]
    I = findall(x->x==headers[1],dat[:,1])
    SFP = Dict()
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
"""
    JR(P, exp_info)

    Calculate strength contribution of copper jacket for corrections

    #Arguments
    * 'P' : dictionary containing processed data from the 'Murrell'
    * exp_info : details of experiment
"""
function JR!(P, exp_info)
    EaJ = 197000 # Copper activation enthalpy
    n1J = 4.8 # Empirical factor 1
    n2J = 0.22 # Empirical factor 2
    ε̇j = exp_info[:εr]*(exp_info[:d_mm]/exp_info[:L_mm]) # Jacket strain rate
    Jr = 2*10^((log10(ε̇j*exp(EaJ/(8.3145*(exp_info[:T]+278)))))/n1J-n2J) # Copper flow stress at experiment conditions
    Ja = π*exp_info[:d_mm]*exp_info[:L_mm]*1e-6 # Jacket area
    JR = Jr.*Ja.*(P[:ε]*(exp_info[:d_mm]/exp_info[:L_mm])).*1e-3 # Force due to jacket assuming linear increase due to incremental strain
end

function M_interp!(P, t_UT)
    P[:F_kN_i] = lininterp(P[:t_s],P[:F_kN_j], t_UT)
    P[:U_mm_i] = lininterp(P[:t_s],P[:U_mm_c], t_UT)
    P[:σ_MPa_i] = lininterp(P[:t_s],P[:σ_MPa_j], t_UT)
    P[:σ3_MPa_i] = lininterp(P[:t_s],P[:Pc2_MPa], t_UT)
    P[:ε_i] = lininterp(P[:t_s],P[:ε], t_UT)
end
