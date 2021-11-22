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
    if tdmsIN.groups["Numeric"]["TimeStamp"].data[1] < 3.7120499126540936e9
        P[:t_s] = tdmsIN.groups["Numeric"]["TimeStamp"].data
        P[:F_kN] = tdmsIN.groups["Numeric"]["Load"].data
        P[:Ua_mm] = tdmsIN.groups["Numeric"]["Displacement"].data
        P[:U1_mm] = tdmsIN.groups["Numeric"]["LVDT 1"].data
        P[:U2_mm] = tdmsIN.groups["Numeric"]["LVDT 2"].data
        # P[:Pc1_MPa] = tdmsIN.groups["Numeric"]["Pc 700MPa"].data
        P[:Pc2_MPa] = tdmsIN.groups["Numeric"]["PC 1400 MPa"].data
        # P[:Pf_MPa] = tdmsIN.groups["Numeric"]["Pore Pressure"].data
        # P[:PpVol_mm3] = tdmsIN.groups["Numeric"]["PP vol"].data
    else
        A = zeros(6,length(tdmsIN.groups["Numeric"]["TimeStamp"].data))
        A[1,:] = tdmsIN.groups["Numeric"]["TimeStamp"].data
        A[2,:] = tdmsIN.groups["Numeric"]["Load"].data
        A[3,:] = tdmsIN.groups["Numeric"]["Displacement"].data
        A[4,:] = tdmsIN.groups["Numeric"]["LVDT1"].data
        A[5,:] = tdmsIN.groups["Numeric"]["LVDT2"].data
        # P[:Pc1_MPa] = tdmsIN.groups["Numeric"]["PC700"].data
        A[6,:] = tdmsIN.groups["Numeric"]["PC1400"].data
        A = unique(A, dims=2)
        P[:t_s] = A[1,:]
        P[:F_kN] = A[2,:]
        P[:Ua_mm] = A[3,:]
        P[:U1_mm] = A[4,:]
        P[:U2_mm] = A[5,:]
        # P[:Pc1_MPa] = tdmsIN.groups["Numeric"]["PC700"].data
        P[:Pc2_MPa] = A[6,:]
        # P[:Pf_MPa] = tdmsIN.groups["Numeric"]["PF700"].data
        # P[:PpVol_mm3] = tdmsIN.groups["Numeric"]["PFVol"].data
    end
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
    Pc_corr = (P[:Pc2_MPa].-P[:Pc2_MPa][I1]).*0.06
    P[:t_s_c] = P[:t_s].-P[:t_s][I1]
    P[:F_kN_c] = P[:F_kN] .-P[:F_kN][I1].-Pc_corr
    P[:U_mm_c] = ((P[:U1_mm].+P[:U2_mm]).-(P[:U1_mm][I1]+P[:U2_mm][I1]))./2
    P[:U_mm_fc] = P[:U_mm_c] .-(P[:F_kN_c]*exp_info[:K_mm_kN])
    P[:ε] = P[:U_mm_fc]./exp_info[:L_mm]
    P[:Jr] = JR!(P, exp_info)
    P[:F_kN_j] = P[:F_kN_c] .-P[:Jr]
    P[:σ_MPa] = P[:F_kN_c]./((π*0.5e-3*exp_info[:d_mm] .*(1 .+P[:ε])) .^2) .*1e-4
    P[:σ_MPa_j] = P[:F_kN_j]./((π*0.5e-3*exp_info[:d_mm] .*(1 .+P[:ε])) .^2) .*1e-4
end


"""
    SFP_read()

    Read and return data from program SFP_furnace_control_V2.vi

    #Arguments
    - 'SFP' : empty dictionary to store extracted data
    - 'fil::string' : path to file
"""
function SFP_read(fid)
    dat = readdlm(fid)
    headers = dat[:,1]
    I = findall(x->x==headers[1],dat[:,1])
    SFP = Dict()
    SFP[:t] = dat[:,1]
    SFP[:PF] = dat[:,3]
    SFP[:S1] = dat[:,4]
    SFP[:S2] = dat[:,5]
    SFP[:S3] = dat[:,6]
    SFP[:T_lw] = dat[:,7]
    SFP[:T_uw] = dat[:,8]
    SFP[:T_sf] = dat[:,9]
    deleteat!(SFP[:t],I)
    deleteat!(SFP[:PF],I)
    deleteat!(SFP[:S1],I)
    deleteat!(SFP[:S2],I)
    deleteat!(SFP[:S3],I)
    deleteat!(SFP[:T_lw],I)
    deleteat!(SFP[:T_uw],I)
    deleteat!(SFP[:T_sf],I)
    return SFP
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
    Ja = π*exp_info[:d_mm]*1e-3*0.1e-3 # Jacket area
    JR = Jr.*Ja.*(1 .+P[:ε]).*1e-3 # Force due to jacket assuming linear increase due to incremental strain
end
"""
    M_interp!(P, t_UT)

    Interpolate mechanical data to match ultrasonic data pulses
    #Arguments
    * 'P' : dictionary containing processed data from the 'Murrell'
    * exp_info : details of experiment
"""
function M_interp!(P, t_UT)
    P[:F_kN_i] = lininterp(P[:t_s],P[:F_kN_j], t_UT)
    P[:U_mm_i] = lininterp(P[:t_s],P[:U_mm_c], t_UT)
    P[:σ_MPa_i] = lininterp(P[:t_s],P[:σ_MPa_j], t_UT)
    P[:σ3_MPa_i] = lininterp(P[:t_s],P[:Pc2_MPa], t_UT)
    P[:ε_i] = lininterp(P[:t_s],P[:ε], t_UT)
end

"""
    SFP_interp!(P, t_UT)

    Interpolate pump and furnace data to match ultrasonic data pulses
    #Arguments
    * SFP : dictionary containing processed data from the pump/furnace logging program
    * exp_info : details of experiment
"""
function SFP_interp!(SFP, t_UT)
    P[:T_lw_i] = lininterp(SFP[:t],SFP[:T_lw], t_UT)
    P[:T_uw_i] = lininterp(SFP[:t],P[:T_uw], t_UT)
    P[:T_sf_i] = lininterp(SFP[:t],P[:T_sf], t_UT)
end
