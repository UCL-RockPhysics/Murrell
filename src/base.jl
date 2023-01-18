"""
    M_read(P)

    Open and read contents of .tdms file
    #Arguments
    - fid : file path string
    return P a dictionary containing:
    * :t_s : time in seconds, LabView time base (see https://tinyurl.com/56yhfrw4)
    * :F_kN : raw record of external load
    * :Ua_mm : actuator displacement (mm)
    * :U1_mm : LVDT1 displacement (mm)
    * :U2_mm : LVDT2 displacement (mm)
    * :Pc2_MPa : gas pressure (MPa)
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
        P[:Pc2_MPa] = tdmsIN.groups["Numeric"]["PC 1400 MPa"].data
    elseif tdmsIN.groups["Numeric"]["TimeStamp"].data[1] >= 3.7120499126540936e9 && tdmsIN.groups["Numeric"]["TimeStamp"].data[1] < 3.756725042945276e9# This is here since the LabView software was updated and some of the variable names changed
        A = zeros(6,length(tdmsIN.groups["Numeric"]["TimeStamp"].data))
        A[1,:] = tdmsIN.groups["Numeric"]["TimeStamp"].data
        A[2,:] = tdmsIN.groups["Numeric"]["Load"].data
        A[3,:] = tdmsIN.groups["Numeric"]["Displacement"].data
        A[4,:] = tdmsIN.groups["Numeric"]["LVDT1"].data
        A[5,:] = tdmsIN.groups["Numeric"]["LVDT2"].data
        A[6,:] = tdmsIN.groups["Numeric"]["PC1400"].data
        A = unique(A, dims=2) # There is an occassional bug where data is written to the file twice
        P[:t_s] = A[1,:]
        P[:F_kN] = A[2,:]
        P[:Ua_mm] = A[3,:]
        P[:U1_mm] = A[4,:]
        P[:U2_mm] = A[5,:]
        P[:Pc2_MPa] = A[6,:]
    else
        A = zeros(7,length(tdmsIN.groups["Numeric"]["TimeStamp"].data))
        A[1,:] = tdmsIN.groups["Numeric"]["TimeStamp"].data
        A[2,:] = tdmsIN.groups["Numeric"]["Load"].data
        A[3,:] = tdmsIN.groups["Numeric"]["Displacement"].data
        A[4,:] = tdmsIN.groups["Numeric"]["LVDT1"].data
        A[5,:] = tdmsIN.groups["Numeric"]["LVDT2"].data
        A[6,:] = tdmsIN.groups["Numeric"]["PC1400"].data
        A[7,:] = tdmsIN.groups["Numeric"]["Internal Load"]
        A = unique(A, dims=2) # There is an occassional bug where data is written to the file twice
        P[:t_s] = A[1,:]
        P[:F_kN] = A[2,:]
        P[:Ua_mm] = A[3,:]
        P[:U1_mm] = A[4,:]
        P[:U2_mm] = A[5,:]
        P[:Pc2_MPa] = A[6,:]
        P[:ILC_V] = A[7,:]
    end
    return P
end
# """
# function t_conv!(P)
#     convert LabView timestamp to Unix timestamp
#     #Arguments
#     * P, dictionary containg key :t_s extracted from .tdms file
#     #Returns
#     * :t_s : overwritten value of converted time, s
# """
# function t_conv!(P)
#     P[:t_s] = unix2datetime(datetime2unix(DateTime(1904,1,1,0,0,0)) .+(P[:t_s]))
# end

"""
    M_reduce!(P, exp_info; stresscorr=true, unloadcorr=true)
    Reduce Murrell data to calculate stress, strain etc.
    #Arguments
    * P : raw data extracted from data file
    * exp_info : dictionary containing information about test conditions, must include indices of hit point and sample information, as follows
        exp_info = Dict()
        exp_info[:I]=[1, 2, 3] # Indices at HP, start of unloading and end of unloading
        exp_info[:IE]=[1, 2, 3, 4] # Indices to estimate moduli
        exp_info[:P_exp] = 200 # Pressure of experiment
        exp_info[:T] = 200 # Temperature of experiment
        exp_info[:εr] = 1e-5 # Strain rate of test
        exp_info[:K_mm_kN] = K_D2 # Machine compliance [mm kN^-1]
        exp_info[:L_mm] = 22.00 # Sample length
        exp_info[:d_mm] = 10.00 # Diameter of sample [m]
    * exp_info : dictionary containing experimental parameters
    * stresscorr : correct computed differential stress assuming constant volume deformation
    * unloadcorr : routine to estimate stress on sample during decompression as a result of seal friction
    returns P : dictionary:
        * :σ3_MPa : confining pressure, MPa
        * :t_s_c : time relative to start of experiment, s
        * :U_mm_c : LVDT recorded displacement, averaged between LVDT1 and LVDT2
        * :U_mm_fc : displacement corrected for load and (optionally) confining pressure
        * :F_kN_c : force relative to hitpoint, kN
        * :ε : natural strain, -
        * :σ_MPa_j : differential stress corrected for jacket strength, MPa
"""

function M_reduce!(P, exp_info; stresscorr=true, unloadcorr=true)
    dt = P[:t_s][2].-P[:t_s][1]
    sf = Int(ceil(1/dt))
    I1 = exp_info[:I][1]
    I2 = exp_info[:I][2]
    I3 = exp_info[:IE][3]
    I4 = exp_info[:I][3]
    P[:σ3_MPa] = movingaverage((P[:Pc2_MPa]),sf)
    P[:Pc_corr] = zeros(length(P[:σ3_MPa]))
    for i = 1:length(P[:Pc_corr])
        if P[:σ3_MPa][i] < 250       # Correction to seal friction for confining pressure change relative to hitpoint
            P[:Pc_corr][i] = 0.0219*P[:σ3_MPa][i] +0.8763
        else
            P[:Pc_corr][i] = 8.569e-5*P[:σ3_MPa][i]^2 -0.017*P[:σ3_MPa][i] +5.3616868619086855
        end
    end
    P[:Pc_corr] .-= P[:Pc_corr][I1]
    P[:t_s_c] = P[:t_s].-P[:t_s][I1] #
    P[:U_mm_c] = movingaverage((P[:U1_mm].+P[:U2_mm]),sf)./2
    if unloadcorr == true
        P[:F_kN_c] = zeros(length(P[:F_kN]),1)
        F_smooth = movingaverage(P[:F_kN] .-P[:F_kN][I1],sf)
        ## Force corrections accounting for loading sense and unloading behaviour
        P[:F_kN_c][1:I2] = F_smooth[1:I2].-P[:Pc_corr][1:I2] # Correct force measurements where load is positive sense
        P[:F_kN_c][I2+1:I3] =   P[:F_kN_c][I2] .-
                                (P[:Pc_corr][I2].-
                                P[:Pc_corr][I2+1:I3]) # Correct for piston locked by friction
        P[:F_kN_c][I3+1:I4] =   P[:F_kN_c][I3] .+
                                (F_smooth[I3+1:I4].-F_smooth[I3]) .-
                                (P[:Pc_corr][I3].-P[:Pc_corr][I3+1:I4]) # Correct unload for offset resulting from friction until out of contact
        P[:F_kN_c][I4+1:end] =  P[:F_kN_c][I4] .+
                                (F_smooth[I4+1:end].-F_smooth[I4]) .-
                                (P[:Pc_corr][I4].-P[:Pc_corr][I4+1:end]) # Now differential load is given by seal friction, correct for this
        ## Displacement corrections based on load sense
        P[:U_mm_fc] = zeros(length(P[:U_mm_c]))
        P[:U_mm_fc][1:I1-1]     =   (P[:σ3_MPa][1:I1-1] ./45.8e3).*exp_info[:L_mm]
        P[:U_mm_fc][I1:I2]      =   P[:U_mm_fc][I1-1].+
                                (P[:U_mm_c][I1:I2].-P[:U_mm_c][I1-1]).-
                                (P[:F_kN_c][I1:I2]*exp_info[:K_mm_kN]).-
                                ((P[:σ3_MPa][I1-1].-P[:σ3_MPa][I1:I2])*1.26e-3) # Correct displacement for machine stiffness and confining pressure applying a correction of 1.26 µm/MPa
        P[:U_mm_fc][I2+1:I3]    =   P[:U_mm_fc][I2].-
                                ((P[:σ3_MPa][I2+1:I3].-P[:σ3_MPa][I2]).*1.26e-3) # Correct displacement for machine stiffness and confining pressure applying a correction of 1.26 µm/MPa
        P[:U_mm_fc][I3+1:end]   =   P[:U_mm_fc][I3].+
                                (P[:U_mm_c][I3+1:end].-P[:U_mm_c][I3]).-
                                ((P[:F_kN_c][I3+1:end].-P[:F_kN_c][I3])*exp_info[:K_mm_kN]).-
                                ((P[:σ3_MPa][I3+1:end].-P[:σ3_MPa][I3])*1.26e-3)
        P[:U_mm_c] .-= P[:U_mm_c][I1]
        P[:U_mm_fc] .-= P[:U_mm_fc][1]
        P[:U_mm_fc][1:I1-1] .= P[:U_mm_fc][I1]
    else
        P[:F_kN_c] = P[:F_kN].-P[:Pc_corr] .-P[:F_kN][I1]
        P[:U_mm_c] .-= P[:U_mm_c][I1]
        P[:U_mm_fc] = P[:U_mm_c].-P[:F_kN_c].*exp_info[:K_mm_kN]
        P[:U_mm_fc][1:I1-1] .= 0
    end
    ## Strain computation
    P[:ε] = -log.(1 .-(P[:U_mm_fc].-P[:U_mm_fc][I1])./exp_info[:L_mm]) # Compute natural strain
    P[:Jr] = JR!(P, exp_info) # Get force resulting from jacket
    P[:F_kN_j] = P[:F_kN_c] .-P[:Jr] # Correct force for jacket rheology
    if stresscorr == true
        P[:σ_MPa_j] = P[:F_kN_c]./(π.*(0.5e-3*exp_info[:d_mm].*(1 .+P[:ε])).^2).*1e-3 # Apply surface area correction to stress
    else
        P[:σ_MPa_j] = P[:F_kN_c]./(π.*(0.5e-3*exp_info[:d_mm]).^2).*1e-3 # Do not apply surface area correction
    end
end


"""
    SFP_read(fid)

    Read and return data from program SFP_furnace_control_V2.vi

    #Arguments
    - 'fid::string' : path to file
"""
function SFP_read(fid::String)
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
    ε̇j = exp_info[:εr] # Jacket strain rate
    Jr = 2*10^((log10(ε̇j*exp(EaJ/(8.3145*(exp_info[:T]+278)))))/n1J-n2J) # Copper flow stress at experiment conditions
    Ja = π*(5.25e-3^2-5e-3^2)  # Jacket area
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
    P[:F_kN_i] = lininterp(P[:t_s],vec(P[:F_kN_j]), t_UT)
    P[:U_mm_i] = lininterp(P[:t_s],vec(P[:U_mm_c]), t_UT)
    P[:σ_MPa_i] = lininterp(P[:t_s],vec(P[:σ_MPa_j]), t_UT)
    P[:σ3_MPa_i] = lininterp(P[:t_s],vec(P[:Pc2_MPa]), t_UT)
    P[:ε_i] = lininterp(P[:t_s],vec(P[:ε]), t_UT)
end

"""
    SFP_interp!(P, t_UT)

    Interpolate pump and furnace data to match ultrasonic data pulses
    #Arguments
    * SFP : dictionary containing processed data from the pump/furnace logging program
    * exp_info : details of experiment
"""
function SFP_interp!(SFP, t_UT)
    SFP[:T_lw_i] = lininterp(float.(SFP[:t]),float.(SFP[:T_lw]), t_UT)
    SFP[:T_uw_i] = lininterp(float.(SFP[:t]),float.(SFP[:T_uw]), t_UT)
    SFP[:T_sf_i] = lininterp(float.(SFP[:t]),float.(SFP[:T_sf]), t_UT)
end
