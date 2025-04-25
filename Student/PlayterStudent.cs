//============================================================================
// PlayterStudent.cs  –  student‑editable dynamics for the Playter Doll project
//-----------------------------------------------------------------------------
// Final consolidated version (2025‑04‑25)
//   • Exact rotational dynamics (PDF §4) and shoulder spring‑damper (§7)
//   • Quaternion renormalisation each step (keeps |q| ≈ 1)
//   • Linear‑momentum balance solved with LinSysEq (3×3) every call so the
//     centre‑of‑mass (COM) velocity remains zero when no external forces act.
//   • dbgVal ↔ HUD mapping:
//         0 ωx   1 ωy   2 ωz   3 |q|–1
//         4 aG.x 5 aG.y 6 aG.z 7 |vG|
//============================================================================
using System;

public partial class PlayterSim : Simulator
{
    // 3×3 system for COM acceleration each RHS call
    LinSysEq sys3;

    //------------------------------------------------------------------------
    // One‑time allocation
    //------------------------------------------------------------------------
    private void StudentInit()
    {
        sys3 = new LinSysEq(3);
    }

    //------------------------------------------------------------------------
    // Right‑hand side  ẋ = f(x,t)
    //------------------------------------------------------------------------
    private void RHSFuncPlayter(double[] xState, double t, double[] fOut)
    {
        /* ===== 1. UNPACK STATE ==========================================*/
        //  ── generalized speeds
        omegaX  = xState[0];  omegaY  = xState[1];  omegaZ  = xState[2];
        omegaFL = xState[3];  omegaFR = xState[4];
        vx      = xState[5];  vy      = xState[6];  vz      = xState[7];
        //  ── generalized coordinates
        q0 = xState[8];   q1 = xState[9];   q2 = xState[10];  q3 = xState[11];
        thetaL = xState[12];  thetaR = xState[13];

        /* ===== 2. ROTATIONAL DYNAMICS (PDF §4 Eq 14) ==================*/
        double ρ2 = rho * rho;
        Vex ω = new Vex(omegaX, omegaY, omegaZ);
        Vex H = new Vex(ρ2*omegaX, ρ2*gammaY*omegaY, ρ2*gammaZ*omegaZ);
        Vex ω×H = Vex.Cross(ω, H);
        fOut[0] = -ω×H.x / ρ2;               // ω̇x
        fOut[1] = -ω×H.y / (ρ2*gammaY);      // ω̇y
        fOut[2] = -ω×H.z / (ρ2*gammaZ);      // ω̇z

        /* ===== 3. QUATERNION KINEMATICS (PDF §3 Eq 5) =================*/
        fOut[8]  = 0.5*(-q1*omegaX - q2*omegaY - q3*omegaZ);
        fOut[9]  = 0.5*( q0*omegaX - q3*omegaY + q2*omegaZ);
        fOut[10] = 0.5*( q3*omegaX + q0*omegaY - q1*omegaZ);
        fOut[11] = 0.5*(-q2*omegaX + q1*omegaY + q0*omegaZ);
        // ── renormalise quaternion every step
        double qNorm = Math.Sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
        double qErr  = qNorm - 1.0;
        if(Math.Abs(qErr) > 1e-10)
        {
            double inv = 1.0 / qNorm;
            q0*=inv; q1*=inv; q2*=inv; q3*=inv;
        }

        /* ===== 4. SHOULDER TORQUES & ARM ANGULAR ACC ===================*/
        double Iarm = mA * L * L;
        double TL = -Iarm * (k*thetaL + c*omegaFL);
        double TR = -Iarm * (k*thetaR + c*omegaFR);
        fOut[3] = TL / Iarm;     // ω̇FL
        fOut[4] = TR / Iarm;     // ω̇FR
        // reaction couples on body
        fOut[2] +=  (TL-TR)*cosPhi / (ρ2*gammaZ);
        fOut[1] += -(TL+TR)*sinPhi / (ρ2*gammaY);

        /* ===== 5. KINEMATIC DERIVATIVES ================================*/
        fOut[12] = omegaFL;   // θ̇L
        fOut[13] = omegaFR;   // θ̇R
        fOut[14] = vx;        // ẋG
        fOut[15] = vy;        // ẏG
        fOut[16] = vz;        // żG

        /* ===== 6. TRANSLATIONAL DYNAMICS (COM balance) =================*/
        // ---- geometry (positions in body frame)
        Vex rSL = new Vex( 1.0, h, 0.0);  // left shoulder → G
        Vex rSR = new Vex(-1.0, h, 0.0);  // right shoulder → G
        Vex rL  = new Vex( L*Math.Cos(thetaL),  L*Math.Sin(thetaL)*cosPhi,  L*Math.Sin(thetaL)*sinPhi);
        Vex rR  = new Vex(-L*Math.Cos(thetaR), -L*Math.Sin(thetaR)*cosPhi, -L*Math.Sin(thetaR)*sinPhi);
        Vex rFL = rL + rSL;   // arm masses rel G
        Vex rFR = rR + rSR;

        // ---- angular helpers
        Vex αB  = new Vex(fOut[0], fOut[1], fOut[2]);
        Vex Sz  = new Vex(0.0, -sinPhi, cosPhi);
        Vex ωFL = omegaFL * Sz;   Vex ωFR = omegaFR * Sz;
        Vex αFL = fOut[3]  * Sz;   Vex αFR = fOut[4]  * Sz;

        // ---- absolute accelerations of arm masses
        Vex aSL = Vex.Cross(αB, rSL) + Vex.Cross(ω, Vex.Cross(ω, rSL));
        Vex aSR = Vex.Cross(αB, rSR) + Vex.Cross(ω, Vex.Cross(ω, rSR));
        Vex aFL = aSL + Vex.Cross(αB, rL) + Vex.Cross(ω, Vex.Cross(ω, rL))
                       + Vex.Cross(αFL, rL) + Vex.Cross(ωFL, Vex.Cross(ωFL, rL))
                       + 2.0*Vex.Cross(ω, Vex.Cross(ωFL, rL));
        Vex aFR = aSR + Vex.Cross(αB, rR) + Vex.Cross(ω, Vex.Cross(ω, rR))
                       + Vex.Cross(αFR, rR) + Vex.Cross(ωFR, Vex.Cross(ωFR, rR))
                       + 2.0*Vex.Cross(ω, Vex.Cross(ωFR, rR));

        // ---- solve Mtot aG + mA(aFL+aFR) = 0  (LinSysEq 3×3)
        double Mtot = 1.0 + 2.0*mA;
        sys3.ZeroA();                 // IMPORTANT: clear any old coeffs
        sys3.SetA(0,0,Mtot); sys3.SetA(1,1,Mtot); sys3.SetA(2,2,Mtot);
        Vex RHS = -mA*(aFL + aFR);
        sys3.SetB(0, RHS.x); sys3.SetB(1, RHS.y); sys3.SetB(2, RHS.z);
        sys3.SolveGauss();
        fOut[5] = sys3.Sol(0);   // aG.x
        fOut[6] = sys3.Sol(1);   // aG.y
        fOut[7] = sys3.Sol(2);   // aG.z

        /* ===== 7. DEBUG OUTPUT (dbgVal array) ===========================*/
        SetDebugVal(0, omegaX);      // ωx
        SetDebugVal(1, omegaY);      // ωy
        SetDebugVal(2, omegaZ);      // ωz
        SetDebugVal(3, qErr);        // |q|–1
        SetDebugVal(4, fOut[5]);     // aG.x
        SetDebugVal(5, fOut[6]);     // aG.y
        SetDebugVal(6, fOut[7]);     // aG.z
        double vMag = Math.Sqrt(vx*vx + vy*vy + vz*vz);
        SetDebugVal(7, vMag);        // |vG|
    }
}

