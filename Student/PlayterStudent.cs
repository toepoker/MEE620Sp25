//============================================================================
// PlayterStudent.cs
// Simulation of a Playter Doll.
// (student-editable dynamics file)
//============================================================================
using System;

public partial class PlayterSim : Simulator
{
    // -------------------------------------------------------------------------
    // Extra working objects ---------------------------------------------------
    // -------------------------------------------------------------------------
    LinSysEq sys;          // linear equation solver (size 8 × 8)
    double[,] Amat;        // mass / inertia matrix  (8 × 8)
    double[]  Bmat;        // RHS vector             (8 × 1)

    //-------------------------------------------------------------------------
    // StudentInit: allocate arrays once before the simulation starts
    //-------------------------------------------------------------------------
    private void StudentInit()
    {
        sys  = new LinSysEq(8);
        Amat = new double[8, 8];
        Bmat = new double[8];
        // (add any other one-time setup here, e.g. pre-compute constants)
    }

    //-------------------------------------------------------------------------
    // RHSFuncPlayter: evaluates the right-hand side ẋ = f(x, t)
    //-------------------------------------------------------------------------
    private void RHSFuncPlayter(double[] xx, double t, double[] ff)
    {
        /* ---- 1. UNPACK STATE ----------------------------------------------*/
        omegaX  = xx[0];  omegaY  = xx[1];  omegaZ  = xx[2];
        omegaFL = xx[3];  omegaFR = xx[4];
        vx      = xx[5];  vy      = xx[6];  vz      = xx[7];

        q0 = xx[8]; q1 = xx[9]; q2 = xx[10]; q3 = xx[11];
        thetaL = xx[12]; thetaR = xx[13];
        xG = xx[14]; yG = xx[15]; zG = xx[16];

        /* ---- 2. BODY ANGULAR ACCELERATION (PDF §4 eq. 14) -----------------*/
        double rho2 = rho * rho;           // ρ²
        Vex angVel = new Vex(omegaX, omegaY, omegaZ);
        Vex angMo  = new Vex(rho2 * omegaX,                // Lx  = ρ² ωx
                            rho2 * gammaY * omegaY,       // Ly  = ρ² γy ωy
                             rho2 * gammaZ * omegaZ);      // Lz  = ρ² γz ωz
        Vex cross  = Vex.Cross(angVel, angMo);              // ω × H

        ff[0] = -cross.x / rho2;                 // ω̇x
        ff[1] = -cross.y / (rho2 * gammaY);      // ω̇y
        ff[2] = -cross.z / (rho2 * gammaZ);      // ω̇z

        /* ---- 3. QUATERNION KINEMATICS (PDF §3 eq. 5) ----------------------*/
        ff[8]  = 0.5 * (-q1 * omegaX - q2 * omegaY - q3 * omegaZ);
        ff[9]  = 0.5 * ( q0 * omegaX - q3 * omegaY + q2 * omegaZ);
        ff[10] = 0.5 * ( q3 * omegaX + q0 * omegaY - q1 * omegaZ);
        ff[11] = 0.5 * (-q2 * omegaX + q1 * omegaY + q0 * omegaZ);

        // Quaternion-norm check (debug)
        double qNorm = Math.Sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
        SetDebugVal(3, qNorm);        // should stay ~1.0

        /* =================== 4. SHOULDER TORQUES (PDF §7 eq. 39) ============*/
        double TL   = -mA * L * L * (k * thetaL + c * omegaFL);
        double TR   = -mA * L * L * (k * thetaR + c * omegaFR);
        double Iarm =  mA * L * L;   // point-mass arm inertia about hinge

        /* ---- 5. ARM ANGULAR ACCELERATIONS --------------------------------*/
        ff[3] = TL / Iarm;           // ω̇FL
        ff[4] = TR / Iarm;           // ω̇FR

        /* ---- 6. REACTION COUPLES ON BODY ----------------------------------*/
        // Couples arise because equal-and-opposite shoulder torques act on body.
        double armCoupleZ = (TL - TR) * cosPhi;          // about b̂z
        double armCoupleY = -(TL + TR) * sinPhi;         // about b̂y

        ff[2] += armCoupleZ / (rho2 * gammaZ);           // add to ω̇z
        ff[1] += armCoupleY / (rho2 * gammaY);           // add to ω̇y

        /* ---- 7. KINEMATIC DERIVATIVES (θ̇, ẋ, ẏ, ż) ----------------------*/
        ff[12] = omegaFL;            // θ̇L
        ff[13] = omegaFR;            // θ̇R
        ff[14] = vx;                 // ẋG
        ff[15] = vy;                 // ẏG
        ff[16] = vz;                 // żG

        /* ---- 8. TRANSLATIONAL ACCELS (free-floating, no gravity) ----------*/
        ff[5] = 0.0;                 // aₓ (N.frame)
        ff[6] = 0.0;                 // a_y
        ff[7] = 0.0;                 // a_z

        /* ---- 9. DEBUG OUTPUT ---------------------------------------------*/
        SetDebugVal(0, omegaX);
        SetDebugVal(1, omegaY);
        SetDebugVal(2, omegaZ);
        SetDebugVal(4, ff[0]);       // ω̇x (diagnostic)
    }
} // end partial class PlayterSim


