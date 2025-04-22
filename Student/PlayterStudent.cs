//============================================================================
// PlayterSim.cs
// Simulation of a Playter Doll.
//============================================================================
using System;

public partial class PlayterSim : Simulator
{
    // -------------------------------------------------------------------------
    // Parameter fields are declared in another partial file (PlayterParams.cs)
    // -------------------------------------------------------------------------
    // double mA, rho, gammaY, gammaZ, h, L, k, c, phi, cosPhi, sinPhi;
    // double omegaX, omegaY, omegaZ, omegaFL, omegaFR, vx, vy, vz;
    // double q0, q1, q2, q3, thetaL, thetaR, xG, yG, zG;

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
        // (add any other one‑time setup here)
    }

    //-------------------------------------------------------------------------
    // RHSFuncPlayter: evaluates the right‑hand side ẋ = f(x, t)
    //-------------------------------------------------------------------------
    private void RHSFuncPlayter(double[] xx, double t, double[] ff)
    {
        /* ---- 1. UNPACK STATE ------------------------------------------------*/
        omegaX  = xx[0];  omegaY  = xx[1];  omegaZ  = xx[2];
        omegaFL = xx[3];  omegaFR = xx[4];
        vx      = xx[5];  vy      = xx[6];  vz      = xx[7];

        q0 = xx[8]; q1 = xx[9]; q2 = xx[10]; q3 = xx[11];
        thetaL = xx[12]; thetaR = xx[13];
        xG = xx[14]; yG = xx[15]; zG = xx[16];

        /* ---- 2. BODY ANGULAR ACCELERATION (eq. 14) -------------------------*/
        double rho2 = rho * rho;
        Vex angVel = new Vex(omegaX, omegaY, omegaZ);
        Vex angMo  = new Vex(rho2 * omegaX,
                             rho2 * gammaY * omegaY,
                             rho2 * gammaZ * omegaZ);
        Vex cross  = Vex.Cross(angVel, angMo);

        ff[0] = -cross.x / rho2;                 // ω̇X
        ff[1] = -cross.y / (rho2 * gammaY);      // ω̇Y
        ff[2] = -cross.z / (rho2 * gammaZ);      // ω̇Z

        /* ---- 3. QUATERNION KINEMATICS (eq. 5) -----------------------------*/
        ff[8]  = 0.5 * (-q1 * omegaX - q2 * omegaY - q3 * omegaZ);
        ff[9]  = 0.5 * ( q0 * omegaX - q3 * omegaY + q2 * omegaZ);
        ff[10] = 0.5 * ( q3 * omegaX + q0 * omegaY - q1 * omegaZ);
        ff[11] = 0.5 * (-q2 * omegaX + q1 * omegaY + q0 * omegaZ);

        /* =================== 4. SHOULDER TORQUES (eq. 39) ===================*/
        double TL   = -mA * L * L * (k * thetaL + c * omegaFL);
        double TR   = -mA * L * L * (k * thetaR + c * omegaFR);
        double Iarm =  mA * L * L;   // moment of inertia of point mass arm

        /* ---- 5. ARM ANGULAR ACCELERATIONS ---------------------------------*/
        ff[3] = TL / Iarm;           // ω̇FL
        ff[4] = TR / Iarm;           // ω̇FR

        /* ---- 6. REACTION COUPLE ON BODY -----------------------------------*/
        double armCoupleZ = (TL - TR) * cosPhi;
        ff[2] += armCoupleZ / (rho2 * gammaZ);
        ff[1] += (TL + TR) * (-sinPhi) / (rho2 * gammaY);

        /* ---- 7. KINEMATIC DERIVATIVES (θ̇, ẋ, ẏ, ż) -----------------------*/
        ff[12] = omegaFL;            // θ̇L
        ff[13] = omegaFR;            // θ̇R
        ff[14] = vx;                 // ẋG
        ff[15] = vy;                 // ẏG
        ff[16] = vz;                 // żG

        /* ---- 8. TRANSLATIONAL ACCELS (add gravity later) ------------------*/
        ff[5] = ff[6] = ff[7] = 0.0;

        /* ---- 9. DEBUG OUTPUT ----------------------------------------------*/
        SetDebugVal(0, omegaX);
        SetDebugVal(1, omegaY);
        SetDebugVal(2, omegaZ);
    }
} // end class PlayterSim


