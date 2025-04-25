//============================================================================
// PlayterStudent.cs – student‑editable dynamics for the Playter Doll
//-----------------------------------------------------------------------------
// 2025‑04‑25  – step‑2 update
//   • added momentum‑consistent translational accelerations (PDF §5–6)
//   • leaves placeholders (TODO) for future matrix form (PDF §8)
//============================================================================
using System;

public partial class PlayterSim : Simulator
{
    // -------------------------------------------------------------------------
    // Extra working objects ----------------------------------------------------
    // -------------------------------------------------------------------------
    LinSysEq sys;          // 8×8 linear equation solver (for future use)
    double[,] Amat;        // mass / inertia matrix  (8 × 8)
    double[]  Bmat;        // RHS vector             (8 × 1)

    //-------------------------------------------------------------------------
    // StudentInit – allocate arrays once before the simulation starts
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
        /* ---- 1. UNPACK STATE ----------------------------------------------*/
        omegaX  = xx[0];  omegaY  = xx[1];  omegaZ  = xx[2];
        omegaFL = xx[3];  omegaFR = xx[4];
        vx      = xx[5];  vy      = xx[6];  vz      = xx[7];

        q0 = xx[8]; q1 = xx[9]; q2 = xx[10]; q3 = xx[11];
        thetaL = xx[12]; thetaR = xx[13];
        xG = xx[14]; yG = xx[15]; zG = xx[16];

        /* ---- 2. BODY ANGULAR ACCELERATION (PDF §4 eq. 14) -----------------*/
        double rho2 = rho * rho;
        Vex angVel = new Vex(omegaX, omegaY, omegaZ);
        Vex angMo  = new Vex(rho2 * omegaX,
                             rho2 * gammaY * omegaY,
                             rho2 * gammaZ * omegaZ);
        Vex crossH = Vex.Cross(angVel, angMo);

        ff[0] = -crossH.x / rho2;                  // ω̇x
        ff[1] = -crossH.y / (rho2 * gammaY);       // ω̇y
        ff[2] = -crossH.z / (rho2 * gammaZ);       // ω̇z

        /* ---- 3. QUATERNION KINEMATICS (PDF §3 eq. 5) ----------------------*/
        ff[8]  = 0.5 * (-q1 * omegaX - q2 * omegaY - q3 * omegaZ);
        ff[9]  = 0.5 * ( q0 * omegaX - q3 * omegaY + q2 * omegaZ);
        ff[10] = 0.5 * ( q3 * omegaX + q0 * omegaY - q1 * omegaZ);
        ff[11] = 0.5 * (-q2 * omegaX + q1 * omegaY + q0 * omegaZ);

        // Quaternion‑norm debug
        double qNorm = Math.Sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
        SetDebugVal(3, qNorm);

        /* =================== 4. SHOULDER TORQUES (PDF §7 eq. 39) ===========*/
        double TL   = -mA * L * L * (k * thetaL + c * omegaFL);
        double TR   = -mA * L * L * (k * thetaR + c * omegaFR);
        double Iarm =  mA * L * L;   // arm inertia about hinge

        /* ---- 5. ARM ANGULAR ACCELERATIONS --------------------------------*/
        ff[3] = TL / Iarm;           // ω̇FL
        ff[4] = TR / Iarm;           // ω̇FR

        /* ---- 6. REACTION COUPLES ON BODY ----------------------------------*/
        double armCoupleZ = (TL - TR) *  cosPhi;
        double armCoupleY = (TL + TR) * (-sinPhi);
        ff[2] += armCoupleZ / (rho2 * gammaZ);
        ff[1] += armCoupleY / (rho2 * gammaY);

        /* ---- 7. KINEMATIC DERIVATIVES (θ̇, ẋ, ẏ, ż) ----------------------*/
        ff[12] = omegaFL;            // θ̇L
        ff[13] = omegaFR;            // θ̇R
        ff[14] = vx;                 // ẋG
        ff[15] = vy;                 // ẏG
        ff[16] = vz;                 // żG

        /* ===================================================================
         * 8. TRANSLATIONAL ACCELERATIONS – momentum consistency
         *    PDF §5–6  (Eq. 37 simplified)  – no external forces.
         *    We compute arm‑mass accelerations relative to G, then choose aG
         *    so that Σ F_int = 0  ⇒  total linear momentum is conserved.
         * -------------------------------------------------------------------*/
        // Constant helpers
        double Mtot = 1.0 + 2.0 * mA;     // nondim total mass (body + 2 arms)

        // --- geometry (positions relative to body COM) --------------------
        Vex rSLG = new Vex( 1.0,  h, 0.0);   // left shoulder rel G
        Vex rSRG = new Vex(-1.0,  h, 0.0);   // right shoulder rel G

        Vex rFLS = new Vex( L * Math.Cos(thetaL),
                            L * Math.Sin(thetaL) *  cosPhi,
                            L * Math.Sin(thetaL) *  sinPhi);
        Vex rFRS = new Vex(-L * Math.Cos(thetaR),
                           -L * Math.Sin(thetaR) * cosPhi,
                           -L * Math.Sin(thetaR) * sinPhi);

        Vex rFLG = rFLS + rSLG;   // left arm mass rel G
        Vex rFRG = rFRS + rSRG;   // right arm mass rel G

        // --- angular kinematics ------------------------------------------
        Vex omegaNB  = new Vex(omegaX, omegaY, omegaZ);
        Vex alphaNB  = new Vex(ff[0],  ff[1],  ff[2]);

        Vex Sz = new Vex(0.0, -sinPhi, cosPhi);  // ŝz axis
        Vex omegaFLB = omegaFL * Sz;
        Vex omegaFRB = omegaFR * Sz;
        Vex alphaFLB = ff[3]  * Sz;
        Vex alphaFRB = ff[4]  * Sz;

        // --- arm accelerations (α × r + ω × (ω × r) for each rotation) ----
        // Left arm: acceleration of shoulder rel G, then add arm swing terms.
        Vex aSL = Vex.Cross(alphaNB, rSLG) + Vex.Cross(omegaNB, Vex.Cross(omegaNB, rSLG));
        Vex aFLSrel = Vex.Cross(alphaNB, rFLS) + Vex.Cross(omegaNB, Vex.Cross(omegaNB, rFLS))
                      + Vex.Cross(alphaFLB, rFLS) + Vex.Cross(omegaFLB, Vex.Cross(omegaFLB, rFLS))
                      + 2.0 * Vex.Cross(omegaNB, Vex.Cross(omegaFLB, rFLS));
        Vex aFL = aSL + aFLSrel;   // absolute acc of left arm mass

        // Right arm
        Vex aSR = Vex.Cross(alphaNB, rSRG) + Vex.Cross(omegaNB, Vex.Cross(omegaNB, rSRG));
        Vex aFRSrel = Vex.Cross(alphaNB, rFRS) + Vex.Cross(omegaNB, Vex.Cross(omegaNB, rFRS))
                      + Vex.Cross(alphaFRB, rFRS) + Vex.Cross(omegaFRB, Vex.Cross(omegaFRB, rFRS))
                      + 2.0 * Vex.Cross(omegaNB, Vex.Cross(omegaFRB, rFRS));
        Vex aFR = aSR + aFRSrel;

        // --- choose aG so that Σ F_int = 0 (no external force) ------------
        Vex aGvec = (-mA / Mtot) * (aFL + aFR);
        ff[5] = aGvec.x;
        ff[6] = aGvec.y;
        ff[7] = aGvec.z;

        /* ---- 9. DEBUG OUTPUT ---------------------------------------------*/
        SetDebugVal(0, omegaX);
        SetDebugVal(1, omegaY);
        SetDebugVal(2, omegaZ);
        SetDebugVal(4, ff[0]);       // ω̇x
        SetDebugVal(5, ff[5]);       // aₓ (check ~0 if COM stationary)
    }
} // end partial class PlayterSim
