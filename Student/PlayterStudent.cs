//============================================================================
// PlayterSim.cs
// Simulation of a Playter Doll.
//============================================================================
using System;

public partial class PlayterSim : Simulator
{
    /* … existing field declarations stay the same … */

    //-----------------------------------------------------------------------
    // RHSFuncPlayter: Evaluates the right sides of the differential equations
    //-----------------------------------------------------------------------
    private void RHSFuncPlayter(double[] xx, double t, double[] ff)
    {
        /* ----------- 1.  UNPACK STATE (unchanged) ------------------------ */
        omegaX  = xx[0];  omegaY  = xx[1];  omegaZ  = xx[2];
        omegaFL = xx[3];  omegaFR = xx[4];
        vx      = xx[5];  vy      = xx[6];  vz      = xx[7];

        q0 = xx[8];  q1 = xx[9];  q2 = xx[10];  q3 = xx[11];
        thetaL = xx[12];  thetaR = xx[13];
        xG = xx[14];  yG = xx[15];  zG = xx[16];

        /* ----------- 2.  BODY ANGULAR ACCEL (Playter eq. 14) ------------- */
        double rho2 = rho * rho;
        Vex AngVel  = new Vex(omegaX, omegaY, omegaZ);
        Vex AngMo   = new Vex(rho2 * omegaX,
                              rho2 * gammaY * omegaY,
                              rho2 * gammaZ * omegaZ);
        Vex AngVelCrossAngMo = Vex.Cross(AngVel, AngMo);
        ff[0] = -AngVelCrossAngMo.x / rho2;
        ff[1] = -AngVelCrossAngMo.y / (rho2 * gammaY);
        ff[2] = -AngVelCrossAngMo.z / (rho2 * gammaZ);

        /* ----------- 3.  QUATERNION KINEMATICS (eq. 5) ------------------- */
        ff[8]  = .5*(-q1*omegaX - q2*omegaY - q3*omegaZ);
        ff[9]  = .5*( q0*omegaX - q3*omegaY + q2*omegaZ);
        ff[10] = .5*( q3*omegaX + q0*omegaY - q1*omegaZ);
        ff[11] = .5*(-q2*omegaX + q1*omegaY + q0*omegaZ);

        /* =====================  NEW CODE START  ========================= */
        // 4. Shoulder spring‑damper torques  (eq. 39)
        double TL = -mA * L * L * (k * thetaL + c * omegaFL);   // *** NEW ***
        double TR = -mA * L * L * (k * thetaR + c * omegaFR);   // *** NEW ***
        double Iarm = mA * L * L;                               // *** NEW ***

        // 5. Arm angular accelerations
        ff[3] = TL / Iarm;   // ω̇FL                              // *** NEW ***
        ff[4] = TR / Iarm;   // ω̇FR                              // *** NEW ***

        // 6. Couple the arm torques back to body (simple proj. onto b̂z)
        double armCoupleZ = (TL - TR) * cosPhi;                  // *** NEW ***
        ff[2] += armCoupleZ / (rho2 * gammaZ);                   // *** NEW ***
        ff[1] += (TL + TR) * (-sinPhi) / (rho2 * gammaY);        // *** NEW ***

        // 7. Kinematic derivatives for arms and CG
        ff[12] = omegaFL;   // θ̇L                                // *** NEW ***
        ff[13] = omegaFR;   // θ̇R                                // *** NEW ***
        ff[14] = vx;        // ẋG                                // *** NEW ***
        ff[15] = vy;        // ẏG                                // *** NEW ***
        ff[16] = vz;        // żG                                // *** NEW ***

        // 8. No translational accelerations (add gravity later)
        ff[5] = ff[6] = ff[7] = 0.0;                             // *** NEW ***
        /* =====================  NEW CODE END    ========================= */

        /* ----------- 9.  DEBUG VALUES (unchanged) ----------------------- */
        SetDebugVal(0, omegaX);
        SetDebugVal(1, omegaY);
        SetDebugVal(2, omegaZ);
    }
} // end class

