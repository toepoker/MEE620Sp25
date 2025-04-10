//============================================================================
// PlayterSim.cs
// Simulation of a Playter Doll.
//
// This is where students write their code.
//============================================================================
using System;

public partial class PlayterSim : Simulator
{
    // // Parameters
    // double mA;     // dimensionless arm mass
    // double rho;    // dimless radius if gyration for moment of inertia, I_Gx
    // double gammaY; // ratio I_Gy/I_Gx
    // double gammaZ; // ratio I_Gz/I_Gx
    // double h;      // dimless vertical (b_y) distance of shoulder from body CG
    // double L;      // dimless distance of arm mass from shoulder
    // double k;      // dimless torsional stiffness of shoulder spring
    // double c;      // dimless torsional damping coeff in shoulder damper
    // double phi;    // angle of arm swing plane relative to vertical
    // double cosPhi;
    // double sinPhi;

    // // generalized speeds
    // double omegaX;
    // double omegaY;
    // double omegaZ;
    // double omegaFL;
    // double omegaFR;
    // double vx;
    // double vy;
    // double vz;

    // // generalized coordinates
    // double q0;      // quaternion coords
    // double q1;
    // double q2;
    // double q3;
    // double thetaL;  // left arm angle
    // double thetaR;  // right arm angle
    // double xG;      // coordinates of body's center of mass
    // double yG;
    // double zG;

    // Some extra stuff a student might want (feel free to define more.)
    LinSysEq sys;    // linear algebraic equation solver
    double[,] Amat;  // some arrays that might be handy 
    double[] Bmat;



    //------------------------------------------------------------------------
    // StudentInit: Student might want to initialize things before simulation
    //              begins
    //------------------------------------------------------------------------
    private void StudentInit()
    {
        sys = new LinSysEq(8);
        Amat = new double[8,8];  // an 8 by 8 array of doubles
        Bmat = new double[8];    // an 8 by 1 array of doubles
        
    }

    //------------------------------------------------------------------------
    // RHSFuncPlayter:  Evaluates the right sides of the differential
    //                   equations for the Playter Doll
    //------------------------------------------------------------------------
    private void RHSFuncPlayter(double[] xx, double t, double[] ff)
    {
        omegaX  = xx[0];
        omegaY  = xx[1];
        omegaZ  = xx[2];
        omegaFL = xx[3];
        omegaFR = xx[4];
        vx      = xx[5];
        vy      = xx[6];
        vz      = xx[7];

        q0      = xx[8];
        q1      = xx[9];
        q2      = xx[10];
        q3      = xx[11];
        thetaL  = xx[12];
        thetaR  = xx[13];
        xG      = xx[14];
        yG      = xx[15];
        zG      = xx[16];

        ff[8] = .5*(-q1*omegaX - q2*omegaY - q3*omegaZ);
        ff[9] = .5*(q0*omegaX - q3*omegaY + q2*omegaZ);
        ff[10] = .5*(q3*omegaX + q0*omegaY - q1*omegaZ);
        ff[11] = .5*(-q2*omegaX + q1*omegaY + q0*omegaZ);
        

        // COMMENT THESE OUT OR REMOVE WHEN READY
        ff[0] = ff[1] = ff[2] = 0.0;   // derivs of body angular velocities set to zero
        ff[3] = ff[4] = 0.0;           // derivs of arm angular velocities
        ff[5] = ff[6] = ff[7] = 0.0;   // derivs of cener of mass velocities
        ff[8] = ff[9] = ff[10] = ff[11] = 0.0;  // derivs of quaternion coords
        ff[12] = ff[13] = 0.0;         // derivs of arm angles
        ff[14] = ff[15] = ff[16] = 0.0;  // derivs of CG coordinates
    }

} // end class
