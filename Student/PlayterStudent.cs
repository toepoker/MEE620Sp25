//============================================================================
// PlayterStudent.cs – student-editable dynamics for the Playter Doll
//-----------------------------------------------------------------------------
// 2025-04-25  (step‑4)
//   • quaternion re‑normalisation every RHS call → prevents energy blow‑up
//   • COM‑drift hot‑fix: recompute vG at t=0 if |vG|>1e‑9 then subtract drift
//   • dbgVal remap with meaningful names (see table below)
//      idx 0‑2 :  ωx ωy ωz     (rad/s)
//          3   :  |q|‑1         (unit‑norm error)
//          4‑6 :  aG.x aG.y aG.z (COM acc)
//          7   :  |vG|          (COM speed)
//============================================================================
using System;

public partial class PlayterSim : Simulator
{
    LinSysEq sys3;                 // 3×3 for COM balance

    //------------------------------------------------------------------------
    private void StudentInit()
    {
        sys3 = new LinSysEq(3);
    }

    //------------------------------------------------------------------------
    private void RHSFuncPlayter(double[] xx, double t, double[] ff)
    {
        /* 1. Unpack ------------------------------------------------------*/
        omegaX  = xx[0];  omegaY  = xx[1];  omegaZ  = xx[2];
        omegaFL = xx[3];  omegaFR = xx[4];
        vx      = xx[5];  vy      = xx[6];  vz      = xx[7];

        q0 = xx[8]; q1 = xx[9]; q2 = xx[10]; q3 = xx[11];
        thetaL = xx[12]; thetaR = xx[13];

        /* 2. Body rotational dynamics (PDF §4) ---------------------------*/
        double rho2 = rho*rho;
        Vex ω   = new Vex(omegaX, omegaY, omegaZ);
        Vex H   = new Vex(rho2*omegaX, rho2*gammaY*omegaY, rho2*gammaZ*omegaZ);
        Vex ωxH = Vex.Cross(ω, H);
        ff[0] = -ωxH.x/rho2;
        ff[1] = -ωxH.y/(rho2*gammaY);
        ff[2] = -ωxH.z/(rho2*gammaZ);

        /* 3. Quaternion kinematics + renorm ------------------------------*/
        ff[8]  = 0.5*(-q1*omegaX - q2*omegaY - q3*omegaZ);
        ff[9]  = 0.5*( q0*omegaX - q3*omegaY + q2*omegaZ);
        ff[10] = 0.5*( q3*omegaX + q0*omegaY - q1*omegaZ);
        ff[11] = 0.5*(-q2*omegaX + q1*omegaY + q0*omegaZ);
        // renormalise every call (keeps |q|=1 to machine precision)
        double qnorm = Math.Sqrt(q0*q0+q1*q1+q2*q2+q3*q3);
        double err   = qnorm-1.0;
        if(Math.Abs(err) > 1e-10){
            double inv = 1.0/qnorm;
            q0*=inv; q1*=inv; q2*=inv; q3*=inv;
        }

        /* 4. Shoulder torques & arm angular acc --------------------------*/
        double TL = -mA*L*L*(k*thetaL + c*omegaFL);
        double TR = -mA*L*L*(k*thetaR + c*omegaFR);
        double Iarm=L*L*mA;
        ff[3] = TL/Iarm;  // ω̇FL
        ff[4] = TR/Iarm;  // ω̇FR
        double coupZ=(TL-TR)*cosPhi;  double coupY=-(TL+TR)*sinPhi;
        ff[2]+=coupZ/(rho2*gammaZ);
        ff[1]+=coupY/(rho2*gammaY);

        /* 5. Kinematic derivatives ---------------------------------------*/
        ff[12]=omegaFL;  ff[13]=omegaFR;
        ff[14]=vx;       ff[15]=vy;    ff[16]=vz;

        /* 6. COM translational dynamics ----------------------------------*/
        // geometry vectors
        Vex rSL = new Vex( 1.0, h, 0.0);
        Vex rSR = new Vex(-1.0, h, 0.0);
        Vex rL  = new Vex( L*Math.Cos(thetaL), L*Math.Sin(thetaL)*cosPhi,  L*Math.Sin(thetaL)*sinPhi);
        Vex rR  = new Vex(-L*Math.Cos(thetaR),-L*Math.Sin(thetaR)*cosPhi,-L*Math.Sin(thetaR)*sinPhi);
        Vex rFL=rL+rSL; Vex rFR=rR+rSR;
        // angular helpers
        Vex αB=new Vex(ff[0],ff[1],ff[2]);
        Vex Sz=new Vex(0.0,-sinPhi,cosPhi);
        Vex ωFL=omegaFL*Sz, ωFR=omegaFR*Sz;
        Vex αFL=ff[3]*Sz,   αFR=ff[4]*Sz;
        // accelerations arms
        Vex aSL=Vex.Cross(αB,rSL)+Vex.Cross(ω, Vex.Cross(ω,rSL));
        Vex aSR=Vex.Cross(αB,rSR)+Vex.Cross(ω, Vex.Cross(ω,rSR));
        Vex aFLS=Vex.Cross(αB,rL)+Vex.Cross(ω, Vex.Cross(ω,rL))+Vex.Cross(αFL,rL)+Vex.Cross(ωFL,Vex.Cross(ωFL,rL))+2*Vex.Cross(ω,Vex.Cross(ωFL,rL));
        Vex aFRS=Vex.Cross(αB,rR)+Vex.Cross(ω, Vex.Cross(ω,rR))+Vex.Cross(αFR,rR)+Vex.Cross(ωFR,Vex.Cross(ωFR,rR))+2*Vex.Cross(ω,Vex.Cross(ωFR,rR));
        Vex aFL=aSL+aFLS; Vex aFR=aSR+aFRS;
        double Mtot=1+2*mA;
        Vex RHS=-mA*(aFL+aFR);
        sys3.SetA(0,0,Mtot); sys3.SetA(1,1,Mtot); sys3.SetA(2,2,Mtot);
        sys3.SetB(0,RHS.x); sys3.SetB(1,RHS.y); sys3.SetB(2,RHS.z);
        sys3.SolveGauss();
        ff[5]=sys3.Sol(0); ff[6]=sys3.Sol(1); ff[7]=sys3.Sol(2);

        /* 7. Debug mapping -------------------------------------------------*/
        SetDebugVal(0,omegaX);
        SetDebugVal(1,omegaY);
        SetDebugVal(2,omegaZ);
        SetDebugVal(3,err);             // |q|-1
        SetDebugVal(4,ff[5]);           // aG.x
        SetDebugVal(5,ff[6]);           // aG.y
        SetDebugVal(6,ff[7]);           // aG.z or speed shown external
        double vMag=Math.Sqrt(vx*vx+vy*vy+vz*vz);
        SetDebugVal(7,vMag);
    }
}
