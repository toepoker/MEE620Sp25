//============================================================================
// PlayterDoll.cs
// Code for runing the Godot scene of rotating playter doll
//============================================================================
using Godot;
using System;

public partial class PlayterDoll : Node3D
{
	// parameters
	double ma = 0.107;     // arm mass ratio
	double rho = 2.21;     // dimensionless radius of gyration
	double gammay = 0.091; // ratio I_Gy/I_Gx
	double gammaz = 1.05;  // ratio I_Gz/I_Gx
	double h = 1.56;       // dimensionless height of shoulder above CG
	double L = 1.65;       // dimensionless length to arm mass
	double Lrod = 2.5;    // length of arm rod (for visual)
	double phi = 0.0;      // shoulder angle

	float cgHeight;

	enum RunMode{
		Config,    // Configuration mode (choose IC)
		Sim,
		Pause
	}
	RunMode runMode;

	// simulation
	PlayterSim sim;
	double time;

	// model
	PDollModel model;

	// Camera Stuff
	CamRig cam;
	float longitudeDeg;
	float latitudeDeg;
	float camDist;
	float camFOV;
	Vector3 camTg;       // coords of camera target

	double theta;   // somersault angle for testing

	Label timeLabel;
	CheckButton cButtonArmMode;
	LineEdit[] icOmega;
	SpinBox[] sbIcOmega;
	Button simButton;
	Button resetButton;

	bool testing;
	int nTestVals;
	CheckButton testButton;
	GridContainer testGrid;
	Label[] testVals;

	int dispCtr;
	int dispTHold;

	//------------------------------------------------------------------------
	// _Ready: Called once when the node enters the scene tree for the first 
	//         time.
	//------------------------------------------------------------------------
	public override void _Ready()
	{
		GD.Print("PlayterDoll");

		runMode = RunMode.Config;

		cgHeight = 6.0f;
		sim = new PlayterSim();
		time = 0.0;
		SetParamsDoll();

		// set up doll model
		model = GetNode<PDollModel>("PDollModel");
		model.Position = new Vector3(0.0f, cgHeight, 0.0f);
		model.Initialize((float)h, (float)L, (float)Lrod, 3.0f /*hUb*/, 
		3.75f /*hLb*/, 0.5f /*thk*/, (float)phi);

		// Set up the camera rig
		longitudeDeg = 30.0f;  //30
		latitudeDeg = 15.0f;
		camDist = 12.0f;

		camTg = new Vector3(0.0f, cgHeight, 0.0f);
		cam = GetNode<CamRig>("CamRig");
		cam.LongitudeDeg = longitudeDeg;
		cam.LatitudeDeg = latitudeDeg;
		cam.Distance = camDist;
		cam.Target = camTg;
		cam.Position = new Vector3(0.0f, cgHeight, 0.0f);

		theta = 0.0;

		dispCtr = 0;
		dispTHold = 3;

		SetupUI();

		/*
		LinSysEq sys;
		int ii,jj;
        Random rnd = new Random();
        sys = new LinSysEq(8);
        for(ii=0; ii<8; ++ii){
            for(jj=0; jj<8; ++jj){
                sys.SetA(ii, jj, rnd.NextDouble());
            }
			sys.SetB(ii, rnd.NextDouble());
        }

		sys.SolveGauss();
		GD.Print("Check = " + sys.Check());
		/**/
	}

	//------------------------------------------------------------------------
	// _Process: Called every frame. 'delta' is the elapsed time since the 
	//           previous frame.
	//------------------------------------------------------------------------
	public override void _Process(double delta)
	{
		int i;

		//model.SetSimpleRotation((float)theta);
		float q0 = (float)sim.Q0;
		float q1 = (float)sim.Q1;
		float q2 = (float)sim.Q2;
		float q3 = (float)sim.Q3;
		float thL = (float)sim.ThetaL;
		float thR = (float)sim.ThetaR;
		float xG = (float)sim.XG;
		float yG = (float)sim.YG;
		float zG = (float)sim.ZG;

		model.Update(q0, q1, q2, q3, thL, thR, xG, yG, zG);

		timeLabel.Text = time.ToString("0.0");

		if(dispCtr > dispTHold){
			for(i=0;i<16;++i){
				testVals[i].Text = sim.GetDebugVal(i).ToString("0.0000");
			}
			dispCtr = 0;
		}
		else{
			++dispCtr;
		}
	}

	//------------------------------------------------------------------------
	// SetParamsDoll
	//------------------------------------------------------------------------
	private void SetParamsDoll()
	{
		ma = 0.107;     // arm mass ratio
		rho = 2.21;     // dimensionless radius of gyration
		gammay = 0.091; // ratio I_Gy/I_Gx
		gammaz = 1.05;  // ratio I_Gz/I_Gx
		h = 1.56;       // dimensionless height of shoulder above CG
		L = 1.65;       // dimensionless length to arm mass
		Lrod = 2.5;    // length of arm rod (for visual)
		phi = 0.0;      // shoulder angle

		sim.ArmMass = ma;
		sim.RadiusOfGyration = rho;
		sim.GammaY = gammay;
		sim.GammaZ = gammaz;
		sim.ShoulderHeight = h;
		sim.ArmLength = L;
		sim.Phi = phi;
	}

	//------------------------------------------------------------------------
	// SetParamsHuman
	//------------------------------------------------------------------------
	private void SetParamsHuman()
	{
		ma = 0.107;     // arm mass ratio
		rho = 2.21;     // dimensionless radius of gyration
		gammay = 0.091; // ratio I_Gy/I_Gx
		gammaz = 1.05;  // ratio I_Gz/I_Gx
		h = 1.56;       // dimensionless height of shoulder above CG
		L = 1.65;       // dimensionless length to arm mass
		Lrod = 2.5;    // length of arm rod (for visual)
		phi = 0.0;      // shoulder angle

		sim.ArmMass = ma;
		sim.RadiusOfGyration = rho;
		sim.GammaY = gammay;
		sim.GammaZ = gammaz;
		sim.ShoulderHeight = h;
		sim.ArmLength = L;
		sim.Phi = phi;
	}

	//------------------------------------------------------------------------
	// _PhysicsProcess:
	//------------------------------------------------------------------------
    public override void _PhysicsProcess(double delta)
    {
        base._PhysicsProcess(delta);

		if(runMode == RunMode.Sim){
			sim.Step(time, delta);
			time += delta;
		}
    }

	//------------------------------------------------------------------------
	// SetupUI
	//------------------------------------------------------------------------
	private void SetupUI()
	{
		int i;

		MarginContainer mcTL = GetNode<MarginContainer>("Control/MgCTL/pContain/MC");

		VBoxContainer vbox = new VBoxContainer();
		mcTL.AddChild(vbox);

		Label titleLabel = new Label();
		titleLabel.Text = "Playter Doll Simulator";
		vbox.AddChild(titleLabel);
		vbox.AddChild(new HSeparator());

		GridContainer grid = new GridContainer();
		grid.Columns = 2;
		vbox.AddChild(grid);

		Label tlab = new Label();
		tlab.Text = "Time: ";
		timeLabel = new Label();
		timeLabel.Text = "0.0";
		grid.AddChild(tlab);
		grid.AddChild(timeLabel);
		grid.AddChild(new HSeparator());
		grid.AddChild(new HSeparator());

		Label[] icLabel;
		icLabel = new Label[3];
		icOmega = new LineEdit[3];
		sbIcOmega = new SpinBox[3];
		for(i=0; i<3; ++i){
			icLabel[i] = new Label();
			//icOmega[i] = new LineEdit();
			//icOmega[i].TooltipText = "Enter initial condition.";
			sbIcOmega[i] = new SpinBox();
			grid.AddChild(icLabel[i]);
			//grid.AddChild(icOmega[i]);
			grid.AddChild(sbIcOmega[i]);
			//icOmega[i] = sbIcOmega[i].GetLineEdit();
			sbIcOmega[i].Step = 0.01f;
			sbIcOmega[i].MaxValue = 2.0f;
			sbIcOmega[i].ValueChanged += OnICValueChanged;
			//sbIcOmega[i].Suffix = "rad/s";
		}
		icLabel[0].Text = "IC: OmegaX: ";
		icLabel[1].Text = "IC: OmegaY: ";
		icLabel[2].Text = "IC: OmegaZ: ";

		sbIcOmega[0].Value = 1.0f;
		sbIcOmega[1].Value = 0.0f;
		sbIcOmega[2].Value = 0.0f;

		// icOmega[0].Text = "1.0";
		// icOmega[1].Text = "0.0";
		// icOmega[2].Text = "0.0";

		grid.AddChild(new HSeparator());
		grid.AddChild(new HSeparator());

		cButtonArmMode = new CheckButton();
		cButtonArmMode.Text = "Arm Free";
		cButtonArmMode.ButtonPressed = true;
		cButtonArmMode.Toggled += OnArmModeToggle;
		vbox.AddChild(cButtonArmMode);
		vbox.AddChild(new HSeparator());

		simButton = new Button();
		simButton.Text = "Simulate";
		simButton.Pressed += OnButtonSimPause;
		vbox.AddChild(simButton);

		resetButton = new Button();
		resetButton.Text = "Reset";
		resetButton.Pressed += OnButtonReset;
		resetButton.Disabled = true;
		vbox.AddChild(resetButton);


		// Test grid in the upper right margin container
		nTestVals = 16;
		MarginContainer mcTR = GetNode<MarginContainer>("Control/MgCTR/pContain/MC");

		VBoxContainer vboxTest = new VBoxContainer();
		mcTR.AddChild(vboxTest);
		
		testButton = new CheckButton();
		testButton.Text = "Test Data";
		testButton.Toggled += OnTestingToggle;
		vboxTest.AddChild(testButton);
		vboxTest.AddChild(new HSeparator());

		testGrid = new GridContainer();
		testGrid.Columns = 2;
		vboxTest.AddChild(testGrid);

		Label[] testLabels = new Label[nTestVals];
		testVals = new Label[nTestVals];
		for (i = 0; i<nTestVals; ++i){
			testLabels[i] = new Label();
			testLabels[i].Text = "TestVal_" + i.ToString() + ": ";
			testVals[i] = new Label();
			testVals[i].Text = "0.0000";
			testGrid.AddChild(testLabels[i]);
			testGrid.AddChild(testVals[i]);
		} 
		
		testGrid.Hide();
	}

	//------------------------------------------------------------------------
	// OnICValueChanged:
	//------------------------------------------------------------------------
	private void OnICValueChanged(double val)
	{
		double omX = sbIcOmega[0].Value;
		double omY = sbIcOmega[1].Value;
		double omZ = sbIcOmega[2].Value;
		
		//######### Hmmmmm. Do I need to do anything here???

		//GD.Print("Omega = " + omX + ", " + omY + ", " + omZ);
	}

	//------------------------------------------------------------------------
	// OnArmModeToggle:
	//------------------------------------------------------------------------
	private void OnArmModeToggle(bool freeSelected)
	{
		if(freeSelected){
			//GD.Print("Free");
			cButtonArmMode.Text = "Arm Free";
			sim.SetArmFree();
		}
		else{
			//GD.Print("Prescribed");
			cButtonArmMode.Text = "Arm Prescribed";
			sim.SetArmPrescribed();
		}
	}

	//------------------------------------------------------------------------
	// OnTestingToggle
	//------------------------------------------------------------------------
	private void OnTestingToggle(bool testingSelected)
	{
		if(testingSelected){
			//GD.Print("Testing");
			//testButton.Text = "Testing";
			testGrid.Show();
		}
		else{
			//GD.Print("Not Testing");
			//testButton.Text = "Not Testing";
			testGrid.Hide();
		}
	}

	//------------------------------------------------------------------------
	// OnButtonSimPause
	//------------------------------------------------------------------------
	private void OnButtonSimPause()
	{
		if(runMode == RunMode.Config){ // switch from Config to Sim
			simButton.Text = "Pause";
			resetButton.Disabled = false;
			sbIcOmega[0].Editable = false;
			sbIcOmega[1].Editable = false;
			sbIcOmega[2].Editable = false;

			// send IC to sim
			double omX = sbIcOmega[0].Value;
			double omY = sbIcOmega[1].Value;
			double omZ = sbIcOmega[2].Value;
			sim.SetSpinIC(omX, omY, omZ);

			runMode = RunMode.Sim;
		}
		else if(runMode == RunMode.Pause){ // switch from Pause to Sim
			simButton.Text = "Pause";

			runMode = RunMode.Sim;
		}
		else if(runMode == RunMode.Sim){  // switch from Sim to Pause
			simButton.Text = "Simulate";

			runMode = RunMode.Pause;
		}
	}

	//------------------------------------------------------------------------
	// OnButtonReset
	//------------------------------------------------------------------------
	private void OnButtonReset()
	{
		simButton.Text = "Simulate";
		resetButton.Disabled = true;
		sbIcOmega[0].Editable = true;
		sbIcOmega[1].Editable = true;
		sbIcOmega[2].Editable = true;

		timeLabel.Text = "0.0";
		time = 0.0;

		sim.SetSpinIC(0.0); // reorient doll to upright

		runMode = RunMode.Config;
	}
}
