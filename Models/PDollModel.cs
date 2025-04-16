//============================================================================
// PDollModel.cs
// Model of the Playter Doll.
//============================================================================
using Godot;
using System;

public partial class PDollModel : Node3D
{
	bool initialized;

	Node3D rootNode;
	Node3D vertShift;
	PDUpperBody ubody;
	PdLowerBody lbody;
	Node3D shoulderFrame;
	Node3D afL;            // left arm frame
	Node3D afR;            // right arm frame

	Vector3 armRotL;
	Vector3 armRotR;
	Vector3 rG;

	Vector3 eulerAngles;

	//------------------------------------------------------------------------
	// _Ready: Called once when the node enters the scene tree for the first 
	//         time.
	//------------------------------------------------------------------------
	public override void _Ready()
	{
		initialized = false;

		rootNode = GetNode<Node3D>("RootNode");
		vertShift = GetNode<Node3D>("RootNode/VertShift");
		ubody = GetNode<PDUpperBody>("RootNode/VertShift/PDUpperBody");
		lbody = GetNode<PdLowerBody>("RootNode/VertShift/PDLowerBody");
		shoulderFrame = vertShift.GetNode<Node3D>("ShoulderFrame");
		afL = shoulderFrame.GetNode<Node3D>("ArmFrameL");
		afR = shoulderFrame.GetNode<Node3D>("ArmFrameR");

		eulerAngles = new Vector3(0.0f, 0.0f, 0.0f);
		rG = new Vector3();
		armRotL = new Vector3();
		armRotR = new Vector3();
	}

	//------------------------------------------------------------------------
	// update  Updates the orientation of the model and the angles of the
	//         arms.
	//------------------------------------------------------------------------
	public void Update(float q0, float q1, float q2, float q3,
		float thetaL, float thetaR, float xG, float yG, float zG)
	{
		this.Quaternion = new Quaternion(q1, q2, q3, q0);

		rG.X = xG;   rG.Y = yG;  rG.Z = zG;
		rootNode.Position = rG;

		armRotL.Z = thetaL;
		afL.Rotation = armRotL;

		armRotR.Z = thetaR;
		afR.Rotation = armRotR;
	}

	//------------------------------------------------------------------------
	// Initialize:
	//------------------------------------------------------------------------
	public void Initialize(float h, float L, float Lrod, float hUB, float hLB,
		float thk, float shAngle)
	{
		float shoulderOffset = 0.5f*thk;

		ubody.GenMesh(1.0f, hUB, 0.5f*thk);
		vertShift.Position = new Vector3(0.0f, h+shoulderOffset, 0.0f);
		lbody.Position = new Vector3(0.0f, -hUB, 0.0f);
		lbody.GenMesh(1.0f, hLB, 0.5f*thk, 0.5f);
		
		shoulderFrame.Rotation = new Vector3(shAngle, 0.0f, 0.0f);
		shoulderFrame.Position = new Vector3(0.0f, -shoulderOffset, 0.0f);

		MeshInstance3D ShCylL = afL.GetNode<MeshInstance3D>("CylL");
		ShCylL.Position = new Vector3(0.0f, 0.0f, 0.0f);
		ShCylL.Rotation = new Vector3(-0.5f*MathF.PI, 0.0f, 0.0f);
		CylinderMesh shMeshL = (CylinderMesh)ShCylL.Mesh;
		shMeshL.TopRadius = shoulderOffset;
		shMeshL.BottomRadius = shoulderOffset;
		shMeshL.Height = 1.2f*thk;

		MeshInstance3D ShCylR = afR.GetNode<MeshInstance3D>("CylR");
		ShCylR.Position = new Vector3(0.0f, 0.0f, 0.0f);
		ShCylR.Rotation = new Vector3(0.5f*MathF.PI, 0.0f, 0.0f);
		CylinderMesh shMeshR = (CylinderMesh)ShCylR.Mesh;
		shMeshR.TopRadius = shoulderOffset;
		shMeshR.BottomRadius = shoulderOffset;
		shMeshR.Height = 1.2f*thk;

		afL.Position = new Vector3(1.0f, 0.0f, 0.0f);
		MeshInstance3D aRodL = afL.GetNode<MeshInstance3D>("Rod");
		aRodL.Position = new Vector3(0.5f*Lrod, 0.0f, 0.0f);
		aRodL.Rotation = new Vector3(0.0f, 0.0f, MathF.PI*0.5f);
		CylinderMesh aRodMeshL = (CylinderMesh)aRodL.Mesh;
		aRodMeshL.TopRadius = 0.1f;
		aRodMeshL.BottomRadius = 0.1f;
		aRodMeshL.Height = Lrod;
		//GD.Print("Lrod = " + Lrod);

		MeshInstance3D aMassL = afL.GetNode<MeshInstance3D>("Mass");
		aMassL.Position = new Vector3(L, 0.0f, 0.0f);
		SphereMesh aSphMeshL = (SphereMesh)aMassL.Mesh;
		aSphMeshL.Radius = 0.4f;
		aSphMeshL.Height = 0.8f;

		afR.Position = new Vector3(-1.0f, 0.0f, 0.0f);
		MeshInstance3D aRodR = afR.GetNode<MeshInstance3D>("Rod");
		aRodR.Position = new Vector3(-0.5f*Lrod, 0.0f, 0.0f);
		aRodR.Rotation = new Vector3(0.0f, 0.0f, MathF.PI*0.5f);
		CylinderMesh aRodMeshR = (CylinderMesh)aRodR.Mesh;
		aRodMeshR.TopRadius = 0.1f;
		aRodMeshR.BottomRadius = 0.1f;
		aRodMeshR.Height = Lrod;

		MeshInstance3D aMassR = afR.GetNode<MeshInstance3D>("Mass");
		aMassR.Position = new Vector3(-L, 0.0f, 0.0f);
		SphereMesh aSphMeshR = (SphereMesh)aMassR.Mesh;
		aSphMeshR.Radius = 0.4f;
		aSphMeshR.Height = 0.8f;
	}


	//------------------------------------------------------------------------
	// UpdateSimpleRotation
	//------------------------------------------------------------------------
	public void SetSimpleRotation(float th)
	{
		eulerAngles.X = th;
		eulerAngles.Y = 0.0f;
		eulerAngles.Z = 0.0f;

		rootNode.Rotation = eulerAngles;
	}


	// Called every frame. 'delta' is the elapsed time since the previous frame.
	// public override void _Process(double delta)
	// {
	// }
}
