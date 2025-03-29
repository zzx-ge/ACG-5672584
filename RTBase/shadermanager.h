class ShaderManager
{
public:
	bool TM = false; //ToneMapping
	bool BVH = false; //Bounding Volume Hierarchy
	bool MIS = false;	//Multiple Importance Sampling
	bool EM = false; //Environment Mapping
	bool VPL = false; //Virtual Point Light
	bool LT = false; //Light Tracing
	int Filter = 0; // 0 : boxfilter, 1 : gaussianfilter, 2 : MitchellNetravaliFilter
	bool denoiser = false;
	ShaderManager() {}
};