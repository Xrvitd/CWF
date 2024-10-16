#include <iostream>
#include <fstream>
#include <functional>
#include <io.h>
#include <random>
#include <omp.h>

#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>
#include <BGAL/Optimization/LinearSystem/LinearSystem.h>
#include <BGAL/Optimization/ALGLIB/optimization.h>
#include <BGAL/Optimization/LBFGS/LBFGS.h>
#include <BGAL/Integral/Integral.h>
#include <BGAL/Model/ManifoldModel.h>
#include <BGAL/Model/Model_Iterator.h>

#include <BGAL/Tessellation3D/Tessellation3D.h>

#include <BGAL/CVTLike/CPD.h>
#include <BGAL/CVTLike/CVT.h>




void CWF3DTest()
{
	// your can try: Nums = 600      1000   1000          2000 	   
	//               File = mobius1  block  block_smooth  bunny	 


	int Nums = 600;
	std::string file = "mobius1";
	std::cout << "Now file: " << file << std::endl;
	
	std::cout << Nums << "   " << file << "   \n";
	std::string filepath = "../../data/";
	std::string modelname = file;

	// .obj to .off
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	igl::readOBJ(filepath + modelname + ".obj", V, F);

	igl::writeOFF("Temp.off", V, F);
	igl::writeOBJ("Temp.obj", V, F);

	BGAL::_ManifoldModel model("Temp.obj");

	std::function<double(BGAL::_Point3& p)> rho = [](BGAL::_Point3& p)
		{
			return 1;
		};

	BGAL::_LBFGS::_Parameter para;
	para.is_show = true;
	para.epsilon = 1e-30;
	para.max_iteration = 50;
	BGAL::_CVT3D cvt(model, rho, para);
	int num = Nums;
	cvt.calculate_(num, (char*)modelname.c_str());
	
	


}




int alltest()
{
	

	std::cout << "====================CWF3DTest" << std::endl;
	CWF3DTest();
	std::cout << "successful!" << std::endl;
	return 0;
}

int main()
{
	alltest();
	return 0;
}