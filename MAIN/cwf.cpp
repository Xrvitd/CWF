#include <iostream>
#include <fstream>

#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>

#include <BGAL/CVTLike/CVT.h>


void CWF3D(std::string file, std::string pointsFile)
{
	std::string filepath = "../../data/";
	std::string modelname = file;

	// .obj to .off
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	igl::readOBJ(modelname, V, F);

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
	BGAL::_CVT3D cvt(model, rho, para);
	int num = 0;

	cvt.set_outpath("./");
	std::string filename = std::filesystem::path(modelname).filename().string();
	cvt.calculate_(num, (char*)filename.c_str(), (char*)pointsFile.c_str());
}


int main(int argc, char* argv[])
{
	CWF3D(argv[1], argv[2]);

	return 0;
}