#pragma once
#ifdef _HAS_STD_BYTE
#undef _HAS_STD_BYTE
#endif
#define _HAS_STD_BYTE 0


#include "BGAL/BaseShape/Point.h"
#include "BGAL/BaseShape/Polygon.h"
#include "BGAL/BaseShape/Triangle.h"
#include "BGAL/BaseShape/Line.h"
#include "BGAL/Model/ManifoldModel.h"
#include "BGAL/Model/Model_Iterator.h"
#include "BGAL/Tessellation3D/Tessellation3D.h"
#include "BGAL/Optimization/LBFGS/LBFGS.h"

#include <string>

namespace BGAL
{
	class _CVT3D
	{
		std::string outpath = "../../data/LBFGSOUT/";
	public:
		_CVT3D(const _ManifoldModel& model);
		_CVT3D(const _ManifoldModel& model, std::function<double(_Point3& p)>& rho, _LBFGS::_Parameter para);
		void calculate_(int site_num, char* modelNamee, char* pointsName = nullptr);

		void set_outpath(const std::string& path)
		{
			outpath = path;
		}

		const std::vector<_Point3>& get_sites() const
		{
			return _sites;
		}
		const _Restricted_Tessellation3D& get_RVD() const
		{
			return _RVD;
		}
		const _Restricted_Tessellation3D& get_RVD2() const
		{
			return _RVD2;
		}
	public:
		const _ManifoldModel& _model;
		_Restricted_Tessellation3D _RVD;
		_Restricted_Tessellation3D _RVD2;
		std::vector<_Point3> _sites{};
		std::function<double(_Point3& p)> _rho;
		_LBFGS::_Parameter _para;
	};
} // namespace BGAL
