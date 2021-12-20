#pragma once

//
// Reference from https://en.wikipedia.org/wiki/Sellmeier_equation
// and https://refractiveindex.info/
//

namespace png {
	enum TransparentMaterialType{ 
		Sapphire
		,BK7
		,HighVariance
	};
	double inline refractiveIndex(TransparentMaterialType materialType, double wavelength_nm) {
		double B1, B2, B3, C1, C2, C3;
		switch (materialType)
		{
		case Sapphire:
			B1 = 1.5039759;
			B2 = 0.55069141;
			B3 = 6.5927379;
			C1 = 5.48041129e-3;
			C2 = 1.47994281e-2;
			C3 = 402.89514;
			break;
		case BK7:
			B1 = 1.03961212;
			B2 = 0.231792344;
			B3 = 1.01046945;
			C1 = 0.00600069867;
			C2 = 0.0200179144;
			C3 = 103.560653;
			break;
		case HighVariance:
			B1 = -2.0;
			B2 = B3 = 0.0;
			C1 = C2 = C3 = 1.0;
			break;
		default:
			break;
		}

		if (wavelength_nm >= 0) {
			double x = 0.001 * wavelength_nm; // = wavelength (micro metres)
			double n_sq = 1 
				+ B1 * x * x / (x * x - C1) 
				+ B2 * x * x / (x * x - C2) 
				+ B3 * x * x / (x * x - C3);
			return sqrt(n_sq);
		}
		else {
			return 1;
		}
	}
}