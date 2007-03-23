#ifndef ATOMDATA_H
#define ATOMDATA_H
#include <vector>

namespace chemistry {
class Atomdata {
	public:
		Atomdata();
		Atomdata(int i);
		double rm() const;
		double density(double r) const;
		int number() const;
	private:
		std::vector<double> m_r;
		std::vector<double> m_density;
		double m_rm;
		int m_atomicnumber;
};
};




#endif
