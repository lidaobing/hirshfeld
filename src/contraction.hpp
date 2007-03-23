//: contraction.h
// all units is a.u.

#ifndef CONTRACTION_H
#define CONTRACTION_H

#include <vector>
#include <iosfwd>




namespace chemistry {

        class Atom;

	class Contraction {
		public:
			enum {D5 = -2, SP, S, P, D, F};
			enum {S_0};
			enum {P_x, P_y, P_z};
			enum {D_xx, D_yy, D_zz, D_xy, D_xz, D_yz};
			//D5_0  = d_{3z^2-r^2}
			//D5_p1 = d_{xz}
			//D5_n1 = d_{yz}
			//D5_p2 = d_{x^2-y^2}
			//D5_n2 = d_{xy}
			enum {D5_0, D5_p1, D5_n1, D5_p2, D5_n2};
			enum {D5_z2, D5_xz, D5_yz, D5_x2y2, D5_xy};


			Contraction(Atom & atom, int majortype, int minortype);
			void add(double alpha, double coef);
			void calcN();

			double calc(double x, double y, double z) const;
			friend std::ostream& operator<<(std::ostream& os, const Contraction& con);
		private:

			struct Prim {
				double alpha;
				double N;
				double coef;
			};

			Atom & m_atom;

			int m_majortype;
			int m_minortype;

			double & m_x;
			double & m_y;
			double & m_z;

			std::vector<Prim> m_prims;

			double m_N;

			double int2prim(int index1, int index2) const;
	};
	std::ostream& operator<<(std::ostream& os, const Contraction& con);

};

#endif
