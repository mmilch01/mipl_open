/*******************************************************************************
Copyright (c) 2016, 2017, 2018
Authored by: Mikhail Milchenko, mmilchenko@wustl.edu
This software is free for non-commercial use by acadimic, government, and non-profit/not-for-profit institutions. If this code is modified and/or redistributed, this note must appear in the beginning of each of re-distributed source code files. A commercial version of the software is available, and is licensed through the Office of Technology Management at Washington University School of Medicine. For more information regarding the commercial license, please email your request to otm@dom.wustl.edu.

BY DOWNLOADING THIS SOFTWARE YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS PROVIDED 'AS IS', WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT. IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
**********************************************************************************/

#include "mlib3math.h"
#include "mlib3volume.h"
#include <fstream>
#include <iostream>
#include <iomanip>


//#include <limits>
#include "newm/newmatio.h"
#include "mlib3consoleutil.h"

//////////////////////////////////////////////////////////////////////
//
//	Math class - includes basic math functions
//
//////////////////////////////////////////////////////////////////////
const long_Real MMath::UndefDouble = -1e-36;
const long_Real MMath::PI = 3.1415926535897932384626433832795;
const long_Real MMath::RAD = 180.0/MMath::PI;

const int MMath::PLANE_YZ=0;
const int MMath::PLANE_XZ=1;
const int MMath::PLANE_XY=2;

const int MMath::COMPARE_CORR=0;
const int MMath::COMPARE_NMI=1;
const int MMath::COMPARE_ABSDIFF=2;

void MMath::GaussCoeff(Real* buf, int n)
{
	int r = n + 1, c = 1;
	buf[0] = 1;
	for (int i = 1; i<n; i++)
	{
		buf[i] = (buf[i - 1] * (r - c)) / c; c++;
	}
	buf[n] = 1;
	Real cf = pow(2.0, n);
	for (int i = 0; i <= n; i++) buf[i] = buf[i] / cf;
}
void MMath::VectorToE1(ColumnVector& A, ColumnVector& B, Matrix& T)
{
	ColumnVector v = B - A;

	Real n = v.NormFrobenius();
	Real nxy = sqrt(v(1)*v(1) + v(2)*v(2));
	Real nxz = sqrt(v(1)*v(1) + v(3)*v(3));
	Real x = v(1), y = v(2), z = v(3);
	Real cz, sz;
	if (nxy>0) { cz = x / nxy; sz = y / nxy; }
	else { cz = 1; sz = 0; }
	Real cy = nxy / n;
	Real sy = z / n;
	Matrix M1, M2, M3;
	IdentityTransform(M1); M1(1, 1) = cz; M1(1, 2) = sz; M1(2, 1) = -sz; M1(2, 2) = cz;
	IdentityTransform(M2); M2(1, 1) = cy; M2(1, 3) = sy; M2(3, 1) = -sy; M2(3, 3) = cy;
	IdentityTransform(M3); M3(1, 1) = 1 / n;

	//Translate to origin first.
	ColumnVector Aneg = -A;
	Matrix M0; TranslateMatrix(Aneg, M0);

	T = M3*M2*M1*M0;
}

void MMath::TranslateMatrix(ColumnVector& v, Matrix& T) {
	IdentityTransform(T);
	T(1, 4) = v(1); T(2, 4) = v(2); T(3, 4) = v(3);
}

void MMath::ConvolveSep3D_sm(Volume& v0, Volume& v1, Real* krn, int half_krn_sz, Volume* mask)
{
	v1 = v0; v1 = 0;
	ConvolveVolume1D(v0, v1, 1, krn, half_krn_sz, mask);
	Volume v2; v2.InitMemory(v0.m_dims);
	ConvolveVolume1D(v1, v2, 2, krn, half_krn_sz, mask);
	if (v1.SZ()>1)
		ConvolveVolume1D(v2, v1, 3, krn, half_krn_sz, mask);
}
void MMath::Convolve1D(ColumnVector& v, ColumnVector& res, Real* krn, int half_krn_sz)
{
	int sz = v.size();
	res.resize(sz);

	Real sum = 0;
	//	bool bMaskStop;
	//main part
	for (int x = half_krn_sz; x<sz - half_krn_sz; x++)
	{
		sum = 0;
		for (int i = x - half_krn_sz, j = 0; i <= x + half_krn_sz; i++, j++) sum += v(i + 1)*krn[j];
		res(x + 1) = sum;
	}
	//left
	for (int x = 0; x<half_krn_sz; x++)
	{
		sum = 0;
		for (int i = x - half_krn_sz, j = 0; i <= x + half_krn_sz; i++, j++) sum += (i<0) ? v(1)*krn[j] : v(i + 1)*krn[j];
		res(x + 1) = sum;
	}
	//right
	for (int x = sz - half_krn_sz; x<sz; x++)
	{
		sum = 0;
		for (int i = x - half_krn_sz, j = 0; i <= x + half_krn_sz; i++, j++) sum += (i >= sz) ? v(sz)*krn[j] : v(i + 1)*krn[j];
		res(x + 1) = sum;
	}
}

void MMath::ConvolveVolume1D(Volume& v0, Volume& v1, int dim, Real* krn, int half_krn_sz, Volume* mask)
{
	int i1, i2, i3, d1, d2, d3;
	v1.InitMemory(v0.m_dims);
	if (dim == 1) { i1 = 2; i2 = 1; i3 = 0; }
	else if (dim == 2) { i1 = 2; i2 = 0; i3 = 1; }
	else if (dim == 3) { i1 = 1; i2 = 0; i3 = 2; }
	d1 = v0.m_dims[i1]; d2 = v0.m_dims[i2]; d3 = v0.m_dims[i3];
	int x[3];
	ColumnVector r0, r1, mask1D;
	r0.resize(d3); r1.resize(d3); mask1D.resize(d3);
	x[i3] = 0;
	for (x[i1] = 0; x[i1]<d1; x[i1]++)
		for (x[i2] = 0; x[i2]<d2; x[i2]++)
		{
			v0.Get1D(r0, x[0], x[1], x[2], dim);
			if (mask)
			{
				mask->Get1D(mask1D, x[0], x[1], x[2], dim);
				ConvolveBinom1D_3(r0, &mask1D);
			}
			else Convolve1D(r0, r1, krn, half_krn_sz);

			v1.Set1D((mask) ? r0 : r1, x[0], x[1], x[2], dim);
		}
}
void MMath::ConvolveBinom1D_3(ColumnVector& v, ColumnVector* mask)
{
	int sz = v.size();

	Real sum = 0, v1, v2, v3, m1, m2, m3;
	//main part
	v1 = 0; v2 = v(1); v3 = v2;
	m1 = 0; m2 = (*mask)(1); m3 = m2;
	for (int x = 1; x<sz - 1; x++)
	{
		if ((*mask)(x + 1) == 0) { m1 = 0; continue; }
		m2 = (*mask)(x + 1); m3 = (*mask)(x + 2);
		v2 = v(x + 1); v3 = v(x + 2);
		if (m1)
		{
			if (m3) v(x + 1) = (v1 + v2 + v2 + v3) / 4;
			else   v(x + 1) = (v1 + v2 + v2) / 3;
		}
		else if (m3) v(x + 1) = (v2 + v2 + v3) / 3;
		v1 = v2; m1 = m2;
	}
	if ((*mask)(sz) && m1) v(sz) = (v1 + 2 * v(sz)) / 3;
}
//v1 is moving, v2 is target
bool MMath::MutualInfoNeib(Volume& v1, Volume& v2, Volume& vout, int neib_sz)
{
	if (v1.m_nElements != v2.m_nElements || neib_sz > v1.SX() / 8) return false;

	Matrix hist1, hist2, jhist;
	int nbins = 50;
	Range3 sub;
	Real res = 0, h1 = 0, h2 = 0;
	Real eps = 1e-14;
	int i, j;
	vout = v1; vout = 0;


	for (int z = neib_sz; z < v1.SZ() - neib_sz; z++) {
		cout << "slice " << z << "/" << v1.SZ() - neib_sz << "\r" << std::flush;
		for (int y = neib_sz; y < v1.SY() - neib_sz; y++) {
			for (int x = neib_sz; x < v1.SX() - neib_sz; x++) {
				sub.Set(x - neib_sz, x + neib_sz, y - neib_sz, y + neib_sz, z - neib_sz, z + neib_sz);
				v1.HistSubvol(hist1, nbins, sub); hist1 /= hist1.Sum();
				v2.HistSubvol(hist2, nbins, sub); hist2 /= hist2.Sum();
				v2.JointHistSubvol(jhist, v1, nbins, sub); jhist /= jhist.Sum();

				res = h1 = h2 = 0;
				for (i = 1; i <= nbins; i++)
				{
					h1 -= hist1(i, 1) * log(hist1(i, 1) + eps);
					h2 -= hist2(i, 1) * log(hist2(i, 1) + eps);
				}
				for (i = 1; i <= nbins; i++)
					for (j = 1; j <= nbins; j++)
						res -= jhist(i, j) * log(jhist(i, j) + eps);
				vout(x, y, z) = (h1 + h2) / res;
			}
		}
	}
	return true;
}

//return conjugate gradient computed in fixed size neighborhood, for the entire volume
bool MMath::EtaNeib(Volume& v1, Volume& v2, Volume& vout, int neib_sz)
{
	if (v1.m_nElements != v1.m_nElements || neib_sz > v1.SX() / 4) return false;
	//helper functions; potential further optimization to compute partial sums.
	struct EX {
		static void EtaLocalSum(Volume& u11, Volume& u12, Volume& u22, int x0, int y0, int z0, int neib_sz, Real& s11, Real& s12, Real& s22, int& n)
		{
			int x, y, z; n = 0; s11 = s12 = s22 = 0;
			for (z = z0 - neib_sz; z < z0 + neib_sz; z++) {
				for (y = y0 - neib_sz; y < y0 + neib_sz; y++) {
					for (x = x0 - neib_sz; x < x0 + neib_sz; x++) {
						if (u11(x, y, z) * u22(x, y, z) <= 0) continue;
						s11 += u11(x, y, z); s12 += u12(x, y, z); s22 += u22(x, y, z); n++;
					}
				}
			}
		}
	};
	Volume g1x, g1y, g1z, g2x, g2y, g2z, u11, u22, u12;
	Gradient1D(v1, g1x, 1); Gradient1D(v1, g1y, 2); Gradient1D(v1, g1z, 3);
	Gradient1D(v2, g2x, 1); Gradient1D(v2, g2y, 2); Gradient1D(v2, g2z, 3);
	int x, y, z;
	u11 = v1; u11 = 0; u22 = v1; u22 = 0; u12 = v1; u12 = 0;
	Real a11, a22, a12, q;

	//phase 1. pre-compute conjugate gradients.
	for (z = 0; z < v1.SZ(); z++) {
		for (y = 0; y < v1.SY(); y++) {
			for (x = 0; x < v1.SX(); x++) {
				a11 = g1x(x, y, z) * g1x(x, y, z) + g1y(x, y, z) * g1y(x, y, z) + g1z(x, y, z) * g1z(x, y, z);
				a22 = g2x(x, y, z) * g2x(x, y, z) + g2y(x, y, z) * g2y(x, y, z) + g2z(x, y, z) * g2z(x, y, z);
				q = a11 * a22;
				if (q <= 0) continue;
				a12 = g1x(x, y, z) * g2x(x, y, z) + g1y(x, y, z) * g2y(x, y, z) + g1z(x, y, z) * g2z(x, y, z);
				u11(x, y, z) = a11; u22(x, y, z) = a22; u12(x, y, z) = (a12 * a12) / sqrt(q);
			}
		}
	}

	Real s11 = 0, s12 = 0, s22 = 0;
	int n = 0;
	/*
		//phase 2. compute value for the first neighborhood.
		for (z = 0; z < neib_sz; z++) {
			for (y = 0; y < neib_sz; y++) {
				for (x = 0; x < neib_sz; x++) {

					if (u11(x,y,z)*u22(x,y,z)<=0) continue;
					s11+=u11(x,y,z); s12+=u12(x,y,z); s22+=u22(x,y,z); n++;
				}
			}
		}
	*/
	vout = v1; vout = 0;
	//phase 3. compute value for the remaining neighborhoods.
	for (z = neib_sz; z < v1.SZ() - neib_sz; z++) {
		cout << "slice " << z << "/" << v1.SZ() - neib_sz << "\r" << std::flush;
		for (y = neib_sz; y < v1.SY() - neib_sz; y++) {
			for (x = neib_sz; x < v1.SX() - neib_sz; x++) {
				EX::EtaLocalSum(u11, u12, u22, x, y, z, neib_sz, s11, s12, s22, n);
				vout(x, y, z) = (n > 0) ? s12 / sqrt(s11 * s22) : 0;
			}
		}
	}

	return true;
}
//return correlation coefficient computed in fixed size neighborhood for entire volume
bool MMath::CorrCoeffNeib(Volume& v1, Volume& v2, Volume& vout, int neib_sz)
{
	if (v1.m_nElements != v1.m_nElements || neib_sz > v1.SX() / 4) return false;

	//this structure is for potential further optimization, to compute partial sums.
	struct EX {
		static Real CCLocal(Volume& v1, Volume& v2, int x0, int y0, int z0, int neib_sz)
		{
			int x, y, z;
			Real n = 8 * neib_sz * neib_sz * neib_sz, s1 = 0, s2 = 0, t1, t2, acc1 = 0, acc2 = 0, acc3 = 0;
			//pass 1, sums.
			for (z = z0 - neib_sz; z < z0 + neib_sz; z++) {
				for (y = y0 - neib_sz; y < y0 + neib_sz; y++) {
					for (x = x0 - neib_sz; x < x0 + neib_sz; x++) {
						s1 += v1(x, y, z); s2 += v2(x, y, z);
					}
				}
			}
			s1 /= n; s2 /= n;
			//pass 2, cc accumulators
			for (z = z0 - neib_sz; z < z0 + neib_sz; z++) {
				for (y = y0 - neib_sz; y < y0 + neib_sz; y++) {
					for (x = x0 - neib_sz; x < x0 + neib_sz; x++) {
						t1 = v1(x, y, z) - s1; t2 = v2(x, y, z) - s2; acc1 += t1 * t2; acc2 += t1 * t1; acc3 += t2 * t2;
					}
				}
			}
			if (acc2 * acc3 == 0) return 0;
			return fabs(acc1 / sqrt(acc2 * acc3));
		}
	};

	vout = v1; vout = 0; int x, y, z;
	//phase 3. compute value for the remaining neighborhoods.
	for (z = neib_sz; z < v1.SZ() - neib_sz; z++) {
		cout << "slice " << z << "/" << v1.SZ() - neib_sz << "\r" << std::flush;
		for (y = neib_sz; y < v1.SY() - neib_sz; y++) {
			for (x = neib_sz; x < v1.SX() - neib_sz; x++) {
				vout(x, y, z) = EX::CCLocal(v1, v2, x, y, z, neib_sz);
			}
		}
	}
	return true;
}
//Refer to [1], pp. 14-15
//1. Roche A, Malandain G, Ayache N, Pennec X. Multimodal image registration by maximization of the correlation ratio: INRIA; 1998.
// in ref, Y is moving image and X is target.
// here, v1 is moving image and v2 is target.

bool MMath::CorrRatioNeib(Volume& v1, Volume& v2, Volume& vout, int neib_sz)
{
	if (v1.m_nElements != v1.m_nElements || neib_sz > v1.SX() / 4) return false;
	Range3 sub;
	int nbins = 50;
	Real npts = 8 * neib_sz * neib_sz * neib_sz;
	ColumnVector mi(nbins), si(nbins), pxi(nbins);

	Matrix hist1, hist2, hist12;
	int i, j;
	Real s2, eta, m;
	vout = v1; vout = 0;

	for (int z = neib_sz; z < v1.SZ() - neib_sz; z++) {
		cout << "slice " << z << "/" << v1.SZ() - neib_sz << "\r" << std::flush;
		for (int y = neib_sz; y < v1.SY() - neib_sz; y++) {
			for (int x = neib_sz; x < v1.SX() - neib_sz; x++) {
				sub.Set(x - neib_sz, x + neib_sz, y - neib_sz, y + neib_sz, z - neib_sz, z + neib_sz);
				v1.HistSubvol(hist1, nbins, sub); hist1 /= hist1.Sum();
				v2.HistSubvol(hist2, nbins, sub); hist2 /= hist2.Sum();
				v1.JointHistSubvol(hist12, v2, nbins, sub); hist12 /= hist12.Sum();
				mi = 0; si = 0; s2 = 0; pxi = 0; eta = 0; m = 0;

				//compute first order moments, m and s2
				for (i = 1; i <= nbins; i++)
				{
					for (j = 1; j <= nbins; j++)
					{
						mi(i) += j * hist12(i, j);
						pxi(j) += hist12(i, j);
					}
					mi(i) = (hist2(i, 1) > 0) ? mi(i) / hist2(i, 1) : 0;
					m += i * hist1(i, 1);
					s2 += i * i * hist1(i, 1);
				}
				s2 -= m * m;
				//second order moments.
				for (i = 1; i <= nbins; i++)
				{
					for (j = 1; j <= nbins; j++)
						si(i) += j * j * hist12(i, j);
					si(i) = (hist2(i, 1) > 0) ? (si(i) / hist2(i, 1) - mi(i) * mi(i)) : 0;
				}
				for (i = 1; i <= nbins; i++)
				{
					eta += si(i) * pxi(i);
				}
				vout(x, y, z) = (s2 > 0) ? (1.0 - eta / s2) : 0;
			}
		}
	}
	return true;
}

void MMath::Gradient1D(Volume& v, Volume& grad, int axis)
{
	grad = v; grad = 0;
	Real krn[3] = { -0.5, 0, 0.5 };
	for (int i = 0; i < 3; i++) krn[i] /= v.m_voxel_dims[axis - 1];
	ConvolveVolume1D(v, grad, axis, krn, 1);
}

//for the definition of eta, see:
//Rowland, D. J., Garbow, J. R., Laforest, R., & Snyder, A. Z. (2005). 
//Registration of [18F]FDG microPET and small-animal MRI. 
//Nuclear Medicine and Biology, 32(6), 567–72. 
//http://doi.org/10.1016/j.nucmedbio.2005.05.002

Real MMath::Eta(Volume& v1, Volume& v2, Volume* pmask)
{
	Volume g1x, g1y, g1z, g2x, g2y, g2z;
	Gradient1D(v1, g1x, 1); Gradient1D(v1, g1y, 2); Gradient1D(v1, g1z, 3);
	Gradient1D(v2, g2x, 1); Gradient1D(v2, g2y, 2); Gradient1D(v2, g2z, 3);
	Real q, a11, a22, a12, u11 = 0, u12 = 0, u22 = 0;
	int n = 0;

	for (int i = 0; i < v1.m_nElements; i++)
	{
		if (pmask && !(*pmask)(i)) continue;
		a11 = g1x(i) * g1x(i) + g1y(i) * g1y(i) + g1z(i) * g1z(i);
		a22 = g2x(i) * g2x(i) + g2y(i) * g2y(i) + g2z(i) * g2z(i);
		a12 = g1x(i) * g2x(i) + g1y(i) * g2y(i) + g1z(i) * g2z(i);
		q = a11 * a22;
		if (q <= 0) continue;
		u12 += (a12 * a12) / sqrt(q);
		u11 += a11;
		u22 += a22;
		n++;
	}
	return (n > 0) ? u12 / sqrt(u11 * u22) : 0;
}
