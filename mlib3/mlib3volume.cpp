/*******************************************************************************
Copyright (c) 2016, 2017, 2018
Authored by: Mikhail Milchenko, mmilchenko@wustl.edu
This software is free for non-commercial use by acadimic, government, and non-profit/not-for-profit institutions. If this code is modified and/or redistributed, this note must appear in the beginning of each of re-distributed source code files. A commercial version of the software is available, and is licensed through the Office of Technology Management at Washington University School of Medicine. For more information regarding the commercial license, please email your request to otm@dom.wustl.edu.

BY DOWNLOADING THIS SOFTWARE YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS PROVIDED 'AS IS', WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT. IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
**********************************************************************************/

#include "mlib3volume.h"
#include "mlib3consoleutil.h"
#ifdef _4DFP
extern "C"
{
	#include <Getifh.h>
	#include <rec.h>
}
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

Volume::INTERP_METHOD Volume::s_interp_method=Volume::TRILINEAR;
Volume::VOXEL_CENSORING_METHOD Volume::s_vox_cens_method=Volume::NONE;
Volume::VERBOSITY_LEVEL Volume::s_verbosity = Volume::VERB_MINIMAL;

void Volume::InitVars()
{
#ifdef _4DFP
	m_pIFH=NULL;
#endif
	m_bVoxDimInited=true;
	m_voxel_dims[0]=m_voxel_dims[1]=m_voxel_dims[2]=1;
	m_pBuf=NULL;
	m_coord_format=FSL;
	m_ElCur=0;
	m_nElements=0;
	m_offset[0]=m_offset[1]=m_offset[2]=0;
}
/******************************************************
* Create volume with variable number of dimensions.
*******************************************************/
void Volume::InitMemory(int dim1, int dim2, int dim3)
{
	int dims[]={dim1, dim2, dim3}, nEl=1, i;
	//memory already allocated
	if(m_dims[0]==dim1 && m_dims[1]==dim2 && m_dims[2]==dim3 && m_pBuf) return;
	memcpy(m_dims,dims,3*sizeof(int));
	for (i=0; i<3; i++) nEl*=dims[i];
	bool bChangeSize=true;
	if(m_nElements==nEl && m_pBuf!=NULL) bChangeSize=false;		
	m_nElements=nEl;
	if(bChangeSize)
	{
		if(m_pBuf) delete[] m_pBuf;
		m_pBuf=new Real[m_nElements];
	}
	m_ElCur=0;
	(*this)=0;
}
// newW/H/D - new dimensions
// x/y/z0 - translate vector from old origin. 
// fill - fill value
void	Volume::ReSize(int newW, int newH, int newD, int x0, int y0, int z0, Real fill, Volume &resized)
{
	resized.InitMemory(newW,newH,newD);
	resized.SetVoxelDims(m_voxel_dims);
	int cx,cy,cz;
	bool wby,wbz,wbyz;
//	Real val;
	for(int uz=0; uz<newD; uz++)
	{
		cz=uz+z0;
		wbz=(cz>=0 && cz<SZ());
		for(int uy=0; uy<newH; uy++)
		{
			cy=uy+y0;
			wby=(cy>=0 && cy<SY());
			wbyz=(wbz && wby);
			if(!wbyz)
				for(int ux=0; ux<newW; ux++)	resized(ux,uy,uz)=fill;
			else 
				for(int ux=0; ux<newW; ux++)	
				{
					cx=ux+x0;
					if(cx>=0 && cx<SX())
					{
//						val=(*this)(cx,cy,cz);
						resized(ux,uy,uz)=(*this)(cx,cy,cz);
					}
					else
						resized(ux,uy,uz)=fill;
				}
		}
	}
}

Real Volume::Interp3(ColumnVector& X)
{
		if (s_interp_method==NNEIB)
		{
			if (!IsInsideBoundary(round(X(1)), round(X(2)), round(X(3)))) return 0;
			return (*this)(round(X(1)),round(X(2)),round(X(3)));
		}
		static Real z_old, y_old, x_old, z_new,x_new,y_new,x0,x1,y0,y1,z0,z1,xd,yd,zd,xd1,yd1,zd1;
		static Real c00,c10,c01,c11,c0,c1;
		z_old=X(3);
		z0=floor(z_old); z1=_mn(SZ()-1,_ceil(z_old));
		zd=(z_old-z0); zd1=1-zd;
		y_old=X(2);
		y0=floor(y_old); y1=_mn(SY()-1,_ceil(y_old));
		yd=(y_old-y0); yd1=1-yd;
		x_old=X(1);
		x0=floor(x_old); x1=_mn(SX()-1,_ceil(x_old));
		xd=(x_old-x0);xd1=1-xd;
		c00=(*this)(x0,y0,z0)*xd1+(*this)(x1,y0,z0)*xd;
		c10=(*this)(x0,y1,z0)*xd1+(*this)(x1,y1,z0)*xd;
		c01=(*this)(x0,y0,z1)*xd1+(*this)(x1,y0,z1)*xd;
		c11=(*this)(x0,y1,z1)*xd1+(*this)(x1,y1,z1)*xd;
		c0=c00*yd1+c10*yd;
		c1=c01*yd1+c11*yd;
		return c0*zd1+c1*zd;
	};

bool	Volume::ReadAnalyzeHeader(char *fname)
{
	Analyze an;
	if ( ! an.ReadHeaderOnly(fname) ) return false;
	InitMemory(an,NULL);
	return true;
}
bool	Volume::ReadAnalyze(char* fname)
{
	Analyze an;
	unsigned char *imbuf;
	if(!(imbuf=an.ReadAll(fname))) 
	{
		if(imbuf) delete[] imbuf;
		return false;
	}
	InitMemory(an,imbuf);
	delete[] imbuf;
	return true;
}
void Volume::Swap4(char* a)
{
	char t;
	t = a[0]; a[0] = a[3]; a[3] = t;
	t = a[1]; a[1] = a[2]; a[2] = t;
}

void	Volume::Get1D(Real* buf, int x, int y, int z, int dim)
{
	if(dim==1) //x
	{
		memcpy(buf,m_pBuf+m_dims[0]*m_dims[1]*z+m_dims[0]*y,m_dims[0]*sizeof(Real));						
	}
	else if(dim==2) //y
	{
		Real* st=m_pBuf+m_dims[0]*m_dims[1]*z+x, *p;
		int i;
		for(i=0, p=st; i<m_dims[1]; i++, p+=m_dims[0]) buf[i]=p[0];
	}
	else if(dim==3) //z
	{
		Real* st=m_pBuf+m_dims[0]*y+x, *p;
		int pl=m_dims[0]*m_dims[1], i;
		for(i=0, p=st; i<m_dims[2]; i++, p+=pl) buf[i]=p[0];
	}
}

void	Volume::Set1D(Real* buf, int x, int y, int z, int dim)
{
	if(dim==1) //x
	{
		memcpy(m_pBuf+m_dims[0]*m_dims[1]*z+m_dims[0]*y,buf,m_dims[0]*sizeof(Real));
	}
	else if(dim==2) //y
	{
		Real* st=m_pBuf+m_dims[0]*m_dims[1]*z+x, *p;
		int i;
		for(i=0, p=st; i<m_dims[1]; i++, p+=m_dims[0]) p[0]=buf[i];
	}
	else if(dim==3) //z
	{
		Real* st=m_pBuf+m_dims[0]*y+x, *p;
		int pl=m_dims[0]*m_dims[1],i;
		for(i=0, p=st; i<m_dims[2]; i++, p+=pl) p[0]=buf[i];
	}
}

#ifdef _4DFP
void Volume::GenerateIFH(IFH* pIFH, string& program)
{
	m_pIFH = (pIFH) ? new IFH((IFH&)(*pIFH)) : new IFH();
	strcpy(m_pIFH->conversion_program, program.c_str());
	m_pIFH->matrix_size[0] = SX(); m_pIFH->matrix_size[1] = SY(); m_pIFH->matrix_size[2] = SZ();
	m_pIFH->scaling_factor[0] = m_voxel_dims[0];
	m_pIFH->scaling_factor[1] = m_voxel_dims[1];
	m_pIFH->scaling_factor[2] = m_voxel_dims[2];
}
bool	Volume::Read4dfp(char* root0)
{
//	if(cc.m_bAHeader) return ReadAnalyze(cc.m_root);

	char name[MAXL];
	char root[MAXL];
	char* occ=strstr(root0,".4dfp");

	if ( occ && strlen(occ)==5 )
	{
		strncpy(root, root0, strlen(root0)-5);
		root[strlen(root0) - 5]=0;
	}
	else strcpy(root, root0);

	sprintf (name,"%s.4dfp.ifh",root);
	if (!m_pIFH) m_pIFH=new IFH;
	if(Getifh(name, m_pIFH) != 0) { delete m_pIFH; m_pIFH=NULL; return false; }
	bool isbig=(m_pIFH->imagedata_byte_order[0]=='b');
	InitMemory(m_pIFH->matrix_size);
	if(m_nElements<1) return false;
	bool bSwap=(!IsBigEndian() && isbig) || (IsBigEndian() && !isbig);

	sprintf(name,"%s.4dfp.img",root);
	FILE* fp;
	if(!(fp=fopen(name,"rb"))) return false;
	if(sizeof(Real)==4)
	{
		if(!(fread(m_pBuf,4,m_nElements,fp)==m_nElements)) {fclose(fp); return false;}
		if(bSwap) for(int i=0; i<m_nElements; i++) Swap4((char*)(&m_pBuf[i]));
	}
	else if(sizeof(Real)==8)
	{
		float *buf=new float[m_nElements];
		if(!(fread(buf,4,m_nElements,fp)==m_nElements)) {fclose(fp); return false;}
		if(bSwap) for(int i=0; i<m_nElements; i++) Swap4((char*)(&buf[i]));
		for(int i=0; i<m_nElements; i++)
		{
			m_pBuf[i]=buf[i];
			i=i;
		}
		delete[] buf;
	}
	fclose(fp);	
	for(int i=0; i<3; i++) m_voxel_dims[i]=fabs((Real)(m_pIFH->mmppix[i]));
	return true;
}
bool Volume::Write4dfp(char* root, int argc, char *argv[])
{
	char name[MAXL];
	bool bSwap=false, bBigEndian=IsBigEndian(), bRes=true;
	if (!m_pIFH) return false;

	if(((m_pIFH->imagedata_byte_order[0]=='b') && !bBigEndian) ||
		((m_pIFH->imagedata_byte_order[0]=='l') && bBigEndian)) bSwap=true;
	sprintf (name,"%s.4dfp",root);
	bRes=WriteAnalyze(name, bSwap);
	if (bRes) bRes &= (0==Writeifh(argv[0],name,m_pIFH,IsBigEndian()? 'b':'l'));
	char rcsid[]="$Id: " __FILE__ ", " __DATE__ ", " __TIME__;
	bRes=(0==startrece(name,argc,argv,rcsid,IsBigEndian()? 'b':'l'));
	if (bRes) bRes &= (0==endrec());
	return bRes;
}
#endif
bool	Volume::WriteAnalyze(char* fname, bool bSwap/*=false*/)
{
	Analyze an;
	an.SetDefaultHeader();

	an.m_dsr.dime.dim[0]=4;
	an.m_dsr.dime.dim[1]=m_dims[0];
	an.m_dsr.dime.dim[2]=m_dims[1];
	an.m_dsr.dime.dim[3]=m_dims[2];
	an.m_dsr.dime.dim[4]=1;
	an.m_dsr.dime.datatype=Analyze::DT_FLOAT;
	an.m_dsr.dime.bitpix=32;

	an.m_dsr.dime.pixdim[1]=(float)m_voxel_dims[0];
	an.m_dsr.dime.pixdim[2]=(float)m_voxel_dims[1];
	an.m_dsr.dime.pixdim[3]=(float)m_voxel_dims[2];

	if(sizeof(Real)==8) //double, 8-byte
	{
		float* tmp=new float[m_nElements];
		for(int i=0; i< m_nElements; i++) 
		{
			tmp[i]=(float)(m_pBuf[i]);
			if(bSwap) Swap4((char*)(&tmp[i]));
		}
		bool res=an.WriteAll(fname, (unsigned char*)tmp, bSwap);
		delete[] tmp;
		return res;
	}
	else //float, 4-byte 
	{
		if(bSwap) for (int i=0; i<m_nElements; i++) Swap4((char*)(&m_pBuf[i]));
		return an.WriteAll(fname, (unsigned char*)m_pBuf, bSwap);
	}
}

/******************************************************
* Create an array from Analyze header and pixel buffer.
*******************************************************/
void Volume::InitMemory(Analyze& an, void* pBuf/*=NULL*/, int nSlice/*=-1*/)
{
	int dims[3];
	int i;
	for (i=1; i<4; i++) dims[i-1]=an.m_dsr.dime.dim[i];
	bool bFull=(nSlice<0);
	if(bFull && pBuf) InitMemory(dims);
	else {m_dims[0]=dims[0]; m_dims[1]=dims[1]; m_dims[2]=dims[2];}

	int ind1,ind2,ind3;

	if(an.m_dsr.dime.pixdim[1]*an.m_dsr.dime.pixdim[2]*an.m_dsr.dime.pixdim[3]>0)
	{
		m_bVoxDimInited=true;
		m_voxel_dims[0]=an.m_dsr.dime.pixdim[1];
		m_voxel_dims[1]=an.m_dsr.dime.pixdim[2];
		m_voxel_dims[2]=an.m_dsr.dime.pixdim[3];
	}
	if(!pBuf) return;
	if(an.m_dsr.dime.datatype==Analyze::DT_FLOAT)
	{
		float* tmpB=(float*) pBuf;
		for (int z=(bFull? 0:nSlice);z<(bFull?m_dims[2]:(nSlice+1)); z++)
		{
			ind2=m_dims[1]*m_dims[0]*z;
			for(int y=0; y<m_dims[1]; y++)
			{
				ind1=ind2+m_dims[0]*y;
				ind3=(bFull?ind2:0)+m_dims[0]*y;
				for(int x=0; x<m_dims[0]; x++)
				{
					m_pBuf[ind1+x]=(Real)(tmpB[ind3+x]);
				}
			}
		}
	}
	else if(an.m_dsr.dime.datatype==Analyze::DT_SIGNED_INT)
	{
		int* tmpB=(int*) pBuf;
		for (int z=(bFull? 0:nSlice);z<(bFull?m_dims[2]:(nSlice+1)); z++)
		{
			ind2=m_dims[1]*m_dims[0]*z;
			for(int y=0; y<m_dims[1]; y++)
			{
				ind1=ind2+m_dims[0]*y;
				ind3=(bFull?ind2:0)+m_dims[0]*y;
				for(int x=0; x<m_dims[0]; x++)
				{
					m_pBuf[ind1+x]=(Real)(tmpB[ind3+x]);
				}
			}
		}
	}

	else if(an.m_dsr.dime.datatype==Analyze::DT_SIGNED_SHORT)
	{
		short* tmpB=(short*) pBuf;
		for (int z=(bFull? 0:nSlice);z<(bFull?m_dims[2]:(nSlice+1)); z++)
		{
			ind2=m_dims[1]*m_dims[0]*z;
			for(int y=0; y<m_dims[1]; y++)
			{
				ind1=ind2+m_dims[0]*y;
				ind3=(bFull?ind2:0)+m_dims[0]*y;
				for(int x=0; x<m_dims[0]; x++)
				{
					m_pBuf[ind1+x]=(Real)(tmpB[ind3+x]);
				}
			}
		}
	}
	else if(an.m_dsr.dime.datatype==Analyze::DT_UNSIGNED_CHAR)
	{
		unsigned char* tmpB=(unsigned char*) pBuf;
		for (int z=(bFull? 0:nSlice);z<(bFull?m_dims[2]:(nSlice+1)); z++)
		{
			ind2=m_dims[1]*m_dims[0]*z;
			for(int y=0; y<m_dims[1]; y++)
			{
				ind1=ind2+m_dims[0]*y;
				ind3=(bFull?ind2:0)+m_dims[0]*y;
				for(int x=0; x<m_dims[0]; x++)
				{
					m_pBuf[ind1+x]=(Real)(tmpB[ind3+x]);
				}
			}
		}
	}
	else if(an.m_dsr.dime.datatype==Analyze::DT_DOUBLE)
	{
		double* tmpB=(double*) pBuf;
		for (int z=(bFull? 0:nSlice);z<(bFull?m_dims[2]:(nSlice+1)); z++)
		{
			ind2=m_dims[1]*m_dims[0]*z;
			for(int y=0; y<m_dims[1]; y++)
			{
				ind1=ind2+m_dims[0]*y;
				ind3=(bFull?ind2:0)+m_dims[0]*y;
				for(int x=0; x<m_dims[0]; x++)
				{
					m_pBuf[ind1+x]=double(tmpB[ind3+x]);
				}
			}
		}
	}
	m_ElCur=0;
}
Volume& Volume::operator = (Volume& a)
{
	InitMemory(a.m_dims);
	memcpy(m_pBuf,a.m_pBuf,sizeof(Real)*m_nElements);
	memcpy(m_dims,a.m_dims,sizeof(int)*3);
	m_ElCur=a.m_ElCur;
	m_bVoxDimInited=a.m_bVoxDimInited;
	memcpy(m_voxel_dims,a.m_voxel_dims,sizeof(double)*3);
	memcpy(m_offset,a.m_offset,sizeof(double)*3);
	m_nElements=a.m_nElements;
	m_coord_format=a.m_coord_format;
#ifdef _4DFP
	if (a.m_pIFH)
	{
		m_pIFH=new IFH;
		memcpy(m_pIFH,a.m_pIFH,sizeof(IFH));
	}
#endif
	return (*this);
}
#define _r(x) (int)(x+0.5)

void Volume::AddPt(ColumnVector& pt, Real val, Real bg)
{
#define _n(x,y,z) sqrt((x)*(x)+(y)*(y)+(z)*(z))
#define _rr(x) (int)((x)+.5)
	static ColumnVector X[3], DX[3];
	for(int i=0; i<3; i++){
		X[i].ReSize(3); DX[i].ReSize(3);
	}
	static int neib[27][3]={
		{0,0,0},{1,0,0},{2,0,0},{0,1,0},{1,1,0},{2,1,0},{0,2,0},{1,2,0},{2,2,0},
		{0,0,1},{1,0,1},{2,0,1},{0,1,1},{1,1,1},{2,1,1},{0,2,1},{1,2,1},{2,2,1},
		{0,0,2},{1,0,2},{2,0,2},{0,1,2},{1,1,2},{2,1,2},{0,2,2},{1,2,2},{2,2,2}
	};

	for(int i=1;i<=3; i++){
		X[0](i)=_rr(pt(i))-1;
		X[1](i)=_rr(pt(i));
		X[2](i)=_rr(pt(i))+1;
		DX[0](i)=pt(i)-X[0](i);
		DX[1](i)=pt(i)-X[1](i);
		DX[2](i)=pt(i)-X[2](i);
	}
	Real D, dVal=val-bg;
	int x,y,z;
	for (int i=0; i<27; i++){
		x=neib[i][0];
		y=neib[i][1];
		z=neib[i][2];
		D=_n(DX[x](1),DX[y](2),DX[z](3));
		(*this)(X[x](1),X[y](2),X[z](3))+=exp(-D*D)*dVal+bg;
	}
}

void Volume::AddLine(ColumnVector& A, ColumnVector& B, Real val)
{
	ColumnVector Vab=B-A;
	int i,j=0;
	Real L=_r(Vab.MaximumAbsoluteValue1(i));
	ColumnVector X(3),Xst(3);
	X=A;
	Xst=Vab/L;
	do{
		AddPt(X,val,0);
		X+=Xst;
		j++;
	}while(j<=L);
}
#ifdef _4DFP
bool	Volume::Vrtflip()
{
	if (!m_pIFH) return false;
	int k=m_pIFH->orientation-1;
	if ( k < 0 ) return false;
	Matrix flips(3,3);
	float *mmpixt=(m_pIFH->mmppix), *centert=(m_pIFH->center);
	flips << -1 << 1 << -1 << -1 << 1 << 1 << 1 << 1 << 1;
	for (int i=1; i<=3; i++)
	{
		mmpixt[i-1]*=flips(k,i);
		centert[i-1]*=flips(k,i);
		if (flips(k,i)<0) 
//			centert[i-1]=mmpixt[i-1]*m_pIFH->matrix_size[i-1]+1-centert[i-1];
			centert[i-1]=mmpixt[i-1]*(m_pIFH->matrix_size[i-1]+1)-centert[i-1];
	}
	return true;
}
#endif
/////////////////////////////////////////////
// Minimum and maximum
void Volume::Stats(Matrix& stats,Volume* mask)
{
	Real mn, mx;
	mn=mx=m_pBuf[0];
	int indmn=0, indmx=0;
	for (int i=0; i<m_nElements; i++) 
	{
		if (mask && !(*mask)(i)) continue;
		if(m_pBuf[i]<mn){ mn=m_pBuf[i]; indmn=i;};
		if(m_pBuf[i]>mx){ mx=m_pBuf[i]; indmx=i;};
	}
	stats.resize(2,2);
	stats(1,1)=mn;
	stats(1,2)=(Real)indmn;
	stats(2,1)=mx;
	stats(2,2)=(Real)indmx;
}
Real	Volume::CorrCoeff(Volume& v1, Volume* mask/*=NULL*/)
{
	if (m_nElements != v1.m_nElements) return 0;
	Real nEl = (Real)(m_nElements);
	if (mask) nEl = ((*mask) != 0);

	Real s1 = Sum(mask) / nEl, s2 = v1.Sum(mask) / nEl;
	Real acc1 = 0, acc2 = 0, acc3 = 0, tmp1, tmp2;
	for (int i = 0; i < m_nElements; i++)
	{
		if (mask && !(*mask)(i)) continue;
		tmp1 = El(i) - s1;
		tmp2 = v1.El(i) - s2;
		acc1 += tmp1 * tmp2;
		acc2 += tmp1 * tmp1;
		acc3 += tmp2 * tmp2;
	}
	if (acc2 * acc3 == 0) return 0;
	return acc1 / sqrt(acc2 * acc3);
}

Real	Volume::Sorensen(Volume& v1, Real thresh/*=0*/, Volume* mask/*=NULL*/)
{
	int nMatch = 0, matchA = 0, matchB = 0;


	for (int i = 0; i < m_nElements; i++)
	{
		if (mask && !(*mask)(i)) continue;
		if (m_pBuf[i] > thresh && v1.m_pBuf[i] > thresh) nMatch++;
		if (m_pBuf[i] > thresh) matchA++;
		if (v1.m_pBuf[i] > thresh) matchB++;
	}
	return 2.0 * (Real)nMatch / (Real)(matchA + matchB);
}

Real	Volume::Sorensen1(Volume& v1, Real thresh)
{
	int nMatch = 0, matchA = 0, matchB = 0;
	for (int i = 0; i < m_nElements; i++)
	{
		if (m_pBuf[i] == thresh && v1.m_pBuf[i] == thresh) nMatch++;
		if (m_pBuf[i] == thresh) matchA++;
		if (v1.m_pBuf[i] == thresh) matchB++;
	}
	return 2.0 * (Real)nMatch / (Real)(matchA + matchB);
}

Real	Volume::NMI(Volume& v, int nbins/*=100*/, Volume* mask/*=NULL*/)
{
	if (v.m_nElements != m_nElements) return -1;
	Matrix hist1, hist2, jhist;
	Hist(hist1, nbins, mask);
	v.Hist(hist2, nbins, mask);
	JointHist(jhist, v, nbins, mask);
	hist1 /= hist1.Sum();
	//	cout << hist1 << endl;
	hist2 /= hist2.Sum();
	//	cout << hist2 << endl;
	jhist /= jhist.Sum();
	//	cout << jhist << endl;
	Real res = 0, h1 = 0, h2 = 0;
	Real eps = 1e-14;

	for (int i = 1; i <= nbins; i++)
	{
		h1 -= hist1(i, 1) * log(hist1(i, 1) + eps);
		h2 -= hist2(i, 1) * log(hist2(i, 1) + eps);
	}

	for (int i = 1; i <= nbins; i++)
		for (int j = 1; j <= nbins; j++)
		{
			//			t=hist1(i,1)*hist2(j,1);
			//			res+=jhist(i,j)*log(eps+.5*jhist(i,j)/(eps+hist1(i,1)*hist2(j,1)));
			res -= jhist(i, j) * log(jhist(i, j) + eps);
		}
	return (h1 + h2) / res;
	//		return res;
}

void Volume::Hist(Matrix& hist, int nbins/*=100*/, Volume* mask/*=NULL*/)
{
	Matrix stats(2, 2);
	Stats(stats, mask);
	Real vmin = stats(1, 1), vmax = stats(2, 1);
	hist.ReSize(nbins, 1);
	hist = 0;
	Real span = (nbins - 1) / (vmax - vmin);
	long val;
	for (int i = 0; i < m_nElements; i++)
	{
		if (mask && !(*mask)(i)) continue;
		val = (long)((m_pBuf[i] - vmin) * span + .5) + 1;
		if (val < 1) val = 1;
		if (val > nbins) val = nbins;
		hist(val, 1)++;
	}
}
void Volume::JointHist(Matrix& hist, Volume& vol, int nbins/*=100*/, Volume* mask)
{
	Matrix stats1(2, 2), stats2(2, 2);
	Stats(stats1, mask);
	vol.Stats(stats2, mask);
	Real vmin1 = stats1(1, 1), vmax1 = stats1(2, 1),
		vmin2 = stats2(1, 1), vmax2 = stats2(2, 1);
	hist.ReSize(nbins, nbins);
	hist = 0;
	Real span1 = (nbins - 1) / (vmax1 - vmin1),
		span2 = (nbins - 1) / (vmax2 - vmin2);
	long val1, val2;
	for (int i = 0; i < m_nElements; i++)
	{
		if (mask && !(*mask)(i)) continue;
		val1 = (long)((m_pBuf[i] - vmin1) * span1 + .5) + 1;
		if (val1 < 1) val1 = 1;
		if (val1 > nbins) val1 = nbins;

		val2 = (long)((vol[i] - vmin2) * span2 + .5) + 1;
		if (val2 < 1) val2 = 1;
		if (val2 > nbins) val2 = nbins;

		hist(val1, val2)++;
	}
}

bool	Volume::Uncrop(Volume& cropped, int& rx, int& ry, int& rz)
{
#define VU_SAMPLE_NUM 100

	int spts = VU_SAMPLE_NUM;
	//	int spts[VU_SAMPLE_NUM],ind,inc=m_nElements/(VU_SAMPLE_NUM+1);
	//	for (int i=0; i<VU_SAMPLE_NUM; i+=inc)
	//	{
	//		spts[inc]=inc*i+inc/2;
	//	}
	int dx = SX() - cropped.SX(), dy = SY() - cropped.SY(), dz = SZ() - cropped.SZ();
	for (int z = 0; z <= dz; z++)
		for (int y = 0; y <= dy; y++)
			for (int x = 0; x <= dx; x++)
			{
				if (CheckMatch(cropped, x, y, z, spts, false))
				{
					if (CheckMatch(cropped, x, y, z, spts, true))
					{
						//						res.InitMemory(m_dims);
						rx = -x; ry = -y; rz = -z;
						return true;
					}
				}
			}
	return false;
}
void Volume::Permute(int* direct, int* inverse)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (direct[j] == i)
			{
				inverse[i] = j;
				break;
			}
		}
	}
}
void Volume::Flip(int axis)
{
	int ax[3], rx[3];
	switch (axis)
	{
	case 1:
		ax[0] = 0; ax[1] = 1; ax[2] = 2;
		//			rx[0]=0;rx[1]=1;rx[2]=2;
		break;
	case 2:
		ax[0] = 1; ax[1] = 0; ax[2] = 2;
		//			rx[0]=1;rx[1]=0;rx[2]=2;
		break;
	case 3:
		ax[0] = 2; ax[1] = 0; ax[2] = 1;
		//			rx[0]=1;rx[1]=2;rx[2]=0;
		break;
	}
	Permute(ax, rx);

	Real temp;
	int c[3], ic[3];
	for (c[2] = 0, ic[2] = 0; c[2] < m_dims[ax[2]]; c[2]++, ic[2]++)
		for (c[1] = 0, ic[1] = 0; c[1] < m_dims[ax[1]]; c[1]++, ic[1]++)
			for (c[0] = 0, ic[0] = m_dims[ax[0]] - 1; c[0] < m_dims[ax[0]] / 2; c[0]++, ic[0]--)
			{
				temp = (*this)(ic[rx[0]], ic[rx[1]], ic[rx[2]]);
				(*this)(ic[rx[0]], ic[rx[1]], ic[rx[2]]) = (*this)(c[rx[0]], c[rx[1]], c[rx[2]]);
				(*this)(c[rx[0]], c[rx[1]], c[rx[2]]) = temp;
			}
}

// cropped - cropped volume
// dx,dy,dz - distance from origin, in corr. direction
// spts - number of sample points
// bFull - test complete volume, not just sample points.

bool	Volume::CheckMatch(Volume& cropped, int dx, int dy, int dz, int spts, bool bFull/*=false*/)
{
	if (!bFull)
	{
		int dimpts = spts / 3;
		int incx = max(1, cropped.SX() / (dimpts + 1)), incy = max(1, cropped.SY() / (dimpts + 1)), incz = max(1, cropped.SZ() / (dimpts + 1));
		//		Real cv,tv;
		for (int z = incz / 2; z < cropped.SZ(); z = z + incz)
			for (int y = incy / 2; y < cropped.SY(); y = y + incy)
				for (int x = incx / 2; x < cropped.SX(); x = x + incx)
				{
					//					cv=cropped(x,y,z),tv=(*this)(x+dx,y+dy,z+dz);
					if (cropped(x, y, z) != (*this)(x + dx, y + dy, z + dz)) return false;
				}
		return true;
	}
	else
	{
		//		Real cv,tv;
		for (int z = 0; z < cropped.SZ(); z++)
			for (int y = 0; y < cropped.SY(); y++)
				for (int x = 0; x < cropped.SX(); x++)
				{
					//					cv=cropped(x,y,z),tv=(*this)(x+dx,y+dy,z+dz);
					if (cropped(x, y, z) != (*this)(x + dx, y + dy, z + dz)) return false;
				}
		return true;
	}
}

void Volume::HistSubvol(Matrix& hist, int nbins, Range3& sub)
{
	Matrix stats(2, 4);
	StatsSubvol(stats, sub);
	Real vmin = stats(1, 1), vmax = stats(2, 1);
	hist.ReSize(nbins, 1);
	hist = 0;
	Real span = (nbins - 1) / (vmax - vmin);
	long val;
	for (int z = (int)sub.r[2].st; z < (int)sub.r[2].en; z++)
		for (int y = (int)sub.r[1].st; y < (int)sub.r[1].en; y++)
			for (int x = (int)sub.r[0].st; x < (int)sub.r[0].en; x++)
			{
				val = (long)(((*this)(x, y, z) - vmin) * span + .5) + 1;
				if (val < 1) val = 1;
				if (val > nbins) val = nbins;
				hist(val, 1)++;
			}
}
void Volume::JointHistSubvol(Matrix& hist, Volume& vol, int nbins, Range3& sub)
{
	Matrix stats1(2, 2), stats2(2, 2);
	StatsSubvol(stats1, sub);
	vol.StatsSubvol(stats2, sub);
	Real vmin1 = stats1(1, 1), vmax1 = stats1(2, 1),
		vmin2 = stats2(1, 1), vmax2 = stats2(2, 1);
	hist.ReSize(nbins, nbins);
	hist = 0;
	Real span1 = (nbins - 1) / (vmax1 - vmin1),
		span2 = (nbins - 1) / (vmax2 - vmin2);
	long val1, val2;
	for (int z = (int)sub.r[2].st; z < sub.r[2].en; z++)
		for (int y = (int)sub.r[1].st; y < sub.r[1].en; y++)
			for (int x = (int)sub.r[0].st; x < sub.r[0].en; x++)
			{
				val1 = (long)(((*this)(x, y, z) - vmin1) * span1 + .5) + 1;
				if (val1 < 1) val1 = 1;
				if (val1 > nbins) val1 = nbins;

				val2 = (long)((vol(x, y, z) - vmin2) * span2 + .5) + 1;
				if (val2 < 1) val2 = 1;
				if (val2 > nbins) val2 = nbins;

				hist(val1, val2)++;
			}
}

/////////////////////////////////////////////
// Minimum and maximum
void Volume::StatsSubvol(Matrix& stats, Range3& sub)
{
#define _FLAT(x,y,z) (x)+(m_dims[0]*(y)) + (z)*m_dims[1]*m_dims[0]
	Real val, mn, mx, avg = 0, cnt = 1e-6;
	mn = mx = m_pBuf[0];
	Index3 indmn, indmx; indmn.Set(0, 0, 0); indmx.Set(0, 0, 0);


	for (int z = (int)sub[2].st; z <= sub[2].en; z++)
		for (int y = (int)sub[1].st; y <= sub[1].en; y++)
			for (int x = (int)sub[0].st; x <= sub[0].en; x++)
			{
				val = (*this)(x, y, z);
				if (val < mn) { mn = val; indmn.Set(x, y, z); }
				if (val > mx) { mx = val; indmx.Set(x, y, z); }
				avg += val; cnt++;
			}
	stats.resize(3, 4);
	stats(1, 1) = mn;
	stats(1, 2) = (Real)indmn.d[0];	stats(1, 3) = (Real)indmn.d[1]; stats(1, 4) = (Real)indmn.d[2];
	stats(2, 1) = mx;
	stats(2, 2) = (Real)indmx.d[0]; stats(2, 3) = (Real)indmx.d[1]; stats(2, 4) = (Real)indmx.d[2];
	stats(3, 1) = avg / cnt;
}


