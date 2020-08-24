#include "CNDArray.h"
#include <complex>
using namespace std;

// NRvector definitions
template <class T>
NRvector<T>::NRvector() : nn(0), v(NULL) {}

template <class T>
NRvector<T>::NRvector(size_t n) : nn(n), v(n>0 ? new T[n] : NULL) {}

template <class T>
NRvector<T>::NRvector(size_t n, const T& a) : nn(n), v(n>0 ? new T[n] : NULL)
{
	for(size_t i=0; i<n; i++) v[i] = a;
}

template <class T>
NRvector<T>::NRvector(size_t n, const T *a) : nn(n), v(n>0 ? new T[n] : NULL)
{
	for(size_t i=0; i<n; i++) v[i] = *a++;
}

template <class T>
NRvector<T>::NRvector(const NRvector<T> &rhs) : nn(rhs.nn), v(nn>0 ? new T[nn] : NULL)
{
	for(size_t i=0; i<nn; i++) v[i] = rhs[i];
}

template <class T>
NRvector<T> & NRvector<T>::operator=(const NRvector<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
	if (this != &rhs)
	{
		if (nn != rhs.nn) {
			if (v != NULL) delete [] (v);
			nn=rhs.nn;
			v= nn>0 ? new T[nn] : NULL;
		}
		for (size_t i=0; i<nn; i++)
			v[i]=rhs[i];
	}
	return *this;
}

template <class T>
void NRvector<T>::resize(size_t newn)
{
	if (newn != nn) {
		if (v != NULL) delete[] (v);
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
}

template <class T>
void NRvector<T>::assign(size_t newn, const T& a)
{
	if (newn != nn) {
		if (v != NULL) delete[] (v);
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
	for (size_t i=0;i<nn;i++) v[i] = a;
}

template <class T>
NRvector<T>::~NRvector()
{
	if (v != NULL) delete[] (v);
}

// end of NRvector definitions

// NRmatrix definitions
template <class T>
NRmatrix<T>::NRmatrix() : nn(0), mm(0), v(NULL) {}

template <class T>
NRmatrix<T>::NRmatrix(size_t n, size_t m) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	size_t i,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1;i<n;i++) v[i] = v[i-1] + m;
}

template <class T>
NRmatrix<T>::NRmatrix(size_t n, size_t m, const T &a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	size_t i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
	for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = a;
}

template <class T>
NRmatrix<T>::NRmatrix(size_t n, size_t m, const T *a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	size_t i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
	for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = *a++;
}

template <class T>
NRmatrix<T>::NRmatrix(const NRmatrix &rhs) : nn(rhs.nn), mm(rhs.mm), v(nn>0 ? new T*[nn] : NULL)
{
	size_t i,j,nel=mm*nn;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
}

template <class T>
NRmatrix<T> & NRmatrix<T>::operator=(const NRmatrix<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
	if (this != &rhs) {
		size_t i,j,nel;
		if (nn != rhs.nn || mm != rhs.mm) {
			if (v != NULL) {
				delete[] (v[0]);
				delete[] (v);
			}
			nn=rhs.nn;
			mm=rhs.mm;
			v = nn>0 ? new T*[nn] : NULL;
			nel = mm*nn;
			if (v) v[0] = nel>0 ? new T[nel] : NULL;
			for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
		}
		for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
	}
	return *this;
}

template <class T>
void NRmatrix<T>::transpose()
{
	size_t i,j;
	T* vt=new T[nn*mm];
	T* vtt=vt;
	for(j=0;j<mm;j++)
		for(i=0;i<nn;i++)
			*(vtt++)=v[i][j];
	memcpy((void *)v[0],(void *)vt,nn*mm*sizeof(T));
	delete[] (vt);
}

template <class T>
void NRmatrix<T>::transpose(T* dest)
/* 
* Use this method with extremely care since 
* the size of the destination memory (dest) is
* presumed to be the same as the size of original NRmatrix object.
* 
* Recommend usage:
*   NRmatrix<T> mat(n,m)
*	// ... do something here to mat
*	NRmatrix<T> mat2(n,m);
*	mat.transpose(mat2.getRaw()); // copy the transpose of mat to mat2
*/
{
	size_t i,j;
	for(j=0;j<mm;j++)
		for(i=0;i<nn;i++)
			*(dest++)=v[i][j];
}

template <class T>
void NRmatrix<T>::transpose(NRmatrix<T> &lhs)
{
	size_t i,j;
	T* praw=lhs.getRaw();
	if(nn!=lhs.nrows() || mm!=lhs.ncols()){
		return;
	}
	for(j=0;j<mm;j++)
		for(i=0;i<nn;i++)
			*(praw++)=v[i][j];
}

template <class T>
void NRmatrix<T>::resize(size_t newn, size_t newm)
{
	size_t i,nel;
	if (newn != nn || newm != mm) {
		if (v != NULL) {
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		v = nn>0 ? new T*[nn] : NULL;
		nel = mm*nn;
		if (v) v[0] = nel>0 ? new T[nel] : NULL;
		for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	}
}

template <class T>
void NRmatrix<T>::assign(size_t newn, size_t newm, const T& a)
{
	size_t i,j,nel;
	if (newn != nn || newm != mm) {
		if (v != NULL) {
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		v = nn>0 ? new T*[nn] : NULL;
		nel = mm*nn;
		if (v) v[0] = nel>0 ? new T[nel] : NULL;
		for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	}
	for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = a;
}

template <class T>
NRmatrix<T>::~NRmatrix()
{
	if (v != NULL) {
		delete[] (v[0]);
		delete[] (v);
	}
}
// end of NRmatrix definitions

// NRMat3d definitions 
template <class T>
NRMat3d<T>::NRMat3d(): nn(0), mm(0), kk(0), v(NULL) {}

template <class T>
NRMat3d<T>::NRMat3d(size_t n, size_t m, size_t k) : nn(n), mm(m), kk(k), v(new T**[n])
{
	size_t i,j;
	v[0] = new T*[n*m];
	v[0][0] = new T[n*m*k];
	for(j=1; j<m; j++) v[0][j] = v[0][j-1] + k;
	for(i=1; i<n; i++) {
		v[i] = v[i-1] + m;
		v[i][0] = v[i-1][0] + m*k;
		for(j=1; j<m; j++) v[i][j] = v[i][j-1] + k;
	}
}

template <class T>
NRMat3d<T>::NRMat3d(size_t n, size_t m, size_t k, const T &a) : nn(n), mm(m), kk(k), v(new T**[n])
{
	size_t i,j,l;
	v[0] = new T*[n*m];
	v[0][0] = new T[n*m*k];
	for(j=1; j<m; j++) v[0][j] = v[0][j-1] + k;
	for(i=1; i<n; i++) {
		v[i] = v[i-1] + m;
		v[i][0] = v[i-1][0] + m*k;
		for(j=1; j<m; j++) v[i][j] = v[i][j-1] + k;
	}
    for(i=0;i<n;i++)
        for(j=0;j<m;j++)
            for(l=0;l<k;l++)
              v[i][j][l] = a;
}

template <class T>
NRMat3d<T>::NRMat3d(size_t n, size_t m, size_t k, const T *a) : nn(n), mm(m), kk(k), v(new T**[n])
{
	size_t i,j,l;
	v[0] = new T*[n*m];
	v[0][0] = new T[n*m*k];
	for(j=1; j<m; j++) v[0][j] = v[0][j-1] + k;
	for(i=1; i<n; i++) {
		v[i] = v[i-1] + m;
		v[i][0] = v[i-1][0] + m*k;
		for(j=1; j<m; j++) v[i][j] = v[i][j-1] + k;
	}
    for(i=0;i<n;i++)
        for(j=0;j<m;j++)
            for(l=0;l<k;l++)
              v[i][j][l] = *a++;
}

template <class T>
NRMat3d<T> & NRMat3d<T>::operator=(const NRMat3d<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
	if (this != &rhs) {
		size_t i,j,l,nel;
		if (nn !=rhs.nn || mm!=rhs.mm || kk!=rhs.kk) {
			if (v != NULL) {
			    delete[] (v[0][0]);
				delete[] (v[0]);
				delete[] (v);
			}
			nn=rhs.nn;
			mm=rhs.mm;
			kk=rhs.kk;
		    v = nn>0 ? new T**[nn] : NULL;
		    v[0] = new T*[nn*mm];
		    v[0][0] = new T[nn*mm*kk];
		    for(j=1; j<mm; j++) v[0][j] = v[0][j-1] + kk;
		    for(i=1; i<nn; i++) {
			    v[i] = v[i-1] + mm;
			    v[i][0] = v[i-1][0] + mm*kk;
			    for(j=1; j<mm; j++) v[i][j] = v[i][j-1] + kk;
		    }
        }
        for(i=0;i<nn;i++)
            for(j=0;j<mm;j++)
                for(l=0;l<kk;l++)
                    v[i][j][l] = rhs[i][j][l];
	}
	return *this;
}

template <class T>
void NRMat3d<T>::transpose()
// Cautious: the content of current object will be modified.
{
	size_t i,j,l;
	T* vt=new T[nn*mm*kk];
	T* vtt=vt;
	for(l=0;l<kk;l++)
		for(j=0;j<mm;j++)
			for(i=0;i<nn;i++)
				*(vtt++)=v[i][j][l];
	memcpy((void *)v[0][0],(void *)vt,nn*mm*kk*sizeof(T));
	delete[] (vt);
}

template <class T>
void NRMat3d<T>::transpose(T* dest)
/* 
* Use this method with extremely care since 
* the size of the destination memory (dest) is
* presumed to be the same as the size of original NRmatrix object.
* 
* Recommend usage:
*   NRMat3d<T> mat(n,m,k)
*	// ... do something here to mat
*	NRMat<T> mat2(n,m,k);
*	mat.transpose(mat2.getRaw()); // copy the transpose of mat to mat2
*/
{
	size_t i,j,l;
	for(l=0;l<kk;l++)
		for(j=0;j<mm;j++)
			for(i=0;i<nn;i++)
				*(dest++)=v[i][j][l];
}

template <class T>
void NRMat3d<T>::transpose(NRMat3d<T> &lhs)
{
	size_t i,j,l;
	T* praw=lhs.getRaw();
	if(nn!=lhs.dim1() || mm!=lhs.dim2() || kk!=lhs.dim3()){
		return;
	}
	for(l=0;l<kk;l++)
		for(j=0;j<mm;j++)
			for(i=0;i<nn;i++)
				*(praw++)=v[i][j][l];
}

template <class T>
void NRMat3d<T>::resize(size_t newn, size_t newm, size_t newk)
{
	size_t i,j;
	if (newn!=nn || newm!=mm || newk!=kk) {
		if (v != NULL) {
			delete[] (v[0][0]);
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		kk = newk;
		v = nn>0 ? new T**[nn] : NULL;
		v[0] = new T*[nn*mm];
		v[0][0] = new T[nn*mm*kk];
		for(j=1; j<mm; j++) v[0][j] = v[0][j-1] + kk;
		for(i=1; i<nn; i++) {
			v[i] = v[i-1] + mm;
			v[i][0] = v[i-1][0] + mm*kk;
			for(j=1; j<mm; j++) v[i][j] = v[i][j-1] + kk;
		}
	}
}

template <class T>
NRMat3d<T>::~NRMat3d()
{
	if (v != NULL) {
		delete[] (v[0][0]);
		delete[] (v[0]);
		delete[] (v);
	}
}
// end of NRMat3d definitions

// lyxMat4d definitions
template <class T>
lyxMat4d<T>::lyxMat4d(): nn(0), mm(0), kk(0), ww(0), v(NULL) {}

template <class T>
lyxMat4d<T>::lyxMat4d(size_t n, size_t m, size_t k, size_t w) : nn(n), mm(m), kk(k), ww(w), v(new T***[n])
{
	size_t i,j,l;
	size_t mk=m*k;
	size_t kw=k*w;
	size_t mkw=mk*w;
	
	v[0] = new T**[n*m];
	v[0][0] = new T*[n*mk];
	v[0][0][0]=new T[n*mkw];

	/* optimized version */
	/*
	for(j=1; j<m; j++){
		v[0][j] = v[0][j-1] + k;
		v[0][j][0] = v[0][j-1][0] + kw;
		for(l=1; l<k; l++) v[0][j][l]=v[0][j][l-1] + w;
	}
	for(l=1; l<k; l++) v[0][0][l] = v[0][0][l-1] + w;
	for(i=1; i<n; i++) {
		v[i] = v[i-1] + m;
		v[i][0] = v[i-1][0] + mk;
		v[i][0][0] = v[i-1][0][0] + mkw;
		for(l=1; l<k; l++) v[i][0][l] = v[i][0][l-1] + w;
		for(j=1; j<m; j++){
			v[i][j] = v[i][j-1] + k;
			v[i][j][0] = v[i][j-1][0] + kw;
			for(l=1; l<k; l++)	v[i][j][l] = v[i][j][l-1] + w;
		}
	}
	*/
	/* Much clear version */
	/* 
	*  1. To reduce the number of multiple operations, 
	*	  following constant can be calculated first:
	*			size_t mk=m*k;
	*			size_t kw=k*w;
	*			size_t mkw=mk*w;
	*  2. To further reduce the number of multiple and adding operations,
	*	  an iteration scheme should be introduced for calculating the position
	*     that the pointers point to. Namely, instead of calculating offset from
	*     the beginning position at each time:
	*				v[i]=v[0]+i*m;
	*	  the position of the pointer should be calculated based on the previous one:
	*				v[i]=v[i-1]+m;
	*	  Thus, one multiple operation can be saved. However, extremely care must be
	*     taken to make sure all starting positions are initialized before being
	*     referenced, especially for cases of more than one dimensions involved. 
	*  3. The optimized version with consideration of 1 and 2 is presented above.
	*  4. Actual timing test shows that these two methods have no significant difference.
	*		Test environment: Core i5 + DDRIII 8G, Windown 7 enterprise x64
	*		Matrix size 100*100*100*80
	*		Compile mode: Release x64
	*						Time consuming (sec)
	*		Optimized version:		0.015
	*		Much clear version:		0.015
	*		Compile mode: Debug x64
	*						Time consuming (sec)
	*		Optimized version:		0.169
	*		Much clear version:		0.165		
	*/
	/**/
	for(i=0;i<n;i++){
		// each v[i] points to a set of pointers 
		// {v[i][0], v[i][1], ..., v[i][m-1]}
		v[i]=v[0]+i*m;
		for(j=0;j<m;j++){
			// each v[i][j] points to a set of pointers 
			// {v[i][j][0], v[i][j][1], ..., v[i][j][k-1]}
			v[i][j]=v[0][0]+i*mk+j*k; 
			for(l=0;l<k;l++){
				// each v[i][j][l] points to a row of T values
				// {v[i][j][l][0], v[i][j][l][1], ..., v[i][j][l][w-1]}
				v[i][j][l]=v[0][0][0]+i*mkw+j*kw+l*w;
			}
		}
	}
	/**/
}

template <class T>
void lyxMat4d<T>::transpose()
// Cautious: the content of current object will be modified.
{
	size_t i,j,l,s;
	T* vt=new T[nn*mm*kk*ww];
	T* vtt=vt;
	for(s=0;s<ww;s++)
		for(l=0;l<kk;l++)
			for(j=0;j<mm;j++)
				for(i=0;i<nn;i++)
					*(vtt++)=v[i][j][l][s];
	memcpy((void *)v[0][0][0],(void *)vt,nn*mm*kk*ww*sizeof(T));
	delete[] (vt);
}

template <class T>
void lyxMat4d<T>::transpose(T* dest)
/* 
* Use this method with extremely care since 
* the size of the destination memory (dest) is
* presumed to be the same as the size of original NRmatrix object.
* 
* Recommend usage:
*   lyxMat4d<T> mat(n,m,k,w)
*	// ... do something here to mat
*	lyxMat4d<T> mat2(n,m,k,w);
*	mat.transpose(mat2.getRaw()); // copy the transpose of mat to mat2
*/
{
	size_t i,j,l,s;
	for(s=0;s<ww;s++)
		for(l=0;l<kk;l++)
			for(j=0;j<mm;j++)
				for(i=0;i<nn;i++)
					*(dest++)=v[i][j][l][s];
}

template <class T>
void lyxMat4d<T>::transpose(lyxMat4d<T> &lhs)
{
	size_t i,j,l,s;
	T* praw=lhs.getRaw();
	if(nn!=lhs.dim1() || mm!=lhs.dim2() || kk!=lhs.dim3() || ww!=lhs.dim4()) return;
	for(s=0;s<ww;s++)
		for(l=0;l<kk;l++)
			for(j=0;j<mm;j++)
				for(i=0;i<nn;i++)
					*(praw++)=v[i][j][l][s];
}

template <class T>
void lyxMat4d<T>::resize(size_t newn, size_t newm, size_t newk, size_t neww)
{
	size_t i,j,l;
	size_t mk,kw,mkw;
	if (newn!=nn || newm!=mm || newk!=kk || neww!=ww) {
		if (v != NULL) {
			delete[] (v[0][0][0]);
			delete[] (v[0][0]);
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		kk = newk;
		ww = neww;
		v = nn>0 ? new T***[nn] : NULL;
		mk=mm*kk;
		kw=kk*ww;
		mkw=mk*ww;
		v[0] = new T**[nn*mm];
		v[0][0] = new T*[nn*mk];
		v[0][0][0]=new T[nn*mkw];
		// much clear version
		for(i=0;i<nn;i++){
			// each v[i] points to a set of pointers 
			// {v[i][0], v[i][1], ..., v[i][m-1]}
			v[i]=v[0]+i*mm;
			for(j=0;j<mm;j++){
				// each v[i][j] points to a set of pointers 
				// {v[i][j][0], v[i][j][1], ..., v[i][j][k-1]}
				v[i][j]=v[0][0]+i*mk+j*kk; 
				for(l=0;l<kk;l++){
					// each v[i][j][l] points to a row of T values
					// {v[i][j][l][0], v[i][j][l][1], ..., v[i][j][l][w-1]}
					v[i][j][l]=v[0][0][0]+i*mkw+j*kw+l*ww;
				}
			}
		}
	}
}

template <class T>
lyxMat4d<T>::~lyxMat4d()
{
	if (v != NULL) {
		delete[] (v[0][0][0]);
		delete[] (v[0][0]);
		delete[] (v[0]);
		delete[] (v);
	}
}
/*
 * Following work around the lack of separation compilation model
 * of g++
 *
*/
typedef complex<double> dcomplx;

template class NRvector<int>;
template class NRvector<double>;
template class NRvector<dcomplx>;
template class NRvector<NRvector<Doub> *>;
template class NRvector<NRmatrix<Doub> *>;
template class NRvector<NRMat3d<Doub> *>;
//template class NRvector<fftw_complex>;

template class NRmatrix<int>;
template class NRmatrix<double>;
template class NRmatrix<float>;
template class NRmatrix<dcomplx>;
//template class NRmatrix<fftw_complex>;

template class NRMat3d<int>;
template class NRMat3d<double>;
template class NRMat3d<float>;
template class NRMat3d<dcomplx>;
//template class NRMat3d<fftw_complex>;

template class lyxMat4d<int>;
template class lyxMat4d<double>;
template class lyxMat4d<float>;
template class lyxMat4d<dcomplx>;
//template class lyxMat4d<fftw_complex>;

