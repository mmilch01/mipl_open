#ifndef MLIB3_ARRAY_H
#define MLIB3_ARRAY_H

#include <stdlib.h>
#include <time.h>
#include <vector>
#include <algorithm>

const int nML3GrowBy = 64;
const int nML3SearchThreshold = 8;
const int nML3SortThreshold = 8;

template <class TYPE> class ML3Array
{
public:
	ML3Array();
	ML3Array(const TYPE* pData, const int nSize, const int nGrowBy = nML3GrowBy);
	virtual ~ML3Array();
//	int Locate(const TYPE& tData) const;
	int Find(const TYPE& tData) const;
	void SetSize(const int nSize, const int nGrowBy = nML3GrowBy);
	int GetSize() const;
	TYPE* GetDataPtr() const;
	const TYPE& operator()(const int nIndex) const;
	TYPE& operator[](const int nIndex);
	int Add(const TYPE& tData);
	int Insert(const int nIndex, const TYPE& tData);
	int Remove(const int nIndex);
	int SetAt(const int nIndex, const TYPE& tData);
	void Push(const TYPE& tData);
	TYPE Pop();
	void QuickSort(bool (*pSortFunction)(TYPE& , TYPE&));
	void Sort(int (*pSortFunction)(TYPE& , TYPE&), bool bDesc=false);
//	virtual void SortIndex(int* pIndex, const int nSort = 1) const;

private:
	int m_nSize, m_nStoreSize, m_nSort, m_nGrowBy;
	TYPE* m_pData;
};

template <class TYPE>
ML3Array<TYPE>::ML3Array()
{
	m_nSize = m_nStoreSize = 0;
	m_nGrowBy = nML3GrowBy;
	m_pData = (TYPE*)0;
	m_nSort = 0;
};

template <class TYPE>
ML3Array<TYPE>::ML3Array(const TYPE* pData, const int nSize, const int nGrowBy)
{
	m_nSize = m_nStoreSize = 0;
	m_nGrowBy = nML3GrowBy;
	m_pData = (TYPE*)0;
	m_nSort = 0;
	SetSize(nSize,nGrowBy);
	for(int i=0;i<m_nSize;i++) m_pData[i] = pData[i];
}

template <class TYPE>
ML3Array<TYPE>::~ML3Array()
{
	delete []m_pData;
};

template <class TYPE>
void ML3Array<TYPE>::SetSize
(
	const int nSize,
	const int nGrowBy
)
{
	if(nSize<0||nGrowBy<=0) return;

	if(nSize>m_nSize) m_nSort = 0;

	m_nGrowBy = nGrowBy;

	if(nSize>(m_nSize-m_nGrowBy)&&nSize<=m_nStoreSize)
	{
		m_nSize = nSize;
		return;
	}

	int nCount = nSize>m_nSize?m_nSize:nSize;
	m_nSize = nSize;
	m_nStoreSize = 0;
	while(m_nStoreSize<m_nSize) m_nStoreSize += m_nGrowBy;

	TYPE* pData = new TYPE[m_nStoreSize];
	for(int i=0;i<nCount;i++) pData[i] = m_pData[i];
	delete []m_pData;
	m_pData = pData;
};

template <class TYPE>
int ML3Array<TYPE>::GetSize() const
{
	return m_nSize;
}

template <class TYPE>
TYPE* ML3Array<TYPE>::GetDataPtr() const
{
	return m_pData;
}

template <class TYPE>
int ML3Array<TYPE>::Add(const TYPE& tData)
{
//	int nSort = m_nSort;
//	int nIndex = Locate(tData);
	Insert(m_nSize, tData);
//	m_nSort = nSort;
	return m_nSize;
}
/*
template <class TYPE>
int ML3Array<TYPE>::Locate(const TYPE& tData) const
{
	int nIndex = -1;
	int nMin, nMax, nMid;
	if(m_nSort==0||m_nSize==0) nIndex = m_nSize;
	else
	{
		nMin = 0;
		nMax = m_nSize-1;
		while(nMax-nMin>1)
		{
			nMid = (nMax+nMin)/2;
			if(m_nSort>0)
			{
				if(tData<m_pData[nMid]) nMax = nMid;
				else nMin = nMid;
			}
			else
			{
				if(m_pData[nMid]<tData) nMax = nMid;
				else nMin = nMid;
			}
		}
		if(m_nSort>0)
		{
			if(tData<m_pData[nMin]) nIndex = nMin;
			else if(tData<m_pData[nMax]) nIndex = nMax;
			else nIndex = nMax+1;
		}
		else
		{
			if(m_pData[nMin]<tData) nIndex = nMin;
			else if(m_pData[nMax]<tData) nIndex = nMax;
			else nIndex = nMax+1;
		}
	}
	return nIndex;
}
*/
template <class TYPE>
int ML3Array<TYPE>::Find(const TYPE& tData) const
{
	if(m_nSort==0||m_nSize<nML3SearchThreshold)
	{
		for(int i=0;i<m_nSize;i++)
		{
			if(m_pData[i]==tData) return i;
		}
		return -1;
	}
	if(m_pData[0]==tData)
	{
		return 0;
	}
	else if(m_pData[m_nSize-1]==tData)
	{
		return m_nSize-1;
	}
	int nMid, nMin = 0, nMax = m_nSize-1;
	while(nMax-nMin>1)
	{
		nMid = (nMax+nMin)/2;
		if(m_pData[nMid]==tData)
		{
			return nMid;
		}
		if(m_nSort>0)
		{
			if(m_pData[nMid]<tData) nMin = nMid;
			else nMax = nMid;
		}
		else
		{
			if(tData<m_pData[nMid]) nMin = nMid;
			else nMax = nMid;
		}
	}
	return -1;
}

template <class TYPE>
int ML3Array<TYPE>::Insert(const int nIndex, const TYPE& tData)
{
	if(nIndex<0||nIndex>m_nSize)
	{
		return -1;
	}
	m_nSort = 0;
	int nSegmentSize = m_nSize-nIndex;
	SetSize(m_nSize+1,m_nGrowBy);
	if(nSegmentSize>0)
	{
		for(int i=nSegmentSize-1;i>=0;i--)
		{
			m_pData[i+nIndex+1] = m_pData[i+nIndex];
		}
	}
	m_pData[nIndex] = tData;
	return nIndex;
}

template <class TYPE>
int ML3Array<TYPE>::Remove(const int nIndex)
{
	if(nIndex<0||nIndex>=m_nSize)
	{
		return -1;
	}
	int nSegmentSize = m_nSize-nIndex-1;
	if(nSegmentSize>0)
	{
		for(int i=0;i<nSegmentSize;i++)
		{
			m_pData[nIndex+i] = m_pData[nIndex+i+1];
		}
	}
	m_nSize--;
	return m_nSize;
}

template <class TYPE>
int ML3Array<TYPE>::SetAt(const int nIndex, const TYPE& tData)
{
	if(nIndex<0||nIndex>=m_nSize)
	{
		return -1;
	}
	m_nSort = 0;
	m_pData[nIndex] = tData;
	return nIndex;
}

template <class TYPE>
void ML3Array<TYPE>::Push(const TYPE& tData)
{
	m_nSort = 0;
	Add(tData);
}

template <class TYPE>
TYPE ML3Array<TYPE>::Pop()
{
	m_nSize--;
	return m_pData[m_nSize];
}

template <class TYPE>
const TYPE& ML3Array<TYPE>::operator()(const int nIndex) const
{
	return m_pData[nIndex];
}

template <class TYPE>
TYPE& ML3Array<TYPE>::operator[](const int nIndex)
{
	return m_pData[nIndex];
}
/*
template <class TYPE>
void ML3Array<TYPE>::QuickSort(bool (*pSortFunction)(TYPE& , TYPE&))
{
	std::vector<TYPE> v;
	for (int i=0; i<m_nSize; i++)
		v.push_back(m_pData[i]);
	std::sort(v.begin(),v.end(),*pSortFunction);
	for (int i=0; i<m_nSize; i++)
		m_pData[i]=v[i];
}
*/
template <class TYPE>
void ML3Array<TYPE>::Sort(int (*pSortFunction)(TYPE& , TYPE&), bool bDesc/*=false*/)
{
		TYPE tmp;
		int compres;
		for(int i=0; i<m_nSize-1; i++) //bubble sort
			for(int j=i+1;j<m_nSize; j++)
			{
				compres=pSortFunction(m_pData[i],m_pData[j]);
				if((compres<0 && bDesc) || (compres>0 && !bDesc))
				{
					tmp=m_pData[i];
					m_pData[i]=m_pData[j];
					m_pData[j]=tmp;
				}
			}
}
#endif
