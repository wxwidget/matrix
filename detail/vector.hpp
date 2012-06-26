#ifndef ALICE_MATRIX_VECTOR_DETAIL_H
#define ALICE_MATRIX_VECTOR_DETAIL_H
////////////////////////////////////////////////////////////////////////////////
//   constructor
////////////////////////////////////////////////////////////////////////////////
#include <assert.h>

template <typename T>
Vector<T>::Vector() : mDataHolder(), mLength(0), mBegin(NULL) 
{
}

template <typename T>
Vector<T>::Vector(const Vector<T> &A) : mDataHolder(A.mDataHolder), mLength(A.mLength), 
		mBegin(A.mBegin)
{
}

template <typename T>
Vector<T>::Vector(size_t n) : mDataHolder(n), mLength(n), mBegin(mDataHolder.Data())
{
}

template <typename T>
Vector<T>::Vector(size_t n, const T& ini): mDataHolder(n), mLength(n), mBegin(mDataHolder.Data()) 
{
	Fill(mBegin, mBegin + n, ini);
}

template <typename T>
Vector<T>::Vector(size_t n, T *a) : mDataHolder(a), mLength(n) , mBegin(a)
{
}


template <typename T>
inline T& Vector<T>::operator[](uint32_t i) 
{ 
	assert(i >= 0);
	assert(i < mLength);
	return mBegin[i]; 
}

template <typename T>
inline const T& Vector<T>::operator[](uint32_t i) const 
{ 
	assert(i>= 0);
	assert(i < mLength);
	return mBegin[i]; 
}

template <typename T>
inline  T& Vector<T>::operator()(uint32_t i)
{
    return mBegin[i];
}

template <typename T>
inline  const T& Vector<T>::operator()(uint32_t i) const
{
    return this->operator[](i);
}

template <typename T>
Vector<T> & Vector<T>::operator=(const T &a)
{
	Fill(mBegin, mBegin + mLength, a);
	return *this;
}

template <typename T>
Vector<T> Vector<T>::Copy() const
{
	Vector A(mLength);
    memcpy(A.mBegin, mBegin, mLength * sizeof(T)); 
	return A;
}

template <typename T>
Vector<T> & Vector<T>::operator=(const Vector<T> &A)
{
	//auto check mDataHolder.Data() !=  A.mDataHolder.Data()
    mDataHolder = A.mDataHolder; 
    mBegin = A.mBegin;
    mLength = A.mLength;
    return *this;
}

template <typename T>
inline size_t Vector<T>::GetDim() const 
{ 
    return mLength;
} 

template <typename T>
inline   size_t Vector<T>::GetLengh() const
{
    return mLength;
}

template <typename T>
inline bool Vector<T>::IsEmpty() const 
{ 
    return mLength < 1;
} 
template <typename T>
Vector<T>::~Vector() 
{
}


template <typename T>
inline Vector<T> Vector<T>::Sub(uint32_t i0, uint32_t i1)
{
	if (i0 < i1 && i1 <= mLength)
	{
		Vector<T> X(*this);  
		X.mLength = i1-i0;
		X.mBegin += i0;
		return X;
	}
	else
	{
		return Vector<T>();
	}
}

template <typename T>
inline Vector<T> Vector<T>::Prefix(uint32_t end) 
{
    return Sub(0, end);
}

template <typename T>
inline Vector<T> Vector<T>::Suffix(uint32_t begin)
{
    return Sub(begin, mLength);
}

template <typename T>
inline void  Vector<T>::Fill(T* begin, T* end, const T& value)
{
	for (T* p = begin; p < end; ++p)
    {
		*p = value;
    }
}
// iterators
template <typename T>
inline typename Vector<T>::Iterator Vector<T>::Begin()
{
    return mBegin;
}
template <typename T>
inline typename Vector<T>::ConstIterator Vector<T>::Begin() const
{
    return mBegin;
}

template <typename T>
inline typename Vector<T>::Iterator Vector<T>::End()
{
    return mBegin + mLength;
}
template <typename T>
inline typename Vector<T>::ConstIterator Vector<T>::End() const
{
    return mBegin + mLength;
}
template<typename T>
Vector<T>& Vector<T>::operator+= (const T v)
{
    for(T* p = mBegin; p < mBegin + mLength; ++p)
    {
        *p += v;
    }
    return *this;
}
template<typename T>
Vector<T>& Vector<T>::operator-=( const T v)
{
    for(T* p = mBegin; p < mBegin + mLength; ++p)
    {
        *p -= v;
    }
    return *this;
}
template<typename T>
Vector<T>& Vector<T>::operator *= (const T v)
{
    for(T* p = mBegin; p < mBegin + mLength; ++p)
    {
        *p *= v;
    }
    return *this;
}
template<typename T>
Vector<T>& Vector<T>::operator /=( const T v)
{
    for(T* p = mBegin; p < mBegin + mLength; ++p)
    {
        //don't check div zero
        *p /= v;
    }
    return *this;
}
template <typename T>
Vector<T>& Vector<T>::operator+=(const Vector<T>& v)
{
    if (mLength == v.mLength)
    {
        for (int i = 0; i < mLength; ++i)
        {
            mBegin[i] += v[i];
        }
    }
    return *this;
}
template <typename T>
Vector<T>& Vector<T>::operator -= (const Vector<T>& v)
{
    if (mLength == v.mLength)
    {
        for (int i = 0; i < mLength; ++i)
        {
            mBegin[i] -= v[i];
        }
    }
    return *this;
}
template <typename T>
Vector<T>& Vector<T>::operator/=(const Vector<T>& v)
{
    if (mLength == v.mLength)
    {
        for (int i = 0; i < mLength; ++i)
        {
            mBegin[i] /= v[i];
        }
    }
    return *this;
}

template <typename T>
Vector<T>& Vector<T>::operator*=(const Vector<T>& v)
{
    if (mLength == v.mLength)
    {
        for (int i = 0; i < mLength; ++i)
        {
            mBegin[i] *= v[i];
        }
    }
    return *this;
}

template <typename T>
inline void Vector<T>::Set(const T& value)
{
    Fill(mBegin, mBegin + mLength, value);
}

template <typename T>
inline void Vector<T>::Set(uint32_t pos, const T& value)
{
    mBegin[pos] = value;
}

template <typename T>
inline void Vector<T>::Set(uint32_t begin, uint32_t end, const T& value)
{
    Fill(mBegin + begin, mBegin + end, value);
}

template <typename T>
inline void Vector<T>::Set(uint32_t begin, uint32_t end, const T* values, uint32_t start)
{
    if (end < begin) return;
    T* toStart = mBegin + begin;
    const T* fromStart = values + start; 
    for (uint32_t i = 0; i < end - begin; ++i)
    {
        toStart[i] = fromStart[i];
    }
}

#endif
