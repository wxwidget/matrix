/* this class is the reference of a vector, whoes element is stored in a sequence of memory
 * provide 4 constructors:
 *
 *   1: VectorRef() with empty paramator, which do nothing but construct a object
 *   2: VectorRef(size) with a integer. This construct the object by size, the default value is T()
 *   3: VectorRef(T* data) which a memory point.
 *          Note that, this construtor don't allocate memory at all!
 *          which means, the allocate and deallocate of the memory will done by other Object, rather
 *          than this VectorRef object
 *   4: VectorRef(const VectorRef&v)
 *
 * vector refenece will hold share the same memory, if the Object constructor by 1,2,4
 *
 */
#ifndef ALICE_MATRIX_VECTOR_REF_H
#define ALICE_MATRIX_VECTOR_REF_H

#include <iostream>
#include <assert.h>
#include <inttypes.h>
#include "tracer.h"
namespace alice
{
namespace matrix
{

/**
 * a implement of vector reference
 */
template < typename T, typename TracerType = NoneTrace >
class VectorRef
{
public:
    typedef T ValueType;
private:
    T* mData;
    size_t* mRefCount;

public:
    /** the default constructor
     */
    VectorRef();
    explicit VectorRef(size_t n);
    inline	 VectorRef(T* data);
    inline	 VectorRef(const VectorRef& v);

    ~VectorRef();   /*release the reference counter*/


    inline  T* Data();
    inline  const T* Data() const;

    inline  T& operator[](uint32_t i);
    inline  const T& operator[](uint32_t i) const;

    inline 	size_t RefCount() const;
    inline  bool IsNull() const;
    inline  void Destroy();

    inline  VectorRef<T, TracerType> & operator=(const VectorRef<T, TracerType>& V);
};

template <typename T, typename TracerType>
VectorRef<T, TracerType>::VectorRef() : mData(NULL), mRefCount(NULL)
{
}

template <typename T, typename TracerType>
VectorRef<T, TracerType>::VectorRef(size_t n) : mData(NULL), mRefCount(NULL)
{
    if (n >= 1)
    {
        TracerType::Trace("new data storage.\n");
        mData = new T[n];
        mRefCount = new size_t;
        *mRefCount = 1;
    }
}

template <typename T, typename TracerType>
inline	 VectorRef<T, TracerType>::VectorRef(const VectorRef<T, TracerType>& v): mData(v.mData),
        mRefCount(v.mRefCount)
{
    if (v.mRefCount)
    {
        (*(v.mRefCount))++;
    }
}

//data which is not allocated by VectorRec will not be have a reference count
template <typename T, typename TracerType>
VectorRef<T, TracerType>::VectorRef(T* data) : mData(data), mRefCount(NULL)
{
}

template <typename T, typename R>
inline T* VectorRef<T, R>::Data()
{
    return mData;
}

template <class T, typename R>
inline const T& VectorRef<T, R>::operator[](uint32_t i) const
{
    return mData[i];
}

template <class T, typename R>
inline T& VectorRef<T, R>::operator[](uint32_t i)
{
    return mData[i];
}

template <typename T, typename R>
inline const T* VectorRef<T, R>::Data() const
{
    return mData;
}

template <typename T, typename R>
VectorRef<T, R>& VectorRef<T, R>::operator=(const VectorRef<T, R>& V)
{
    if (mData == V.mData)
    {
        return *this;
    }
    if (mRefCount)
    {
        (*mRefCount)--;
        if ((*mRefCount) == 0)
        {
            Destroy();
        }
    }
    mData = V.mData;
    mRefCount = V.mRefCount;
    if (V.mRefCount)
    {
        (*(V.mRefCount))++;
    }
    return *this;
}

template <typename T, typename R>
void VectorRef<T, R>::Destroy()
{
    if (mRefCount)
    {
        R::Trace("Destorying Counter... \n");
        delete mRefCount;
        if (mData != NULL)
        {
            R::Trace("Destorying Data ... \n");
            delete []mData;
        }
        mData = NULL;
    }
}

template<typename T, typename R>
bool VectorRef<T, R>::IsNull() const
{
    return (mData == NULL ? true : false);
}

/*
*  returns -1 if data is external,
*  returns 0 if a is NULL array,
*  otherwise returns the positive number of vectors sharing
*  		this data space.
*/
template <typename T, typename R>
size_t VectorRef<T, R>::RefCount() const
{
    if (mData == NULL)
    {
        return 0;
    }
    else
    {
        return (mRefCount != NULL ? *mRefCount : -1) ;
    }
}

template <class T, class R>
VectorRef<T, R>::~VectorRef()
{
    if (mRefCount)
    {
        (*mRefCount)--;
        if (*mRefCount == 0)
        {
            Destroy();
        }
    }
}


}//end of namespace matrix
}//end of namespace alice





#endif

