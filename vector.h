#ifndef ALICE_MATRIX_VECTOR_H
#define ALICE_MATRIX_VECTOR_H
#include <iostream>
#include <inttypes.h>
#include "tracer.h"
#include "vector_ref.h"
namespace alice
{
namespace matrix
{

template <typename T>
class Vector
{
public:
    typedef T* Iterator;
    typedef const T* ConstIterator;

public:
    //
    Vector();
    explicit Vector(size_t n);
    /**
     * @brief create Vector, and init with ini
     *
     * @param count the lengh of the array
     * @param ini the initial value
     */
    Vector(size_t count, const T& ini);
    /**
     * @brief create a Vector, yet the array is allocated not by the clas
     *
     * @param count the length of the array
     * @param mem the memory which is already allocated, which is not controled by reference counter
     */
    Vector(size_t count,  T* mem);
    Vector(const Vector &lhs);

    ~Vector();

    //assignment
    inline   Vector & operator=(const T &a);
    inline   Vector & operator=(const Vector &lhs);

    //accessors
    inline   T& operator[](uint32_t i);
    inline   const T& operator[](uint32_t i) const;
    inline   T& operator()(uint32_t i);
    inline   const T& operator()(uint32_t i) const;

    inline Vector<T> Sub(uint32_t begin, uint32_t end);
    inline Vector<T> Prefix(uint32_t end);
    inline Vector<T> Suffix(uint32_t begin);
    
    /** Set all element
     *
     * @param value   The value
     */
    inline void Set(const T& value);
    /** Set single element
     *
     * @param pos  Index
     * @param value a[i]
     */
    inline void Set(uint32_t pos, const T& value);
    /** Set some element to a value
     *
     * @param begin     The index of the first element
     * @param end       The index before the end element
     * @param value     The value   
     */
    inline void Set(uint32_t begin, uint32_t end, const T& value);
    /** set some element to value from start index of value
     * 
     * @param begin     The index of the first element
     * @parma end       Then index before the end element
     * @param value     The the start point copying from
     * @param start     copy values[start..end-begin+start] to self
     */
    inline void Set(uint32_t begin, uint32_t end, const T* values, uint32_t start);

    //iterator
    inline Iterator Begin();
    inline ConstIterator Begin() const;
    inline Iterator End();
    inline ConstIterator End() const;
    //query
    const VectorRef<T>&  GetReference() const { return mDataHolder;}
    const T* GetArray() const { mDataHolder.Data();}
    inline   size_t GetDim() const;
    inline   size_t GetLengh() const;
    inline   size_t GetSize() const { return mLength;}
    inline   bool IsEmpty() const;
    inline   operator bool() const { return !IsEmpty();}
    //ops
    inline Vector<T> Copy() const;
    inline Vector<T> Clone() const
    {
        return Copy();
    }
    //basic operation
    inline Vector<T>& operator+=(const T);
    inline Vector<T>& operator-=(const T);
    inline Vector<T>& operator*=(const T);
    inline Vector<T>& operator/=(const T);

    inline Vector<T>& operator+=(const Vector<T>&);
    inline Vector<T>& operator-=(const Vector<T>&);
    inline Vector<T>& operator*=(const Vector<T>&);
    inline Vector<T>& operator/=(const Vector<T>&);
private:
    inline  void Fill(T* begin, T* end, const T& value);
    VectorRef<T> mDataHolder;
    size_t mLength;
    T* mBegin;

};
#include "detail/vector.hpp"
}
}
#include "vector_utils.h"
#endif
