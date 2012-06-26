template<typename T>
class LazyAssign
{
public:
    LazyAssign(T& a):mReal(a)
    {
    }
    inline operator T() 
    {
        return mReal;
    }
    inline operator T() const
    {
        return mReal;
    }
    T& operator = (const T& a)
    {
        mReal = a;
    }
private:
    T& mReal;
};
