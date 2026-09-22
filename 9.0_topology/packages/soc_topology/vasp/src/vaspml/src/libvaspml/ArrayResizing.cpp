#include "ArrayResizing.hpp"

using namespace vaspml;

ArrayResizing1D::ArrayResizing1D(void) : actDim(0), maxDim(0)
{}

ArrayResizing1D::ArrayResizing1D(const Size_t& n) : actDim(n), maxDim(n)
{}

void ArrayResizing1D::init(const Size_t& n)
{
    actDim = n;
    maxDim = n;
}

ArrayResizing2D::ArrayResizing2D(void) : act1Dim(0), max1Dim(0), actSizeRec(0), maxSizeRec(0)
{}

ArrayResizing2D::ArrayResizing2D(const std::size_t& n) :
    act1Dim(n),
    max1Dim(n),
    actSize(n),
    maxSize(n)
{}
