#include "types.hpp"

#include <iostream>
#include <string>

using namespace vaspml;

int main(int /*argc*/, char** /*argv*/)
{
    std::cout << "**********************************************************************************************\n";
    std::cout << "* MultiType (non-shared_ptr)\n";
    std::cout << "**********************************************************************************************\n";
    std::cout << "Number address" << " Number\n";
    MultiType a = 1.23;
    Real& aref = a.get<Real>();
    const Real& acref = a.cget<Real>();
    // a.dget<Real>();  // Should not work, not a shared_ptr!
    // a.dcget<Real>(); // Should not work, not a shared_ptr!
    a = 4.56;
    std::cout << &(std::get<Real>(a)) << " " << std::get<Real>(a) << "\n";
    aref = 7.89;
    // acref = 10.111; // Should not work, const Real&!
    std::cout << &(std::get<Real>(a)) << " " << std::get<Real>(a) << "\n";
    std::cout << &(a.get<Real>()    ) << " " << a.get<Real>()     << "\n";
    std::cout << &(aref             ) << " " << aref              << "\n";
    std::cout << &(acref            ) << " " << acref             << "\n\n";

    std::cout << "**********************************************************************************************\n";
    std::cout << "* const MultiType (non-shared_ptr)\n";
    std::cout << "**********************************************************************************************\n";
    std::cout << "Number address" << " Number\n";
    const MultiType b = 1.23;
    // b.get<Real>(); // Should not work, const MultiType!
    const Real& bcref = b.cget<Real>();
    // b.dget<Real>();  // Should not work, not a shared_ptr!
    // b.dcget<Real>(); // Should not work, not a shared_ptr!
    // b = 4.56; // Should not work, const MultiType!
    // bcref = 7.89; // Should not work, const MultiType!
    std::cout << &(std::get<Real>(b)) << " " << std::get<Real>(b) << "\n";
    std::cout << &(bcref            ) << " " << bcref             << "\n\n";

    std::cout << "**********************************************************************************************\n";
    std::cout << "* MultiType (shared_ptr)\n";
    std::cout << "**********************************************************************************************\n";
    std::cout << "Vector address" << " Data address  " << " Use count\n";
    MultiType v = std::make_shared<Vec1Real>(Vec1Real({1.0, 2.0, 3.0}));
    std::cout << "Setup: .get() vs std::get()\n";
    std::cout << std::get<ShVec1Real>(v) << " " << std::get<ShVec1Real>(v)->data() << " " << std::get<ShVec1Real>(v).use_count() << "\n";
    std::cout << v.get<ShVec1Real>()     << " " << v.get<ShVec1Real>()->data()     << " " << v.get<ShVec1Real>().use_count()     << "\n\n";
    {
        std::cout << "get() -> ShVec\n";
        ShVec1Real vref = v.get<ShVec1Real>();
        std::cout << v.get<ShVec1Real>() << " " << v.get<ShVec1Real>()->data() << " " << v.get<ShVec1Real>().use_count() << "\n";
        std::cout << vref.get()          << " " << vref->data()                << " " << vref.use_count()                << "\n\n";
    }
    {
        std::cout << "get() -> ShVec&\n";
        ShVec1Real& vref = v.get<ShVec1Real>();
        std::cout << v.get<ShVec1Real>() << " " << v.get<ShVec1Real>()->data() << " " << v.get<ShVec1Real>().use_count() << "\n";
        std::cout << vref.get()          << " " << vref->data()                << " " << vref.use_count()                << "\n\n";
    }
    {
        std::cout << "dget() -> Vec&\n";
        Vec1Real& vdref = v.dget<ShVec1Real>();
        std::cout << v.get<ShVec1Real>() << " " << v.get<ShVec1Real>()->data() << " " << v.get<ShVec1Real>().use_count() << "\n";
        std::cout << &(vdref)            << " " << vdref.data()                << "\n\n";
    }
    {
        std::cout << "dget() -> Vec                 CAUTION: makes a copy!\n";
        Vec1Real vdref = v.dget<ShVec1Real>();
        std::cout << v.get<ShVec1Real>() << " " << v.get<ShVec1Real>()->data() << " " << v.get<ShVec1Real>().use_count() << "\n";
        std::cout << &(vdref)            << " " << vdref.data()                << "\n\n";
    }
    {
        std::cout << "cget() -> ShCVec\n";
        ShCVec1Real vcref = v.cget<ShVec1Real>();
        std::cout << v.get<ShVec1Real>() << " " << v.get<ShVec1Real>()->data() << " " << v.get<ShVec1Real>().use_count() << "\n";
        std::cout << vcref.get()         << " " << vcref->data()               << " " << vcref.use_count()               << "\n\n";
    }
    {
        std::cout << "cget() -> const ShCVec&\n";
        const ShCVec1Real& vcref = v.cget<ShVec1Real>();
        std::cout << v.get<ShVec1Real>() << " " << v.get<ShVec1Real>()->data() << " " << v.get<ShVec1Real>().use_count() << "\n";
        std::cout << vcref.get()         << " " << vcref->data()               << " " << vcref.use_count()               << "\n\n";
    }
    {
        std::cout << "dcget() -> const Vec&\n";
        const Vec1Real& vdcref = v.dcget<ShVec1Real>();
        std::cout << v.get<ShVec1Real>() << " " << v.get<ShVec1Real>()->data() << " " << v.get<ShVec1Real>().use_count() << "\n";
        std::cout << &(vdcref)           << " " << vdcref.data()               << "\n\n";
    }
    {
        std::cout << "dcget() -> Vec                CAUTION: makes a copy!\n";
        Vec1Real vdcref = v.dcget<ShVec1Real>();
        std::cout << v.get<ShVec1Real>() << " " << v.get<ShVec1Real>()->data() << " " << v.get<ShVec1Real>().use_count() << "\n";
        std::cout << &(vdcref)           << " " << vdcref.data()               << "\n\n";
    }
    {
        std::cout << "dcget() -> const Vec          CAUTION: makes a copy!\n";
        const Vec1Real vdcref = v.dcget<ShVec1Real>();
        std::cout << v.get<ShVec1Real>() << " " << v.get<ShVec1Real>()->data() << " " << v.get<ShVec1Real>().use_count() << "\n";
        std::cout << &(vdcref)           << " " << vdcref.data()               << "\n\n";
    }

    std::cout << "**********************************************************************************************\n";
    std::cout << "* const MultiType (shared_ptr)\n";
    std::cout << "**********************************************************************************************\n";
    std::cout << "Vector address" << " Data address  " << " Use count\n";
    const MultiType c = std::make_shared<Vec1Real>(Vec1Real({1.0, 2.0, 3.0}));
    std::cout << "Setup: only std::get()\n";
    std::cout << std::get<ShVec1Real>(c) << " " << std::get<ShVec1Real>(c)->data() << " " << std::get<ShVec1Real>(c).use_count() << "\n";
    // std::cout << c.get<ShVec1Real>()     << " " << c.get<ShVec1Real>()->data()     << " " << c.get<ShVec1Real>().use_count()     << "\n\n"; // Should not work, const MultiType!
    //{
    //    std::cout << "get() -> ShVec\n";
    //    ShVec1Real cref = c.get<ShVec1Real>(); // Should not work, const MultiType!
    //}
    //{
    //    std::cout << "get() -> ShVec&\n";
    //    ShVec1Real& cref = c.get<ShVec1Real>(); // Should not work, const MultiType!
    //}
    //{
    //    std::cout << "dget() -> Vec&\n";
    //    Vec1Real& cdref = c.dget<ShVec1Real>(); // Should not work, const MultiType!
    //}
    //{
    //    std::cout << "dget() -> Vec                 CAUTION: makes a copy!\n";
    //    Vec1Real cdref = c.dget<ShVec1Real>(); // Should not work, const MultiType!
    //}
    {
        std::cout << "cget() -> ShCVec\n";
        ShCVec1Real ccref = c.cget<ShVec1Real>();
        std::cout << std::get<ShVec1Real>(c) << " " << std::get<ShVec1Real>(c)->data() << " " << std::get<ShVec1Real>(c).use_count() << "\n";
        std::cout << ccref.get()             << " " << ccref->data()                   << " " << ccref.use_count()                   << "\n\n";
    }
    {
        std::cout << "cget() -> const ShCVec&\n";
        const ShCVec1Real& ccref = c.cget<ShVec1Real>();
        std::cout << std::get<ShVec1Real>(c) << " " << std::get<ShVec1Real>(c)->data() << " " << std::get<ShVec1Real>(c).use_count() << "\n";
        std::cout << ccref.get()             << " " << ccref->data()                   << " " << ccref.use_count()                   << "\n\n";
    }
    {
        std::cout << "dcget() -> const Vec&\n";
        const Vec1Real& cdcref = c.dcget<ShVec1Real>();
        std::cout << std::get<ShVec1Real>(c) << " " << std::get<ShVec1Real>(c)->data() << " " << std::get<ShVec1Real>(c).use_count() << "\n";
        std::cout << &(cdcref)               << " " << cdcref.data()                   << "\n\n";
    }
    {
        std::cout << "dcget() -> Vec                CAUTION: makes a copy!\n";
        Vec1Real cdcref = c.dcget<ShVec1Real>();
        std::cout << std::get<ShVec1Real>(c) << " " << std::get<ShVec1Real>(c)->data() << " " << std::get<ShVec1Real>(c).use_count() << "\n";
        std::cout << &(cdcref)               << " " << cdcref.data()                   << "\n\n";
    }
    {
        std::cout << "dcget() -> const Vec          CAUTION: makes a copy!\n";
        const Vec1Real cdcref = c.dcget<ShVec1Real>();
        std::cout << std::get<ShVec1Real>(c) << " " << std::get<ShVec1Real>(c)->data() << " " << std::get<ShVec1Real>(c).use_count() << "\n";
        std::cout << &(cdcref)               << " " << cdcref.data()                   << "\n\n";
    }

    return 0;
}
