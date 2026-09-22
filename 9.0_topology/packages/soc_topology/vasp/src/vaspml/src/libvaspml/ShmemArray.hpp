#ifndef SHMEMARRAY_HPP
#define SHMEMARRAY_HPP

#include "MlMPI.hpp"
#include "debug.hpp"
#include "math.hpp"
#include "types.hpp"

#include <algorithm>
#include <cstddef>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>

//debug
#include <iostream>

namespace vaspml
{

/*******************************************************************************************
 * @class ShmemArray
 * @brief Can be used to store arrays in a shared memory either with sysv or MPI.
 *
 * The class can be used to allocate an array in shared memory. To allocate the
 * array in shared memory the user has to supply a shared_ptr to an MlMPI class
 * from which a shared memory communicator can be created.
 * The class can for example be created with \n
 * @code
 * ShmemArray<Real> array( n, mpiIn ); \n
 * @endcode
 * And then the shared memory array can be filled with \n
 * @code
 * for ( std::size_t i = 0; i < n; n++ ) \n
 * if ( array[ "mpiShmem" ].get_rank() == 0 ) array.set_value( i, value ); \n
 * array[ "mpiShmem" ].barrier(); \n
 * @endcode
 * @note when the functionality should be used with system V the code has to
 * be compiled with -Dsysv
 *******************************************************************************************/
template<class T>
class ShmemArray
{

  public:
    /*******************************************************************************************
     * default constructor
     *
     * initializes variables either to zero or nullptr
     *******************************************************************************************/
    ShmemArray(void);
    /*******************************************************************************************
     * constructor
     * @param n number of elements of data array
     * @param mpiShared supply MPI class with shared memory support if shared memory should be used
     *******************************************************************************************/
    ShmemArray(const Int n, const std::shared_ptr<MlMPI>& mpiIn = nullptr);
    /*******************************************************************************************
     * deconstruct the class.
     *
     * deallocates and frees the occupied memory
     *******************************************************************************************/
    ~ShmemArray(void);
    // if one of these are needed check folder
    // /fsc/home/jona/DATA/Documents/work/Redesign/Graveyard
    // there are implementations for the serial part
    /// deleted copy constructor
    ShmemArray(const ShmemArray& other) = delete;
    /// deleted copy assignement constructor
    ShmemArray& operator=(const ShmemArray& other) = delete;
    /// deleted move constructor
    ShmemArray(ShmemArray&& other) noexcept = delete;
    /// deleted move assignement operator
    ShmemArray& operator=(ShmemArray&& other) noexcept = delete;
    /*******************************************************************************************
     * write contents of the Shmem array to the screen.
     *
     * array will be written in three column format. First column is rank and
     * second is index in array and the last column is the value stored
     *******************************************************************************************/
    void print(void);
    /*******************************************************************************************
     * get the total size of the shmem array in number of elements
     *
     * @note since the class is parent class for ShmemArray2D and ShememArray2DVariableLen
     * the  value printed will then be the total number of entries.
     *******************************************************************************************/
    std::size_t get_size(void);
    /*******************************************************************************************
     * set value in array
     * @param n index which sould be changed
     * @param value determines to value entry will be set
     * @param rank determines which rank does the writting
     *
     * @note after function call the user has to call the mpi barrier function.
     * This can be done with arrayName[ "mpiShmem" ].barrier(); arrayName denotes the
     * instance of the created ShmemArray
     *******************************************************************************************/
    void set_value(const std::size_t n, T value, Int rank = 0);
    /*******************************************************************************************
     * return value stored in array
     *
     * @param indx denotes the index of the element which should be returned
     *
     * @note value will be copied
     *******************************************************************************************/
    T get_value(std::size_t indx) const;
    /*******************************************************************************************
     * directly access mpi comunicator used in the class from outside the class.
     *
     * @note these are meant for convienience if a mpi barrier is used
     *******************************************************************************************/
    MlMPI* operator->(void);
    /*******************************************************************************************
     * directly access mpi comunicator used in the class from outside the class.
     *
     * @note these are meant for convienience if a mpi barrier is used.
     * @note MlMPI is returned by reference without const. The user can not alter the
     * MlMPI class from outside by the MlMPI is written. The const is neglected since
     * the mlMPI class has to set it's status info variable reporting if an error had
     * happened
     *******************************************************************************************/
    MlMPI& operator[](const std::string& key);
    /*******************************************************************************************
     * distribute the shared memory segement internode. This function has to be called after
     * shared memory was filled
     *******************************************************************************************/
    void distributeInternode(void);
    /*******************************************************************************************
     * return a reference to the stored memory segment.
     *
     * @warning user can change the memory segment from the outside.
     * With great power comes great responsibility
     * @note this function is meant to be used in combination with cblas routines
     *******************************************************************************************/
    T*& get_dataPointer(void);
    /*******************************************************************************************
     * return a pointer to the stored memory segment.
     *
     * can be used as any ordinary c++ const raw pointer
     *
     * @note a copy of the pointer pointing to the data will be made
     *******************************************************************************************/
    const T* get_dataPointer(void) const;

  protected:
    /*******************************************************************************************
     * allocating the data pointer.
     *
     * The data pointer will be allocated as a standard c++ raw pointer when no MlMPI
     * class was supplied to the constructor. When MlMPI is supplied the data
     * will be either allocated with the system V approach or MPI to shared memory
     *******************************************************************************************/
    void allocate(void);
    /*******************************************************************************************
     * data pointer in which data is stored in either shared memory or c++ raw pointer format
     *******************************************************************************************/
    T* data;
    T* dataOffset;
    /*******************************************************************************************
     * total number of elements of type T stored in shared memory array
     *******************************************************************************************/
    std::size_t size;
    /*******************************************************************************************
     * storing mpi communicator
     *
     * This map stores a shared memory communicator and a inter node communicator.
     * The shared mmemory communicator can be accessed with mpiShmem and the internode
     * communicator can be accessed with key mpiInter.
     *******************************************************************************************/
    std::map<std::string, MlMPI> MPIcomms;
    /*******************************************************************************************
     * determines if an MlMPI class was supplied so that data can be stored to shared memory
     *
     * @note variable will also be set to false when MlMPI was supplied but the code
     * was compiled in a way that no shared memory is used.
     *******************************************************************************************/
    bool useMPI;
};

template<class T>
ShmemArray<T>::ShmemArray(void)
{
    size = 0;
    data = nullptr;
    useMPI = false;
}

template<class T>
ShmemArray<T>::ShmemArray(const Int n, const std::shared_ptr<MlMPI>& mpiIn)
{
    if (n < 0)
    {
        throw std::runtime_error(
            "ERROR:ShmemArray<T>::ShmemArray( const Int n, const std::shared_ptr<MlMPI>& mpiIn )\n"
            "array size smaller zero. This does not make sense\n");
    }
    size = n;
    if (mpiIn == nullptr) { useMPI = false; }
    else
    {
        std::tie(MPIcomms["mpiShmem"], MPIcomms["mpiInter"]) = mpiIn->make_splitSharedInternode();
        useMPI = true;
    }
#ifndef use_shmem
    useMPI = false;
#endif
    allocate();
}

template<class T>
ShmemArray<T>::~ShmemArray(void)
{
    if (useMPI)
    {
        if (MPIcomms["mpiShmem"].get_sharedAllocCheck()) { MPIcomms["mpiShmem"].freeWindow(); }
    }
    else { delete[] data; }
}

template<class T>
std::size_t ShmemArray<T>::get_size(void)
{
    return size;
}

template<class T>
void ShmemArray<T>::set_value(const std::size_t indx, T value, Int rank)
{
    if (!useMPI) { data[indx] = value; }
    else
    {
        if (rank == MPIcomms["mpiShmem"].get_rank()) { data[indx] = value; }
    }
}

template<class T>
T ShmemArray<T>::get_value(std::size_t indx) const
{
    return data[indx];
}

template<class T>
void ShmemArray<T>::print(void)
{
    for (std::size_t i = 0; i < size; i++)
    {
        std::cout << " proc = " << MPIcomms["mpiShmem"].get_rank() << "   " << i << "  " << data[i]
                  << std::endl;
    }
}

template<class T>
void ShmemArray<T>::allocate(void)
{
#ifdef use_shmem
    if (MPIcomms["mpiShmem"].get_sharedCommunicatorCheck()
        || MPIcomms["mpiShmem"].get_numberRanks() > 1)
    {
        // allocate array with shared memory
        MPIcomms["mpiShmem"].allocateShmemSpace(data, size);
    }
    else
    {
        if (MPIcomms["mpiShmem"].get_sharedCommunicatorCheck())
        {
            std::cout << "Note you are allocating a ShmemArray<T> without " << std::endl;
            std::cout << "supplying a shared memory communicator" << std::endl;
            std::cout << "Memory will be allocated on all cores" << std::endl;
        }
        // allocate array with non shared memory or single core
        data = new T[size];
    }
#else
    data = new T[size];
#endif
}

template<class T>
MlMPI* ShmemArray<T>::operator->(void)
{
    return &MPIcomms["mpiShmem"];
}

template<class T>
MlMPI& ShmemArray<T>::operator[](const std::string& key)
{
    VASPML_DEBUG_L1(
        auto it = MPIcomms.find(key);
        if (it == MPIcomms.end())
        {
            throw std::runtime_error(
                "ERROR: MlMPI& ShmemArray<T>operator[]( const std::string& key ) \n"
                "key not found in MPIcomms dictionary. Possible values are \n"
                "mpiShmem or mpiInter");
        }
    );
    return MPIcomms[key];
}

template<class T>
void ShmemArray<T>::distributeInternode(void)
{
    if (useMPI)
    {
        // use inter node communicator to send from root to other
        if (MPIcomms["mpiShmem"].get_rank() == 0) MPIcomms["mpiInter"].bcast(data, size, 0);
        MPIcomms["mpiShmem"].barrier();
    }
}

template<class T>
T*& ShmemArray<T>::get_dataPointer(void)
{
    return data;
}

template<class T>
const T* ShmemArray<T>::get_dataPointer(void) const
{
    return data;
}

/*******************************************************************************************
 * @class ShmemArray2D
 * @brief Can be used to store arrays in a shared memory either with sysv or MPI.
 *
 * The class can be used to allocate an array in shared memory. To allocate the
 * array in shared memory the user has to supply a shared_ptr to an MlMPI class
 * from which a shared memory communicator can be created.
 * The class can for example be created with \n
 * @code
 * ShmemArray<Real> array( n0, n1, mpiIn ); \n
 * @endcode
 * And then the shared memory array can be filled with \n
 * @code
 * for ( std::size_t i = 0; i < n0; n++ )
 * for ( std::size_t j = 0; j < n1; n++ )
 * if ( array[ "mpiShmem" ].get_rank() == 0 ) array.set_value( i,j, value );\n
 * array[ "mpiShmem" ].barrier();\n
 * @endcode
 * @note the faster index of the ShmemArray2D will always be the second index
 * @note when the functionality should be used with system V the code has to
 * be compiled with -Dsysv
 *******************************************************************************************/
template<class T>
class ShmemArray2D : public ShmemArray<T>
{

  public:
    /*******************************************************************************************
     * When using in default constructed way data will not be stored in shared memory
     *******************************************************************************************/
    ShmemArray2D(void);
    /*******************************************************************************************
     * construct the ShmemArray2D
     *
     * @param n0 first dimension of the shared memory array
     * @param n1 second dimension of the shared memory array
     * @mpiIn mpi communicator used create the shared memory communicators
     *******************************************************************************************/
    ShmemArray2D(const Int n0, const Int n1, const std::shared_ptr<MlMPI>& mpiIn = nullptr);
    /*******************************************************************************************
     * set value in array
     * @param n index which sould be changed
     * @param value determines to value entry will be set
     * @param rank determines which rank does the writting
     *
     * @note always call an mpi barrier with the proper communicator present in your class
     * after assigning a value
     *******************************************************************************************/
    void set_value(const std::size_t indx0, const std::size_t indx1, T value, Int rank = 0);
    /*******************************************************************************************
     * return value stored in array
     *
     * @param indx index of element which should be returned from array
     * @note a copy of the value will be made
     *******************************************************************************************/
    T get_value(const std::size_t indx0, const std::size_t indx1) const;
    /*******************************************************************************************
     * getting the first dimension of the shared memory array
     *******************************************************************************************/
    std::size_t get_dimension0(void) const;
    /*******************************************************************************************
     * getting the second dimension of the shared memory array
     *******************************************************************************************/
    std::size_t get_dimension1(void) const;

  private:
    /*******************************************************************************************
     * variable to store the first dimension of the shared memory array
     *
     * total size of the shared memory array will be dimension0 * dimension1
     *******************************************************************************************/
    std::size_t dimension0;
    /*******************************************************************************************
     * variable to store the second dimension of the shared memory array
     *
     * total size of the shared memory array will be dimension0 * dimension1
     *******************************************************************************************/
    std::size_t dimension1;
};

template<class T>
ShmemArray2D<T>::ShmemArray2D(void) : ShmemArray<T>()
{
    dimension0 = 0;
    dimension1 = 0;
}

template<class T>
ShmemArray2D<T>::ShmemArray2D(const Int n0, const Int n1, const std::shared_ptr<MlMPI>& mpiIn) :
    ShmemArray<T>(n0 * n1, mpiIn)
{
    if (n0 < 0)
    {
        throw std::runtime_error(
            "ERROR:ShmemArray2D<T>::ShmemArray2D( const Int n0, const Int n1 \n"
            "                                     const std::shared_ptr<MlMPI>& mpiIn )\n"
            "array size smaller zero. This does not make sense\n");
    }
    if (n1 < 0)
    {
        throw std::runtime_error(
            "ERROR:ShmemArray2D<T>::ShmemArray2D( const Int n0, const Int n1 \n"
            "                                     const std::shared_ptr<MlMPI>& mpiIn )\n"
            "array size smaller zero. This does not make sense\n");
    }
    dimension0 = n0;
    dimension1 = n1;
}

template<class T>
void ShmemArray2D<T>::set_value(const std::size_t indx0, const std::size_t indx1, T value, Int rank)
{
    VASPML_DEBUG_L1(
        if (indx0 >= dimension0)
        {
            throw std::out_of_range("ERROR: ShmemArray2D<T>::set_value( const std::size_t indx0, "
                                    "const std::size_t indx1, T value, Int rank )\n"
                                    "indx0 is out of bounds");
        }
        if (indx1 >= dimension1)
        {
            throw std::out_of_range("ERROR: ShmemArray2D<T>::set_value( const std::size_t indx0, "
                                    "const std::size_t indx1, T value, Int rank )\n"
                                    "indx1 is out of bounds");
        }
    );
    if (!this->useMPI)
    {
        std::size_t indx = indx0 * dimension1 + indx1;
        this->data[indx] = value;
    }
    else
    {
        if (rank == this->MPIcomms["mpiShmem"].get_rank())
        {
            std::size_t indx = indx0 * dimension1 + indx1;
            this->data[indx] = value;
        }
        //this -> MPIcomms[ "mpiShmem" ].barrier();
    }
}

template<class T>
T ShmemArray2D<T>::get_value(std::size_t indx0, std::size_t indx1) const
{
    VASPML_DEBUG_L1(
        if (indx0 >= dimension0)
        {
            throw std::out_of_range("ERROR: ShmemArray2D<T>::get_value( const std::size_t indx0, "
                                    "const std::size_t indx1, T value, Int rank )\n"
                                    "indx0 is out of bounds");
        }
        if (indx1 >= dimension1)
        {
            throw std::out_of_range("ERROR: ShmemArray2D<T>::get_value( const std::size_t indx0, "
                                    "const std::size_t indx1, T value, Int rank )\n"
                                    "indx1 is out of bounds");
        }
    );
    std::size_t indx = indx0 * dimension0 + indx1;
    return this->data[indx];
}

template<class T>
std::size_t ShmemArray2D<T>::get_dimension0(void) const
{
    return dimension0;
}

template<class T>
std::size_t ShmemArray2D<T>::get_dimension1(void) const
{
    return dimension1;
}

/*******************************************************************************************
 * @class ShmemArray2DVariableLen
 * @brief Can be used to store arrays in a shared memory either with sysv or MPI.
 *
 * The class can be used to allocate an array in shared memory. To allocate the
 * array in shared memory the user has to supply a shared_ptr to an MlMPI class
 * from which a shared memory communicator can be created.
 * The class can for example be created with
 * @code
 * ShmemArray<Real> array( n0, n1, mpiIn ); \n
 * @endcode
 * And then the shared memory array can be filled with \n
 * @code
 * for ( std::size_t i = 0; i < n0; n++ )
 * for ( std::size_t j = 0; j < n1[ i ]; n++ )
 * if ( array[ "mpiShmem" ].get_rank() == 0 ) array.set_value( i,j, value ); \n
 * array[ "mpiShmem" ].barrier(); \n
 * @endcode
 * @note the faster index of the ShmemArray2DVariableLen will always be the second index
 * @note when the functionality should be used with system V the code has to
 * be compiled with -Dsysv
 *******************************************************************************************/

template<class T>
class ShmemArray2DVariableLen : public ShmemArray<T>
{

  public:
    /*******************************************************************************************
     * When using in default constructed way data will not be stored in shared memory
     *******************************************************************************************/
    ShmemArray2DVariableLen(void);
    /*******************************************************************************************
     * constructing a ShmemArray2DVariableLen
     *
     * @param lengths the lengths of the second dimensions of the ShmemArray. The size of lengths
     *determines the first dimension
     * @param mpiIn mpi class to make shared memory access and internode communication possible
     *
     *******************************************************************************************/
    ShmemArray2DVariableLen(const Vec1Int& lengths, const std::shared_ptr<MlMPI>& mpiIn = nullptr);
    /*******************************************************************************************
     * set single value in array
     *
     * @param indx0 index which should be changed ( slow index )
     * @param indx1 index which should be changed ( fast index )
     * @param value determines to value entry will be set
     * @param shmemRank determines which mpiShmem rank does the writing
     * @param interRank determines which mpiInter rank does the writing
     *
     * @note only the cores for which both shmemRank and interRank are satisfied
     * will do the writing
     *******************************************************************************************/
    void set_value(const std::size_t indx0,
                   const std::size_t indx1,
                   const T           value,
                   const Int         rank = 0);
    /*******************************************************************************************
     * copy data array to the shared memory array
     *
     * @param inputData data array which will be copied to the shared memory array
     * @param shmemRank determines which mpiShmem rank does the writing
     * @param interRank determines which mpiInter rank does the writing
     *
     * @note only the cores for which both shmemRank and interRank are satisfied
     * will do the writing
     *******************************************************************************************/
    void set_value(const std::vector<std::vector<T>>& inputData,
                   const Int                          shemRank = 0,
                   const Int                          interRank = 0);
    /*******************************************************************************************
     * return value stored in array
     *
     * @param indx0 entry which should be returned from array ( slow index )
     * @param indx1 entry which should be returned from array ( fast index )
     *******************************************************************************************/
    T get_value(const std::size_t indx0, const std::size_t indx1) const;
    /*******************************************************************************************
     * get a const pointer to a slice in the second dimension of ShmemArray
     *
     * @param indx0 index of the first dimension for which the slice should be obtained
     *
     * returns a pointer to a slice chosen by indx0
     * as data[ indx0,: ],
     * @note that a pointer has to be created
     *******************************************************************************************/
    const T* get_slice(const std::size_t indx0) const;
    /**
     * get a pointer to a slice in the second dimension of ShmemArray
     *
     * @param indx0 index of the first dimension for which the slice should be obtained
     *
     * returns a pointer to a slice chosen by indx0
     * as data[ indx0,: ], note that a pointer has to be created
     *******************************************************************************************/
    T*  get_slice(const std::size_t indx0);
    T*& get_sliceReference(const std::size_t indx0);
    /*******************************************************************************************
     * returns the number of elements stored in a certain slice
     *
     * @param indx0 index for the first entry of the array
     *******************************************************************************************/
    std::size_t get_lengthSlice(const std::size_t indx0) const;
    /*******************************************************************************************
     * get size of zeroth dimension of the shared memory array
     *******************************************************************************************/
    std::size_t get_size0(void) const;

  private:
    /*******************************************************************************************
     * length of first dimension of ShmemArray
     *******************************************************************************************/
    std::size_t dimension0;
    /*******************************************************************************************
     * lengths of the second dimensions.
     *
     * lengths.size() matches dimension0
     *******************************************************************************************/
    std::vector<std::size_t> lengths;
    /*******************************************************************************************
     * offsets which tell at which entries of the 1D contiguous shared memory array the slow index
     * is increased
     *******************************************************************************************/
    std::vector<std::size_t> offsets;
};

template<class T>
std::size_t ShmemArray2DVariableLen<T>::get_size0(void) const
{
    return dimension0;
}

template<class T>
ShmemArray2DVariableLen<T>::ShmemArray2DVariableLen(void) : ShmemArray<T>()
{
    dimension0 = 0;
    lengths.resize(0);
    offsets.resize(0);
}

template<class T>
ShmemArray2DVariableLen<T>::ShmemArray2DVariableLen(const Vec1Int&                lengths,
                                                    const std::shared_ptr<MlMPI>& mpiIn) :
    ShmemArray<T>(math::sumVector(lengths), mpiIn)
{
    VASPML_DEBUG_L1(
        if (!std::all_of(lengths.begin(), lengths.end(), [](Int i) { return i > 0; }))
        {
            throw std::runtime_error("ERROR:ShmemArray2DVariableLen<T>::ShmemArray2DVariableLen: "
                                     "Some of your vector dimensions are below zero");
        }
    );

    dimension0 = lengths.size();
    offsets.resize(dimension0);
    this->lengths.resize(dimension0);
    this->lengths[0] = lengths[0];
    for (std::size_t i = 1; i < lengths.size(); i++)
    {
        offsets[i] = offsets[i - 1] + lengths[i - 1];
        this->lengths[i] = lengths[i];
    }
}

template<class T>
void ShmemArray2DVariableLen<T>::set_value(const std::size_t indx0,
                                           const std::size_t indx1,
                                           const T           value,
                                           const Int         rank)
{
    VASPML_DEBUG_L1(
        if (indx0 >= dimension0)
        {
            throw std::out_of_range(
                "ERROR: ShmemArray2DVariableLen<T>::set_value( const std::size_t indx0, const "
                "std::size_t indx1, const T value, const Int rank )\n"
                "indx0 is out of bounds");
        }
        if (indx1 >= lengths[indx0])
        {
            throw std::out_of_range(
                "ERROR: ShmemArray2DVariableLen<T>::set_value( const std::size_t indx0, const "
                "std::size_t indx1, const T value, const Int rank )\n"
                "indx1 is out of bounds");
        }
    );
    std::size_t indx = offsets[indx0] + indx1;
    if (!this->useMPI) { this->data[indx] = value; }
    else
    {
        if (rank == this->MPIcomms["mpiShmem"].get_rank()) { this->data[indx] = value; }
        //this->MPIcomms[ "mpiShmem" ].barrier();
    }
}

template<class T>
void ShmemArray2DVariableLen<T>::set_value(const std::vector<std::vector<T>>& inputData,
                                           [[maybe_unused]] const Int         shmemRank,
                                           [[maybe_unused]] const Int         interRank)
{

    VASPML_DEBUG_L1(
        if (this->MPIcomms["mpiShmem"].get_rank() == shmemRank
            && this->MPIcomms["mpiInter"].get_rank() == interRank)
        {

            std::size_t inputSize =
                std::accumulate(inputData.begin(),
                                inputData.end(),
                                0,
                                [](std::size_t sum, const std::vector<T>& inner_vec)
                                { return sum + inner_vec.size(); });

            if (inputSize != this->size)
            {
                throw std::out_of_range("ERROR: ShmemArray2DVariableLen<T>::set_value( const "
                                        "std::vector<std::vector<T>>& inputData, const Int "
                                        "shmemRank, const Int interRank )\n"
                                        "Size of inputData and number of elements in "
                                        "ShmemArray2DVariableLen does not match\n");
            }
        }
    );

    if (!this->useMPI)
    {
        std::size_t counter = 0;
        for (std::size_t indx0 = 0; indx0 < dimension0; indx0++)
        {
            for (std::size_t indx1 = 0; indx1 < lengths[indx0]; indx1++)
            {
                this->data[counter] = inputData[indx0][indx1];
                counter++;
            }
        }
    }
    else
    {
#ifdef use_shmem
        if (this->MPIcomms["mpiShmem"].get_rank() == shmemRank
            and this->MPIcomms["mpiInter"].get_rank() == interRank)
        {
#endif
            std::size_t counter = 0;
            for (std::size_t indx0 = 0; indx0 < dimension0; indx0++)
            {
                for (std::size_t indx1 = 0; indx1 < lengths[indx0]; indx1++)
                {
                    this->data[counter] = inputData[indx0][indx1];
                    counter++;
                }
            }
#ifdef use_shmem
        }
        this->MPIcomms["mpiShmem"].barrier();
        this->distributeInternode();
#endif
    }
}

template<class T>
T ShmemArray2DVariableLen<T>::get_value(std::size_t indx0, std::size_t indx1) const
{
    VASPML_DEBUG_L1(
        if (indx0 >= dimension0)
        {
            throw std::out_of_range(
                "ERROR: ShmemArray2DVariableLen<T>::get_value( const std::size_t indx0, const "
                "std::size_t indx1, T value, Int rank )\n"
                "indx0 is out of bounds");
        }
        if (indx1 > lengths[indx0])
        {
            throw std::out_of_range(
                "ERROR: ShmemArray2DVariableLen<T>::get_value( const std::size_t indx0, const "
                "std::size_t indx1, T value, Int rank )\n"
                "indx1 is out of bounds");
        }
    );
    std::size_t indx = offsets[indx0] + indx1;
    return this->data[indx];
}

template<class T>
const T* ShmemArray2DVariableLen<T>::get_slice(const std::size_t indx0) const
{
    VASPML_DEBUG_L1(
        if (indx0 >= dimension0)
        {
            throw std::out_of_range("ERROR: const T*& ShmemArray2DVariableLen<T>::get_slice( const "
                                    "std::size_t indx0 ) const\n"
                                    "indx0 is out of bounds");
        }
    );
    return this->data + offsets[indx0];
}

template<class T>
T* ShmemArray2DVariableLen<T>::get_slice(const std::size_t indx0)
{
    VASPML_DEBUG_L1(
        if (indx0 >= dimension0)
        {
            throw std::out_of_range(
                "ERROR: T* ShmemArray2DVariableLen<T>::get_slice( const std::size_t indx0 )\n"
                "indx0 is out of bounds");
        }
    );
    return this->data + offsets[indx0];
}

template<class T>
T*& ShmemArray2DVariableLen<T>::get_sliceReference(const std::size_t indx0)
{
    VASPML_DEBUG_L1(
        if (indx0 >= dimension0)
        {
            throw std::out_of_range(
                "ERROR: T* ShmemArray2DVariableLen<T>::get_slice( const std::size_t indx0 )\n"
                "indx0 is out of bounds");
        }
    );
    this->dataOffset = this->data + offsets[indx0];
    return this->dataOffset;
}

template<class T>
std::size_t ShmemArray2DVariableLen<T>::get_lengthSlice(const std::size_t indx0) const
{
    VASPML_DEBUG_L1(
        if (indx0 >= dimension0)
        {
            throw std::out_of_range(
                "ERROR: std::size_t ShmemArray2DVariableLen<T>::get_lengthSlice( const std::size_t "
                "indx0 ) const \n"
                "indx0 is out of bounds");
        }
    );
    return lengths[indx0];
}

} //namespace vaspml

#endif
