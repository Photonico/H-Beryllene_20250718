#ifndef MPI_HPP
#define MPI_HPP

#include "debug.hpp"
#include "types.hpp"

#include <cassert>
#include <functional>
#include <iostream>
#include <tuple>
#include <vector>
// Avoid annoying compiler warnings from OpenMPI's deprecated C++ bindings.
#define OMPI_SKIP_MPICXX 1
#include <cstddef>
#include <mpi.h>
#include <stdexcept>
#include <tuple>

#ifdef sysv
#include <sys/ipc.h>
#include <sys/shm.h>
#endif

namespace vaspml
{

/** determine what MPI_Datatype is used in the template functions
 * used in this header file
 *@param template parameter to be determined
 *@returns MPI_Datatype or throws an assert error
 */
template<typename T>
[[nodiscard]] constexpr MPI_Datatype M_getMPIType(void)
{
    MPI_Datatype mpiType = MPI_DATATYPE_NULL;
    if constexpr (std::is_same_v<T, Real>) { mpiType = MPI_DOUBLE; }
    else if constexpr (std::is_same_v<T, Int>) { mpiType = MPI_INT; }
    else if constexpr (std::is_same_v<T, Int64>) { mpiType = MPI_INT64_T; }
    else if constexpr (std::is_same_v<T, UInt>) { mpiType = MPI_UNSIGNED; }
    else if constexpr (std::is_same_v<T, UInt64>) { mpiType = MPI_UINT64_T; }

    assert((mpiType != MPI_DATATYPE_NULL) && "Error in M_getMPIType; No valid MPI data-type found");

    return mpiType;
}

/*******************************************************************************************
 * @class MlMPI
 * @brief manges MPI instances. Gives ability to broadcast data, make reductions
 *
 * The MlMPI class implements MPI functionality into VASPml. The functionality
 * includes creating new MPI instances. Create shared memory communicators or
 * internode communicators from some given communicator. The class implements
 * MPI data broadcast functionality as member functions like reduce, allreduce, bcast...
 *
 * Usage example: \n
 * create MlMPI world communicator \n
 * <code>
 * MlMPI mpiWorld; \n
 * mpiWorld.make_WorldComm(); \n
 * </code>
 * then let's create a shared memory communicator \n
 * <code>
 * MlMPI sharedMPI = mpiWorld.make_splitShared(); \n
 * </code>
 *******************************************************************************************/
class MlMPI
{

  public:
    /**
     * default constructor
     * initialize
     */
    MlMPI(void);
    /*******************************************************************************************
     * creating a new MPI class
     *
     * @param communicator MPI_Comm from which class will be constructed
     * @param sharedCommunicator determines if the supplied communicator supports shred memory
     *******************************************************************************************/
    MlMPI(MPI_Comm communicator, bool sharedCommunicator = false);
    /*******************************************************************************************
     * desctructing MPI class.
     *
     * deallocates shared memory window if it was used
     * call MPI_Finalize if class was prepared
     * in MPI comm world state
     *******************************************************************************************/
    ~MlMPI(void);
    /*******************************************************************************************
    *******************************************************************************************/
    MlMPI(const MlMPI& other);
    /*******************************************************************************************
    *******************************************************************************************/
    MlMPI& operator=(const MlMPI& other);
    /*******************************************************************************************
     * deleted move constructor
     *******************************************************************************************/
    MlMPI(MlMPI&& other) noexcept;
    /*******************************************************************************************
     * deleting move assignement operator
     *******************************************************************************************/
    MlMPI& operator=(MlMPI&& other) noexcept;
    /**
     * split a given communicator into shared
     * sub communicators
     */
    MlMPI make_splitShared(void) const;
    /**
     * split a given communicator into internode
     * sub communicators
     */
    std::tuple<MlMPI, MlMPI> make_splitSharedInternode(void) const;
    /**
     * create the MPI world communicator from
     * the deafult initialized class
     *
     * note only one MPI world communicator is allowed
     * per program. If user asks for a second
     * MPI world, routine will throw an error
     * and terminate
     */
    void make_WorldComm(void);
    /**
     * free a MPI shared memory window
     *
     * routine will check if there was shared memory
     * window allocated. Otherwise does nothing
     */
    void freeWindow(void);
    /**
     * allocate a shared memory adress
     *
     * @param data is the data pointer of type T to allocate
     * @param size is the size of the data to allocate in units of syzeof(T)
     *        in other words the number of elemnts to allocate for type T.
     * routine allocates shared memory adress either
     * with the MPI functionality or by using the System V
     * library
     */
    template<typename T>
    void allocateShmemSpace(T*& data, std::size_t size);
    //~*******************************************************************************************
    // data distribution
    //~*******************************************************************************************
    /**
     * broadcast vector data from root node to all others
     *@param data  data to distribute
     *@param count size of the send buffer
     *@param root  node from which to send to others
     */
    template<typename T>
    void bcast(T*& data, Int count, Int root);
    /**
     * broadcast scalar data from root node to all others
     *@param data data to distribute
     *@param root node from which to send to others
     */
    template<typename T>
    void bcast(T& data, int root);
    /**
     * scatter data vector from root node to all others and distribute it
     *@param data_send   data to distribute
     *@param data_recive data buffer where to recive
     *@param size_send   size of data chunk to send
     *@param size of data chunk to recive
     *@param root node from which to send to others
     */
    template<typename T, typename U>
    void scatter(T*& data_send, U*& data_recive, Int size_send, Int size_recive, Int root);
    /**
     * collect data from communicator to node root
     *@param data_send   data to send to root
     *@param data_recive buffer where data is recived
     *@param size_send   size of send buffer
     *@param size_recive size of recive buffer
     *@param root        node which collects data
     */
    template<typename T, typename U>
    void gather(T*& data_send, U*& data_recive, Int size_send, Int size_recive, Int root);
    /**
     * collect data from communicator to node root
     *@param data_send   data to send to root
     *@param data_recive buffer where data is recived
     *@param size_send   size of send buffer
     *@param size_recive size of recive buffer
     */
    template<typename T, typename U>
    void allGather(T*& data_send, U*& data_recive, Int size_send, Int size_recive);
    //~*******************************************************************************************
    // end data distribution
    //~*******************************************************************************************
    //~*******************************************************************************************
    // reduce operations
    //~*******************************************************************************************
    /**reduce sum for scalars
     *@param data_in  data to reduce
     *@param data_out memory were reduction is stored to
     *@param reciver  node on wich data is recived
     */
    template<typename T>
    void reduceSum(T& data_in, T& data_out, Int receiver);
    /**reduce sum for vectors
     *@param data_in  data to reduce
     *@param data_out memory were reduction is stored to
     *@param size     number of elements in vector
     *@param reciver  node on wich data is recived
     */
    template<typename T>
    void reduceSum(T*& data_in, T*& data_out, Int size, Int receiver);
    /**reduce sum for scalars
     *@param data_in  data to reduce
     *@param data_out memory were reduction is stored to
     *@param reciver  node on wich data is recived
     */
    template<typename T>
    void reduceProduct(T& data_in, T& data_out, Int receiver);
    /**reduce product for vectors
     *@param data_in  data to reduce
     *@param data_out memory were reduction is stored to
     *@param size     number of elements in vector
     *@param reciver  node on wich data is recived
     */
    template<typename T>
    void reduceProduct(T*& data_in, T*& data_out, Int size, Int receiver);
    /**reduce max for scalars
     *@param data_in  data to reduce
     *@param data_out memory were reduction is stored to
     *@param reciver  node on wich data is recived
     */
    template<typename T>
    void reduceMax(T& data_in, T& data_out, Int receiver);

    /**reduce max for vectors
     *@param data_in  data to reduce
     *@param data_out memory were reduction is stored to
     *@param size     number of elements in vector
     *@param reciver  node on wich data is recived
     */
    template<typename T>
    void reduceMax(T*& data_in, T*& data_out, Int size, Int receiver);
    /**reduce min for scalars
     *@param data_in  data to reduce
     *@param data_out memory were reduction is stored to
     *@param reciver  node on wich data is recived
     */
    template<typename T>
    void reduceMin(T& data_in, T& data_out, Int receiver);
    /**reduce min for vectors
     *@param data_in  data to reduce
     *@param data_out memory were reduction is stored to
     *@param size     number of elements in vector
     *@param reciver  node on wich data is recived
     */
    template<typename T>
    void reduceMin(T*& data_in, T*& data_out, Int size, Int receiver);
    //~*******************************************************************************************
    // end reduce operations
    //~*******************************************************************************************
    //~*******************************************************************************************
    // all reduce operations
    //~*******************************************************************************************
    /**all reduce sum for scalars
     *
     *all nodes have data in the end
     *@param data_in  data to reduce
     *@param data_out memory were reduction is stored to
     */
    template<typename T>
    void allReduceSum(T& data_send, T& data_recv);
    /**all reduce sum for vectors
     *
     *all nodes have data in the end
     *@param data_in  data to reduce
     *@param data_out memory were reduction is stored to
     *@param size     number elements in vector
     */
    template<typename T>
    void allReduceSum(T*& data_send, T*& data_recv, Int size);
    /**all reduce product for scalars
     *
     *all nodes have data in the end
     *@param data_in  data to reduce
     *@param data_out memory were reduction is stored to
     */
    template<typename T>
    void allReduceProduct(T& data_send, T& data_recv);
    /**all reduce product for vectors
     *
     *all nodes have data in the end
     *@param data_in  data to reduce
     *@param data_out memory were reduction is stored to
     *@param size     number elements in vector
     */
    template<typename T>
    void allReduceProduct(T*& data_send, T*& data_recv, Int size);
    /**all reduce max for scalars
     *
     *all nodes have data in the end
     *@param data_in  data to reduce
     *@param data_out memory were reduction is stored to
     */
    template<typename T>
    void allReduceMax(T& data_send, T& data_recv);
    /**all reduce max for vectors
     *
     *all nodes have data in the end
     *@param data_in  data to reduce
     *@param data_out memory were reduction is stored to
     *@param size     number elements in vector
     */
    template<typename T>
    void allReduceMax(T*& data_send, T*& data_recv, Int size);
    /** all reduce min for scalars
     *
     *  all nodes have data in the end
     *@param data_in  data to reduce
     *@param data_out memory were reduction is stored to
     */
    template<typename T>
    void allReduceMin(T& data_send, T& data_recv);
    /** all reduce min function for vectors
     *
     * all nodes have data in the end
     *@param data_in  ->  data to reduce
     *@param data_out ->  memory were reduction is stored to
     *@param size     ->  number elements in vector
     */
    template<typename T>
    void allReduceMin(T*& data_send, T*& data_recv, Int size);
    //~*******************************************************************************************
    // end all reduce operations
    //~*******************************************************************************************
    /// get rank of MPI process
    const int& get_rank(void) const;
    /// get total number of ranks
    const int& get_numberRanks(void) const;
    /*******************************************************************************************
     * get start and end index for block distribution of n elements
     * @param total number of elements
     * returns start and end index as a tuple
     *******************************************************************************************/
    std::tuple<int, int> blockLoop(std::size_t n);
    /// return control varible to check if shared memory was allocated
    bool get_sharedAllocCheck(void) const;
    /// return control variable if communicator is a shared memory communcator
    bool get_sharedCommunicatorCheck(void) const;
    /// check if there was a communicator initialized
    bool get_noCommunicatorCheck(void) const;
    //~*******************************************************************************************
    // synchronization and barriers
    //~*******************************************************************************************
    /**
     * making a MPI_Fence call and the given
     * communicator
     *
     * Should be used when working with MPI_Put
     * and MPI_Get shared memory
     */
    void fence(void);
    /**
     * makes a MPI_Barrier
     *
     * all processes have to wait together here
     */
    void barrier(void);
    //~*******************************************************************************************
    // end synchronization and barriers
    //~*******************************************************************************************
    /*******************************************************************************************
     * return the MPI communicator
     *
     * This routine is meant to be used in combination with the ScaLapack class to couple
     * the scalapack context to an MPI communicator. If no communicator was set the function
     * returns MPI_COMM_NULL
     *
     * @return communicator MPI communicator of class
     *******************************************************************************************/
    MPI_Comm get_communicator(void) const;

  private:
    /*******************************************************************************************
     * check how much shared memory is left. If there is not enough shared memory error is thrown
     *
     * This function checks over all instances of allocated shared memory if there is still
     * enough shared memory left to do the next allocation. If there is not enough
     * shared memory left the function will throw an error and terminate the code.
     *
     * @param sizeType size of the data type to analyse in bytes per element
     * @param numberElements the number of elements of the data type to allocate
     *
     * the product of sizeType * numberElements gives the total size in bytes to allocate
     *******************************************************************************************/
    void checkEnoughSharedMemory(const std::size_t sizeType, const std::size_t numberElements);
    /*******************************************************************************************
     * read the maximal size of available shared memory to the global variable availableShmem
     *******************************************************************************************/
    void readMaximalAvailableSharedMemory(void);
    /// set number of ranks by calling MPI_Comm_size
    void set_numberRanks(void);
    // set rank id by calling MPI_Comm_rank
    void set_rank();
    /// info determines if error occured during MPI calls
    Int info;
    /// stores number of ranks in communicator
    Int numberRanks;
    /// stores rank of process withing the communicator
    Int rank;
    /// storing the MPI communicator
    MPI_Comm communicator;
    /**
     * MPI shared memory window, only used when shared memory is allocated
     *
     * otherwise set to MPI_WIN_NULL
     */
    MPI_Win win;
    /// controll variable if shared memory was allocated
    bool sharedAlloc;
    /// controll variabe if communicator supports shared memory allocation
    bool sharedCommunicator;
    /// checks if a communictor was set. Otherwise communicator would be MPI_COMM_NULL
    bool noCommunicator;
    /// check if file was read which contains the amount of available shared memory
    static bool shmemFileRead;
    /// available shared memory size in bytes
    static unsigned long long int availableShmem;
};

template<typename T>
void MlMPI::bcast(T*& data, Int count, Int root)
{
    info = MPI_Bcast(data, count, M_getMPIType<T>(), root, communicator);
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error(
                "ERROR: MlMP::bcast( T*& data, Int count, Int root ): MPI_Bcast()");
        }
    );
}

template<typename T>
void MlMPI::bcast(T& data, int root)
{
    info = MPI_Bcast(&data, 1, M_getMPIType<T>(), root, communicator);
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error("ERROR: MlMP::bcast( T& data, int root ): MPI_bcast()");
        }
    );
}

template<typename T, typename U>
void MlMPI::scatter(T*& data_send, U*& data_recive, Int size_send, Int size_recive, Int root)
{
    if (&data_send != &data_recive)
    {
        info = MPI_Scatter(data_send,
                           size_send,
                           M_getMPIType<T>(),
                           data_recive,
                           size_recive,
                           M_getMPIType<U>(),
                           root,
                           communicator);
    }
    else
    {
        if (rank == root)
        {
            info = MPI_Scatter(data_send,
                               size_send,
                               M_getMPIType<T>(),
                               MPI_IN_PLACE,
                               size_recive,
                               M_getMPIType<U>(),
                               root,
                               communicator);
        }
        else
        {
            info = MPI_Scatter(NULL,
                               size_send,
                               M_getMPIType<T>(),
                               data_send,
                               size_recive,
                               M_getMPIType<U>(),
                               root,
                               communicator);
        }
    }
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error(
                "ERROR: MlMPI::scatter( T*& data_send, U*& data_recive,\n"
                "                       Int size_send, Int size_recive, Int root ):\n"
                "MPI_Scatter()");
        }
    );
}

template<typename T, typename U>
void MlMPI::gather(T*& data_send, U*& data_recive, Int size_send, Int size_recive, Int root)
{
    info = MPI_Gather(data_send,
                      size_send,
                      M_getMPIType<T>(),
                      data_recive,
                      size_recive,
                      M_getMPIType<U>(),
                      root,
                      communicator);
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error(
                "ERROR: MlMPI::gather(T*& data_send, U*& data_recive,\n"
                "                     Int size_send, Int size_recive, Int root \n"
                "MPI_Gather()");
        }
    );
}

template<typename T, typename U>
void MlMPI::allGather(T*& data_send, U*& data_recive, Int size_send, Int size_recive)
{
    info = MPI_Allgather(data_send,
                         size_send,
                         M_getMPIType<T>(),
                         data_recive,
                         size_recive,
                         M_getMPIType<U>(),
                         communicator);
    if (info != 0)
    {
        throw std::runtime_error("ERROR: MlMPI::allGather( T*& data_send, U*& data_recive,\n"
                                 "                         Int size_send, Int size_recive\n"
                                 "MPI_Allgather()");
    }
}

template<typename T>
void MlMPI::reduceSum(T& data_in, T& data_out, Int receiver)
{
    // differing send and recive buffer
    if (&data_in != &data_out)
    {
        info =
            MPI_Reduce(&data_in, &data_out, 1, M_getMPIType<T>(), MPI_SUM, receiver, communicator);
    }
    // send and recive buffer are the same
    else
    {
        if (this->get_rank() == receiver)
        {
            info = MPI_Reduce(MPI_IN_PLACE,
                              &data_in,
                              1,
                              M_getMPIType<T>(),
                              MPI_SUM,
                              receiver,
                              communicator);
        }
        else
        {
            info = MPI_Reduce(&data_in,
                              &data_out,
                              1,
                              M_getMPIType<T>(),
                              MPI_SUM,
                              receiver,
                              communicator);
        }
    }
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error(
                "ERROR: MlMPI::reduceSum( T& data_in, T& data_out, Int receiver )\n"
                "MPI_Reduce()");
        }
    );
}

template<typename T>
void MlMPI::reduceSum(T*& data_in, T*& data_out, Int size, Int receiver)
{

    if (&data_in != &data_out)
    {
        info =
            MPI_Reduce(data_in, data_out, size, M_getMPIType<T>(), MPI_SUM, receiver, communicator);
    }
    else
    {
        if (this->get_rank() == receiver)
        {
            info = MPI_Reduce(MPI_IN_PLACE,
                              data_in,
                              size,
                              M_getMPIType<T>(),
                              MPI_SUM,
                              receiver,
                              communicator);
        }
        else
        {
            info = MPI_Reduce(data_in,
                              data_out,
                              size,
                              M_getMPIType<T>(),
                              MPI_SUM,
                              receiver,
                              communicator);
        }
    }
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error("ERROR: MlMPI::reduceSum( T*& data_in, T*& data_out,\n"
                                     "                         Int size, Int receiver )\n"
                                     "MPI_Reduce()");
        }
    );
}

template<typename T>
void MlMPI::reduceProduct(T& data_in, T& data_out, Int receiver)
{
    // differing send and recive buffer
    if (&data_in != &data_out)
    {
        info =
            MPI_Reduce(&data_in, &data_out, 1, M_getMPIType<T>(), MPI_PROD, receiver, communicator);
    }
    // send and recive buffer are the same
    else
    {
        if (this->get_rank() == receiver)
        {
            info = MPI_Reduce(MPI_IN_PLACE,
                              &data_in,
                              1,
                              M_getMPIType<T>(),
                              MPI_PROD,
                              receiver,
                              communicator);
        }
        else
        {
            info = MPI_Reduce(&data_in,
                              &data_out,
                              1,
                              M_getMPIType<T>(),
                              MPI_PROD,
                              receiver,
                              communicator);
        }
    }
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error(
                "ERROR: MlMPI::reduceProduct( T& data_in, T& data_out, Int receiver )\n"
                "MPI_Reduce()");
        }
    );
}

template<typename T>
void MlMPI::reduceProduct(T*& data_in, T*& data_out, Int size, Int receiver)
{
    if (&data_in != &data_out)
    {
        info = MPI_Reduce(data_in,
                          data_out,
                          size,
                          M_getMPIType<T>(),
                          MPI_PROD,
                          receiver,
                          communicator);
    }
    else
    {
        if (this->get_rank() == receiver)
        {
            info = MPI_Reduce(MPI_IN_PLACE,
                              data_in,
                              size,
                              M_getMPIType<T>(),
                              MPI_PROD,
                              receiver,
                              communicator);
        }
        else
        {
            info = MPI_Reduce(data_in,
                              data_out,
                              size,
                              M_getMPIType<T>(),
                              MPI_PROD,
                              receiver,
                              communicator);
        }
    }
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error("ERROR: MlMPI::reduceProduct( T*& data_in, T*& data_out,\n"
                                     "                             Int size, Int receiver )\n"
                                     "MPI_Reduce()");
        }
    );
}

template<typename T>
void MlMPI::reduceMax(T& data_in, T& data_out, Int receiver)
{
    // differing send and recive buffer
    if (&data_in != &data_out)
    {
        info =
            MPI_Reduce(&data_in, &data_out, 1, M_getMPIType<T>(), MPI_MAX, receiver, communicator);
    }
    // send and recive buffer are the same
    else
    {
        if (this->get_rank() == receiver)
        {
            info = MPI_Reduce(MPI_IN_PLACE,
                              &data_in,
                              1,
                              M_getMPIType<T>(),
                              MPI_MAX,
                              receiver,
                              communicator);
        }
        else
        {
            info = MPI_Reduce(&data_in,
                              &data_out,
                              1,
                              M_getMPIType<T>(),
                              MPI_MAX,
                              receiver,
                              communicator);
        }
    }
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error(
                "ERROR: MlMPI::reduceMax( T& data_in, T& data_out, Int receiver )\n"
                "MPI_Reduce()");
        }
    );
}

template<typename T>
void MlMPI::reduceMax(T*& data_in, T*& data_out, Int size, Int receiver)
{
    if (&data_in != &data_out)
    {
        info =
            MPI_Reduce(data_in, data_out, size, M_getMPIType<T>(), MPI_MAX, receiver, communicator);
    }
    else
    {
        if (this->get_rank() == receiver)
        {
            info = MPI_Reduce(MPI_IN_PLACE,
                              data_in,
                              size,
                              M_getMPIType<T>(),
                              MPI_MAX,
                              receiver,
                              communicator);
        }
        else
        {
            info = MPI_Reduce(data_in,
                              data_out,
                              size,
                              M_getMPIType<T>(),
                              MPI_MAX,
                              receiver,
                              communicator);
        }
    }
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error("ERROR: MlMPI::reduceMax( T*& data_in, T*& data_out,\n"
                                     "                         Int size, Int receiver )\n"
                                     "MPI_Reduce()");
        }
    );
}

template<typename T>
void MlMPI::reduceMin(T& data_in, T& data_out, Int receiver)
{
    // differing send and recive buffer
    if (&data_in != &data_out)
    {
        info =
            MPI_Reduce(&data_in, &data_out, 1, M_getMPIType<T>(), MPI_MIN, receiver, communicator);
    }
    // send and recive buffer are the same
    else
    {
        if (this->get_rank() == receiver)
        {
            info = MPI_Reduce(MPI_IN_PLACE,
                              &data_in,
                              1,
                              M_getMPIType<T>(),
                              MPI_MIN,
                              receiver,
                              communicator);
        }
        else
        {
            info = MPI_Reduce(&data_in,
                              &data_out,
                              1,
                              M_getMPIType<T>(),
                              MPI_MIN,
                              receiver,
                              communicator);
        }
    }
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error(
                "ERROR: MlMPI::reduceMin( T& data_in, T& data_out, Int receiver )\n"
                "MPI_Reduce()");
        }
    );
}

template<typename T>
void MlMPI::reduceMin(T*& data_in, T*& data_out, Int size, Int receiver)
{
    if (&data_in != &data_out)
    {
        MPI_Reduce(data_in, data_out, size, M_getMPIType<T>(), MPI_MIN, receiver, communicator);
    }
    else
    {
        if (this->get_rank() == receiver)
        {
            MPI_Reduce(MPI_IN_PLACE,
                       data_in,
                       size,
                       M_getMPIType<T>(),
                       MPI_MIN,
                       receiver,
                       communicator);
        }
        else
        {
            MPI_Reduce(data_in, data_out, size, M_getMPIType<T>(), MPI_MIN, receiver, communicator);
        }
    }
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error("ERROR: MlMPI::reduceMax( T*& data_in, T*& data_out,\n"
                                     "                         Int size, Int receiver )\n"
                                     "MPI_Reduce()");
        }
    );
}

template<typename T>
void MlMPI::allReduceSum(T& data_send, T& data_recv)
{
    if (&data_send != &data_recv)
    {
        info = MPI_Allreduce(&data_send, &data_recv, 1, M_getMPIType<T>(), MPI_SUM, communicator);
    }
    else
    {
        info = MPI_Allreduce(MPI_IN_PLACE, &data_send, 1, M_getMPIType<T>(), MPI_SUM, communicator);
    }
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error("ERROR: MlMPI::allReduceSum( T& data_send, T& data_recv )\n"
                                     "MPI_Allreduce()");
        }
    );
}

template<typename T>
void MlMPI::allReduceSum(T*& data_send, T*& data_recv, Int size)
{
    if (&data_send != &data_recv)
    {
        info = MPI_Allreduce(data_send, data_recv, size, M_getMPIType<T>(), MPI_SUM, communicator);
    }
    else
    {
        info =
            MPI_Allreduce(MPI_IN_PLACE, data_send, size, M_getMPIType<T>(), MPI_SUM, communicator);
    }
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error("ERROR: MlMPI::allReduceSum( T*& data_send, T*& data_recv,\n"
                                     "                            Int size )\n"
                                     "MPI_Allreduce()");
        }
    );
}

template<typename T>
void MlMPI::allReduceProduct(T& data_send, T& data_recv)
{
    if (&data_send != &data_recv)
    {
        info = MPI_Allreduce(&data_send, &data_recv, 1, M_getMPIType<T>(), MPI_PROD, communicator);
    }
    else
    {
        info =
            MPI_Allreduce(MPI_IN_PLACE, &data_send, 1, M_getMPIType<T>(), MPI_PROD, communicator);
    }
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error(
                "ERROR: MlMPI::allReduceProduct( T& data_send, T& data_recv )\n"
                "MPI_Allreduce()");
        }
    );
}

template<typename T>
void MlMPI::allReduceProduct(T*& data_send, T*& data_recv, Int size)
{
    if (&data_send != &data_recv)
    {
        info = MPI_Allreduce(data_send, data_recv, size, M_getMPIType<T>(), MPI_PROD, communicator);
    }
    else
    {
        info =
            MPI_Allreduce(MPI_IN_PLACE, data_send, size, M_getMPIType<T>(), MPI_PROD, communicator);
    }
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error(
                "ERROR: MlMPI::allReduceProduct( T*& data_send, T*& data_recv,\n"
                "                                Int size )\n"
                "MPI_Allreduce()");
        }
    );
}

template<typename T>
void MlMPI::allReduceMax(T& data_send, T& data_recv)
{
    if (&data_send != &data_recv)
    {
        info = MPI_Allreduce(&data_send, &data_recv, 1, M_getMPIType<T>(), MPI_MAX, communicator);
    }
    else
    {
        info = MPI_Allreduce(MPI_IN_PLACE, &data_send, 1, M_getMPIType<T>(), MPI_MAX, communicator);
    }
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error("ERROR: MlMPI::allReduceMax( T& data_send, T& data_recv )\n"
                                     "MPI_Allreduce()");
        }
    );
}

template<typename T>
void MlMPI::allReduceMax(T*& data_send, T*& data_recv, Int size)
{
    if (&data_send != &data_recv)
    {
        info = MPI_Allreduce(data_send, data_recv, size, M_getMPIType<T>(), MPI_MAX, communicator);
    }
    else
    {
        info =
            MPI_Allreduce(MPI_IN_PLACE, data_send, size, M_getMPIType<T>(), MPI_MAX, communicator);
    }
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error("ERROR: MlMPI::allReduceMax( T*& data_send, T*& data_recv,\n"
                                     "                            Int size )\n"
                                     "MPI_Allreduce()");
        }
    );
}

template<typename T>
void MlMPI::allReduceMin(T& data_send, T& data_recv)
{
    if (&data_send != &data_recv)
    {
        info = MPI_Allreduce(&data_send, &data_recv, 1, M_getMPIType<T>(), MPI_MIN, communicator);
    }
    else
    {
        info = MPI_Allreduce(MPI_IN_PLACE, &data_send, 1, M_getMPIType<T>(), MPI_MIN, communicator);
    }
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error("ERROR: MlMPI::allReduceMin( T& data_send, T& data_recv )\n"
                                     "MPI_Allreduce()");
        }
    );
}

template<typename T>
void MlMPI::allReduceMin(T*& data_send, T*& data_recv, Int size)
{
    if (&data_send != &data_recv)
    {
        info = MPI_Allreduce(data_send, data_recv, size, M_getMPIType<T>(), MPI_MIN, communicator);
    }
    else
    {
        info =
            MPI_Allreduce(MPI_IN_PLACE, data_send, size, M_getMPIType<T>(), MPI_MIN, communicator);
    }
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error("ERROR: MlMPI::allReduceMin( T*& data_send, T*& data_recv,\n"
                                     "                            Int size )\n"
                                     "MPI_Allreduce()");
        }
    );
}

template<typename T>
void MlMPI::allocateShmemSpace(T*& data, std::size_t size)
{

    if (!shmemFileRead) { readMaximalAvailableSharedMemory(); }

    checkEnoughSharedMemory(sizeof(T), size);

    barrier();
    if (sharedCommunicator)
    {
#ifndef sysv
        MPI_Aint sizeType = sizeof(T);
        MPI_Aint byteSize = size * sizeType;
        //local unit size for displacements, in bytes (positive integer)
        Int disp_unit = 1;
        if (rank == 0)
        {
            MPI_Win_allocate_shared(byteSize, sizeType, MPI_INFO_NULL, communicator, &data, &win);
        }
        else
        {
            MPI_Win_allocate_shared(0, sizeType, MPI_INFO_NULL, communicator, &data, &win);
            MPI_Win_shared_query(win, 0, &sizeType, &disp_unit, &data);
        }
        sharedAlloc = true;
        barrier();
#else
        key_t key;
        Int   shmflg;
        Int   shmid;
        Int   istat = 0;
        key = IPC_PRIVATE;
        if (rank == 0)
        {
            shmflg = IPC_CREAT | IPC_EXCL | 0600 | SHM_NORESERVE;
            shmid = shmget(key, size * sizeof(T), shmflg);
            if (shmid == -1) istat++;
        }
        // summing istat errors
        allReduceSum(istat, istat);
        if (istat != 0)
        {
            std::cout << std::endl;
            throw std::runtime_error("Error in System V shared memory allocation");
        }
        // send shared memory id to all cores in communicator
        bcast(shmid, 0);
        // create void pointer which will be your data pointer
        void* r;
        shmflg = 0;
        // From doc page: If shmaddr is NULL, the system by default chooses the suitable
        // address to attach the segment.
        r = shmat(shmid, NULL, shmflg);
        // cast shared memory ptr to input shmem pointer
        data = static_cast<T*>(r);
        barrier();
        // mark for deletion which has to be done by root rank
        if (rank == 0)
        {
            struct shmid_ds buf;
            istat = shmctl(shmid, IPC_RMID, &buf);
            if (istat == -1)
            {
                std::cout << "Can not detach shared segment \n"
                          << __FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
                exit(1);
            }
        }
        barrier();
#endif
    }
    else
    {
        throw std::runtime_error("ERROR:MlMPI::allocateShmemSpace You are trying to\n"
                                 "allocate shared memory space on a non shared memory\n"
                                 "communicator. Stopping code\n");
    }
    barrier();
}

} //namespace vaspml
#endif
