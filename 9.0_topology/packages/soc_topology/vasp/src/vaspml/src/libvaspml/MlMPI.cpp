#include "MlMPI.hpp"

#include <fstream>
#include <limits>
#include <sstream>

using namespace vaspml;

bool                   MlMPI::shmemFileRead = false;
unsigned long long int MlMPI::availableShmem = 0;

MlMPI::MlMPI(void)
{
    info = 0;
    numberRanks = 0;
    rank = 0;
    sharedAlloc = false;
    noCommunicator = true;
    sharedCommunicator = false;
    communicator = MPI_COMM_NULL;
    win = MPI_WIN_NULL;
}

void MlMPI::make_WorldComm(void)
{

    Int is_initialized = false;

    MPI_Initialized(&is_initialized);
    info = 1;
    if (!is_initialized) { info = MPI_Init(NULL, NULL); }
    else
    {
        throw std::runtime_error("ERROR: MlMPI::MlMPI( void ) \n"
                                 "You are trying to initialize a second MPI_COMM_WORLD\n");
    }

    if (info != 0)
    {
        throw std::runtime_error("ERROR in MlMPI::make_WorldComm(): MPI_Init failed");
    }

    communicator = MPI_COMM_WORLD;
    noCommunicator = false;
    set_numberRanks();
    set_rank();
    sharedAlloc = false;
    sharedCommunicator = false;
}

MlMPI::MlMPI(MPI_Comm communicator, bool sharedCommunicator)
{
    this->communicator = communicator;
    set_numberRanks();
    set_rank();
    this->sharedAlloc = false;
    noCommunicator = false;
    this->sharedCommunicator = sharedCommunicator;
}

MlMPI::MlMPI(const MlMPI& other) :
    info(other.info),
    numberRanks(other.numberRanks),
    rank(other.rank),
    communicator(other.communicator),
    win(other.win),
    sharedAlloc(other.sharedAlloc),
    sharedCommunicator(other.sharedCommunicator),
    noCommunicator(other.noCommunicator)
{}

MlMPI& MlMPI::operator=(const MlMPI& other)
{
    // preventing self assignement;
    if (this == &other) { return *this; }
    info = other.info;
    numberRanks = other.numberRanks;
    rank = other.rank;
    communicator = other.communicator;
    win = other.win;
    sharedAlloc = other.sharedAlloc;
    sharedCommunicator = other.sharedCommunicator;
    noCommunicator = other.noCommunicator;
    return *this;
}

MlMPI::MlMPI(MlMPI&& other) noexcept :
    info(other.info),
    numberRanks(other.numberRanks),
    rank(other.rank),
    communicator(other.communicator),
    win(other.win),
    sharedAlloc(other.sharedAlloc),
    sharedCommunicator(other.sharedCommunicator),
    noCommunicator(other.noCommunicator)
{
    other.info = 0;
    other.numberRanks = 0;
    other.rank = 0;
    other.communicator = MPI_COMM_NULL;
    other.win = MPI_WIN_NULL;
    other.sharedAlloc = false;
    other.sharedCommunicator = false;
    other.noCommunicator = true;
}

MlMPI& MlMPI::operator=(MlMPI&& other) noexcept
{

    if (this == &other) return *this;
    info = other.info;
    numberRanks = other.numberRanks;
    rank = other.rank;
    communicator = other.communicator;
    win = other.win;
    sharedAlloc = other.sharedAlloc;
    sharedCommunicator = other.sharedCommunicator;
    noCommunicator = other.noCommunicator;

    other.info = 0;
    other.numberRanks = 0;
    other.rank = 0;
    other.communicator = MPI_COMM_NULL;
    other.win = MPI_WIN_NULL;
    other.sharedAlloc = false;
    other.sharedCommunicator = false;
    other.noCommunicator = true;
    return *this;
}

MlMPI MlMPI::make_splitShared(void) const
{

    MPI_Comm comm_shmem;
    MPI_Comm_split_type(communicator, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_shmem);
    MlMPI mpi(comm_shmem, true);
    return mpi;
}

std::tuple<MlMPI, MlMPI> MlMPI::make_splitSharedInternode(void) const
{
    MlMPI    shared = make_splitShared();
    MPI_Comm internode_communicator;
    MPI_Comm_split(communicator,
                   shared.get_rank(),
                   shared.get_numberRanks(),
                   &internode_communicator);
    MlMPI internode = MlMPI(internode_communicator, false);
    return std::tie(shared, internode);
}

void MlMPI::set_numberRanks(void)
{
    info = MPI_Comm_size(communicator, &numberRanks);
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error("ERROR: MlMPI::set_numberRanks() MPI_Comm_size failed");
        }
    );
}

void MlMPI::set_rank(void)
{
    info = MPI_Comm_rank(communicator, &rank);
    VASPML_DEBUG_L1(
        if (info != 0)
        {
            throw std::runtime_error("ERROR: MlMPI::set_rank(): MPI_Comm_rank failed");
        }
    );
}

MlMPI::~MlMPI(void)
{
    // free window/shared memory if it was used and not freed
    if (sharedAlloc)
    {
        info = MPI_Win_free(&win);
        VASPML_DEBUG_L1(
            if (info != 0)
            {
                std::cout << "ERROR: MlMP::~MlMPI(): MPI_Win_free( &win ) failed" << std::endl;
            }
        );
    }
}

void MlMPI::freeWindow(void)
{
    if (sharedAlloc)
    {
        info = MPI_Win_free(&win);
        VASPML_DEBUG_L1(
            if (info != 0)
            {
                throw std::runtime_error("ERROR: MlMP::freeWindow() MPI_Comm_compare()");
            }
        );
    }
    sharedAlloc = false;
}

void MlMPI::fence(void)
{
    info = MPI_Win_fence(0, win);
    VASPML_DEBUG_L1(
        if (info != 0) { throw std::runtime_error("ERROR: MlMP::fence() MPI_Win_fence()"); }
    );
}

void MlMPI::barrier(void)
{
    info = MPI_Barrier(communicator);
    VASPML_DEBUG_L1(
        if (info != 0) { throw std::runtime_error("ERROR: MlMP::barrier() MPI_Barrier()"); }
    );
}

const Int& MlMPI::get_rank(void) const
{
    return rank;
}

const Int& MlMPI::get_numberRanks(void) const
{
    return numberRanks;
}

std::tuple<Int, Int> MlMPI::blockLoop(size_t n)
{
    Int iterations_per_process = n / numberRanks;
    Int remainder = n % numberRanks;
    Int start = rank * iterations_per_process + (rank < remainder ? rank : remainder);
    Int end = start + iterations_per_process + (rank < remainder ? 1 : 0);
    return std::tie(start, end);
}

bool MlMPI::get_sharedAllocCheck(void) const
{
    return sharedAlloc;
}

bool MlMPI::get_sharedCommunicatorCheck(void) const
{
    return sharedCommunicator;
}

bool MlMPI::get_noCommunicatorCheck(void) const
{
    return noCommunicator;
}

void MlMPI::readMaximalAvailableSharedMemory(void)
{

    String        fileName = "/proc/sys/kernel/shmmax";
    std::ifstream inStream;
    inStream.open("/proc/sys/kernel/shmmax");
    if (!inStream)
    {
        std::cout << "Warning: MlMPI::readMaximalAvailableSharedMemory " << std::endl;
        std::cout << "Could not open file " << fileName << std::endl;
        availableShmem = std::numeric_limits<unsigned long long int>::max();
    }
    else { inStream >> availableShmem; }
    shmemFileRead = true;
}

void MlMPI::checkEnoughSharedMemory(const std::size_t sizeType, const std::size_t numberElements)
{
    std::size_t size = sizeType * numberElements;
    if (availableShmem >= size) availableShmem -= size;
    else
    {
        std::cout << "WARNING: MlMPI::checkEnoughSharedMemory( const std::size_t sizeType, const "
                     "std::size_t numberElements )!\n Not enough shared memory available"
                  << std::endl;
        //      throw std::runtime_error( "ERROR: MlMPI::checkEnoughSharedMemory( const std::size_t
        //      sizeType, const std::size_t numberElements )\n Not enough shared memory available"
        //      );
    }
}

MPI_Comm MlMPI::get_communicator(void) const
{
    return communicator;
}
