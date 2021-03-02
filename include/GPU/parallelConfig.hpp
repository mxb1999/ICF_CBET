  #ifndef PARALLEL
  #define PARALLEL
    #include <cuda_runtime.h>
    #define CUDA 1
    #define OPENCL 2
    #include <string>
    #include <map>
    class DeviceDataEntry//store and modify pointers for device-host comms
    {
        private:
            void* hostAddress;//store the pointer to the allocated memory on the host system
            void* deviceAddress;//store the pointer to the allocated memory on the device
            double* refPointer;//reference pointer used to find device symbol
            size_t size;//store the size of the data being allocated
            std::string name;
            bool isP;
        public:
            DeviceDataEntry(void* host, void* device, void* ref, const std::string n, size_t s, bool p)
            {
                this->hostAddress = host;
                this->deviceAddress = device;
                this->refPointer = (double*)ref;
                this->name = n;
                this->size = s;
                this->isP = p;
            };
            int getSize()
            {
                return size;
            };
            bool isPointer()
            {
                return isP;
            };
            void* getHostAddress()
            {
                return hostAddress;
            }
            double* getRefAddress()
            {
                return refPointer;
            }
            std::string getName()
            {
                return name;
            }
            void* getDeviceAddress()
            {
                return deviceAddress;
            }
            void changeHostAddress(void* addr)
            {
                hostAddress = addr;
            }
            void changeDeviceAddress(void* addr)
            {
                deviceAddress = addr;
            }
    };
    extern void gpuInit();
    //cuda helper functions
    
    #define entryPair std::string, DeviceDataEntry*
    class GConfig//class for accessing and modifying GPU configuration info
    {
        private:
            short type;//type defined in above macros
            short cus;//number of sms/compute units
            short numProc;//number of cards in system
            char* arch;//architecture generation number, used for compilation and debugging
            std::map<entryPair>* addressHash;
        public:
            void gpuMalloc(void** p, size_t size);//wrapper function for GPU memory allocation
            void gpuMemcpy(void* device, void* host, size_t size, int direction);//wrapper function for GPU memcpy functionality
            void deviceDataTransfer(int dir);//wrapper function for GPU
            void transferConst();
            void addData(void* h, void* d, void* r, const std::string name, size_t size, bool isP);
            DeviceDataEntry* getDataEntry(const std::string name);
            int getType()
            {
                return this->type;
            }
            std::map<entryPair>* getMap()
            {
                return this->addressHash;
            }
            GConfig(int t);
        // void

    };
    extern void CGPUMemCpy(void* dest, void* source, size_t size, int direction);
    extern void CGPUAlloc(void* p, size_t size);
    extern void CGPUDataTrans(GConfig* device, int dir);//function to transfer large quantities of variables from host to GPU, makes bulk transfers easier
    
    
    
    //static int gputype;//stores whether the GPU is utilizing CUDA or OpenCL  //TO BE COMPLETED LATER
    #define NONE 0
    #define WAITING 1
    #define SENDING 2
    //MPI configuration struct, to be kept constant in each system, updated via IO functions
    #define IDLESTAT 0
    #define INITSTAT 1
    #define RAYTRACKSTAT 2
    #define FIELDSOLVESTAT 3
    #define CBETSTAT 4
    #define IOSTAT 5
    struct mpiconfig
    {
        int numProc;//number of processors being used in system
        char* stat;//status of each processor, defined in above macros
        char* modStat;//defines what module each processor is executing
    };
    //Input functions
    extern int importConfigFile();//import user selected parameters, either from GUI or previously exported config file

    //MPI IO handlers
    extern int mpiStatQuery(mpiconfig* mpi);//update mpiconfig struct
    extern void mpiPass(mpiconfig* mpi, int target, void* package);//MPI function to pass data contained in package to target
    extern void mpiSwap(mpiconfig* mpi, int target, void* package, void* receive);//MPI function to simultaneously pass data stored in package and recieve data into recieve
    extern void mpiRecieve(mpiconfig* mpi, int target, void* receive);//MPI function to recieve data expected from target into recieve



    /*
        GPU IO interface definition:
        1. Declare new GPU object-> Define the type (CUDA or OpenCL), Number of Compute units,
        number of cards in system, architecture ID
        2. Use GPUMemcpy and GPUAlloc to generalize interactions, ops depend on struct
        3. io_interface.cpp will take info from struct to determine which device specific functions to call
        4. 
    
    */
    //GPU IO handlers
    extern int GPUTrace();//GPU wrapper function for trace GPU function
    extern int GPUCBET();//GPU wrapper function for CBET GPU calls
#endif