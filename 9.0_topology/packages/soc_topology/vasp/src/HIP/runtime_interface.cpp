#include <hip_runtime.h>

extern "C" {
    void hip_init(int device_num){
        hipInit(0);
        hipSetDevice(device_num);
    }
    void hip_device_get_uuid(char *id, int device_num){
        hipUUID_t uuid;
        hipDevice_t device = device_num;
        hipDeviceGetUuid(&uuid, device);
        id = uuid.bytes;
    }
    void hip_mem_get_info(size_t *free, size_t *total, int *err){
        hipError_t ierr = hipMemGetInfo(free, total);
        *err = ierr;
    }
}
