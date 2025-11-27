#define CL_TARGET_OPENCL_VERSION 300 // optional, but silences that pragma
#include <CL/cl.h>
#include <stdio.h>
#include <stdlib.h>
#define N 100
#define CHECK(err, msg)                                                        \
  if (err != CL_SUCCESS) {                                                     \
    fprintf(stderr, "%s failed with error %d\n", msg, err);                    \
    exit(1);                                                                   \
  }

const char *kernelSoc = "__kernel void multiply_by_two(__global const "
                        "float* in, __global float* out) {\n"
                        "   int i = get_global_id(0);\n"
                        "    out[i] = 2.0f * in[i];\n"
                        "}\n";

int main(void) {
  cl_int err;
  float in[N], out[N];
  for (int i = 0; i < N; i++) {
    in[i] = (float)i;
  }
  cl_uint numPlatforms = 0;
  CHECK(clGetPlatformIDs(0, NULL, &numPlatforms), "clGetPlatformIDs count");
  if (numPlatforms == 0) {
    fprintf(stderr, "No OpenCL platforms found.\n");
    return 1;
  }
  printf("Platforms found: %d\n", numPlatforms);
  cl_platform_id *platforms =
      (cl_platform_id *)malloc(sizeof(cl_platform_id) * numPlatforms);

  CHECK(clGetPlatformIDs(numPlatforms, platforms, NULL),
        "clGetPlatformIDs get");

  cl_uint numDevices = 15;
  cl_int j = 0;
  cl_device_id *Device = NULL;
  for (j = 0; j < numPlatforms; j++) {

    err =
        clGetDeviceIDs(platforms[j], CL_DEVICE_TYPE_GPU, 0, NULL, &numDevices);
    if (err == CL_DEVICE_NOT_FOUND) {
      printf("Platform %d has no devices\n", j);
      printf("Platform %d has no devices, NumDevices = %d\n", j, numDevices);
      continue;
    };
    printf("Numdevices right before mallocing: %d\n", numDevices);
    cl_device_id *devices =
        (cl_device_id *)malloc(sizeof(cl_device_id) * numDevices);

    err = clGetDeviceIDs(platforms[j], CL_DEVICE_TYPE_GPU, numDevices, devices,
                         NULL);
    Device = devices;

    printf("Platform %u has %d devices.\n", j, numDevices);

    char device_name[128];
    clGetDeviceInfo(devices[0], CL_DEVICE_NAME, sizeof(device_name),
                    device_name, NULL);
    printf("OpenCL Device Name: %s\n", device_name);
    clGetDeviceInfo(devices[0], CL_DEVICE_VERSION, sizeof(device_name),
                    device_name, NULL);
    printf("OpenCL Device Version: %s\n", device_name);
    size_t max_size;
    clGetDeviceInfo(devices[0], CL_DEVICE_IMAGE_MAX_ARRAY_SIZE, sizeof(size_t),
                    &max_size, NULL);
    printf("OpenCL Device Array size: %zu\n", max_size);
    break;
  };
  printf("Numdevices after for loop: %d\n", numDevices);

  cl_context context =
      clCreateContext(NULL, numDevices, Device, NULL, NULL, &err);

  if (err != CL_SUCCESS) {
    printf("cl_context failed to create context\n");
    exit(1);
  };

#if CL_TARGET_OPENCL_VERSION >= 200
  cl_command_queue queue =
      clCreateCommandQueueWithProperties(context, *Device, NULL, &err);
#else
  cl_command_queue queue = clCreateCommandQueue(context, *Device, 0, &err);
#endif
  cl_mem bufferIn =
      clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * N, NULL, &err);
  CHECK(err, "bufferIn");
  cl_mem bufferOut =
      clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * N, NULL, &err);
  CHECK(err, "bufferOut");
  // clEnqueueWriteBufferRect For multidimensional data transfers

  err = clEnqueueWriteBuffer(queue, bufferIn, CL_TRUE, 0, sizeof(float) * N, in,
                             0, NULL, NULL);
  CHECK(err, "clEnqueueWriteBuffer");
}
