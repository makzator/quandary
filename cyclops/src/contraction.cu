#include <stdlib.h>
#include <stdio.h>

#include <cuda_runtime.h>
#include <cutensor.h>

#include <unordered_map>
#include <vector>

// Handle cuTENSOR errors
#define HANDLE_ERROR(x) {                                                              \
  const auto err = x;                                                                  \
  if( x != CUTENSOR_STATUS_SUCCESS )                                                   \
  { printf("Error: %s in line %d\n", cutensorGetErrorString(x), __LINE__); exit(-1); } \
}

int main(int argc, char** argv)
{
  // Host element type definition
  typedef float floatTypeA;
  typedef float floatTypeB;
  typedef float floatTypeC;
  typedef float floatTypeCompute;

  // CUDA types
  cudaDataType_t typeA = CUDA_R_32F;
  cudaDataType_t typeB = CUDA_R_32F;
  cudaDataType_t typeC = CUDA_R_32F;
  cutensorComputeType_t typeCompute = CUTENSOR_R_MIN_32F;

  floatTypeCompute alpha = (floatTypeCompute)1.1f;
  floatTypeCompute beta  = (floatTypeCompute)0.9f;

  printf("Include headers and define data types\n");

  /* ***************************** */

  // Create vector of modes
  std::vector<int> modeC{'m','u','n','v'};
  std::vector<int> modeA{'m','h','k','n'};
  std::vector<int> modeB{'u','k','v','h'};
  int nmodeA = modeA.size();
  int nmodeB = modeB.size();
  int nmodeC = modeC.size();

  // Extents
  std::unordered_map<int, int64_t> extent;
  extent['m'] = 96;
  extent['n'] = 96;
  extent['u'] = 96;
  extent['v'] = 64;
  extent['h'] = 64;
  extent['k'] = 64;

  // Create a vector of extents for each tensor
  std::vector<int64_t> extentC;
  for(auto mode : modeC)
      extentC.push_back(extent[mode]);
  std::vector<int64_t> extentA;
  for(auto mode : modeA)
      extentA.push_back(extent[mode]);
  std::vector<int64_t> extentB;
  for(auto mode : modeB)
      extentB.push_back(extent[mode]);

  printf("Define modes and extents\n");

  /* ***************************** */

  // Number of elements of each tensor
  size_t elementsA = 1;
  for(auto mode : modeA)
      elementsA *= extent[mode];
  size_t elementsB = 1;
  for(auto mode : modeB)
      elementsB *= extent[mode];
  size_t elementsC = 1;
  for(auto mode : modeC)
      elementsC *= extent[mode];

  // Size in bytes
  size_t sizeA = sizeof(floatTypeA) * elementsA;
  size_t sizeB = sizeof(floatTypeB) * elementsB;
  size_t sizeC = sizeof(floatTypeC) * elementsC;

  // Allocate on device
  void *A_d, *B_d, *C_d;
  cudaMalloc((void**)&A_d, sizeA);
  cudaMalloc((void**)&B_d, sizeB);
  cudaMalloc((void**)&C_d, sizeC);

  // Allocate on host
  floatTypeA *A = (floatTypeA*) malloc(sizeof(floatTypeA) * elementsA);
  floatTypeB *B = (floatTypeB*) malloc(sizeof(floatTypeB) * elementsB);
  floatTypeC *C = (floatTypeC*) malloc(sizeof(floatTypeC) * elementsC);

  // Initialize data on host
  for(int64_t i = 0; i < elementsA; i++)
      A[i] = (((float) rand())/RAND_MAX - 0.5)*100;
  for(int64_t i = 0; i < elementsB; i++)
      B[i] = (((float) rand())/RAND_MAX - 0.5)*100;
  for(int64_t i = 0; i < elementsC; i++)
      C[i] = (((float) rand())/RAND_MAX - 0.5)*100;

  // Copy to device
  cudaMemcpy(C_d, C, sizeC, cudaMemcpyHostToDevice);
  cudaMemcpy(A_d, A, sizeA, cudaMemcpyHostToDevice);
  cudaMemcpy(B_d, B, sizeB, cudaMemcpyHostToDevice);

  printf("Allocate, initialize and transfer tensors\n");

  /* ***************************** */

  // Initialize cuTENSOR library
  cutensorHandle_t handle;
  cutensorInit(&handle);

  // Create Tensor Descriptors
  cutensorTensorDescriptor_t descA;
  HANDLE_ERROR( cutensorInitTensorDescriptor( &handle,
              &descA,
              nmodeA,
              extentA.data(),
              NULL,/*stride*/
              typeA, CUTENSOR_OP_IDENTITY ) );

  cutensorTensorDescriptor_t descB;
  HANDLE_ERROR( cutensorInitTensorDescriptor( &handle,
              &descB,
              nmodeB,
              extentB.data(),
              NULL,/*stride*/
              typeB, CUTENSOR_OP_IDENTITY ) );

  cutensorTensorDescriptor_t descC;
  HANDLE_ERROR( cutensorInitTensorDescriptor( &handle,
              &descC,
              nmodeC,
              extentC.data(),
              NULL,/*stride*/
              typeC, CUTENSOR_OP_IDENTITY ) );

  printf("Initialize cuTENSOR and tensor descriptors\n");

  /* ***************************** */

   //Retrieve the memory alignment for each tensor
   uint32_t alignmentRequirementA;
   HANDLE_ERROR( cutensorGetAlignmentRequirement( &handle,
              A_d,
              &descA,
              &alignmentRequirementA) );

   uint32_t alignmentRequirementB;
   HANDLE_ERROR( cutensorGetAlignmentRequirement( &handle,
              B_d,
              &descB,
              &alignmentRequirementB) );

   uint32_t alignmentRequirementC;
   HANDLE_ERROR( cutensorGetAlignmentRequirement( &handle,
              C_d,
              &descC,
              &alignmentRequirementC) );

  printf("Query best alignment requirement for our pointers\n");

  /* ***************************** */

  // Create the Contraction Descriptor
  cutensorContractionDescriptor_t desc;
  HANDLE_ERROR( cutensorInitContractionDescriptor( &handle,
              &desc,
              &descA, modeA.data(), alignmentRequirementA,
              &descB, modeB.data(), alignmentRequirementB,
              &descC, modeC.data(), alignmentRequirementC,
              &descC, modeC.data(), alignmentRequirementC,
              typeCompute) );

  printf("Initialize contraction descriptor\n");

  /* ***************************** */

  // Set the algorithm to use
  cutensorContractionFind_t find;
  HANDLE_ERROR( cutensorInitContractionFind(
              &handle, &find,
              CUTENSOR_ALGO_DEFAULT) );

  printf("Initialize settings to find algorithm\n");

  /* ***************************** */

  // Query workspace
  size_t worksize = 0;
  HANDLE_ERROR( cutensorContractionGetWorkspace(&handle,
              &desc,
              &find,
              CUTENSOR_WORKSPACE_RECOMMENDED, &worksize ) );

  // Allocate workspace
  void *work = nullptr;
  if(worksize > 0)
  {
      if( cudaSuccess != cudaMalloc(&work, worksize) ) // This is optional!
      {
          work = nullptr;
          worksize = 0;
      }
  }

  printf("Query recommended workspace size and allocate it\n");

  /* ***************************** */

  // Create Contraction Plan
  cutensorContractionPlan_t plan;
  HANDLE_ERROR( cutensorInitContractionPlan(&handle,
                                            &plan,
                                            &desc,
                                            &find,
                                            worksize) );

  printf("Create plan for contraction\n");

  /* ***************************** */

  cutensorStatus_t err;

  // Execute the tensor contraction
  err = cutensorContraction(&handle,
                            &plan,
                     (void*)&alpha, A_d,
                                    B_d,
                     (void*)&beta,  C_d,
                                    C_d,
                            work, worksize, 0 /* stream */);
  cudaDeviceSynchronize();

  // Check for errors
  if(err != CUTENSOR_STATUS_SUCCESS)
  {
      printf("ERROR: %s\n", cutensorGetErrorString(err));
  }

  printf("Execute contraction from plan\n");

  /* ***************************** */

  if ( A ) free( A );
  if ( B ) free( B );
  if ( C ) free( C );
  if ( A_d ) cudaFree( A_d );
  if ( B_d ) cudaFree( B_d );
  if ( C_d ) cudaFree( C_d );
  if ( work ) cudaFree( work );

  printf("Successful completion\n");

  return 0;
}