# cuFFTAdvisor
Tool for generating optimal setting of the cuFFT library.

This tool is able to:
* benchmark some transformation.
* run fast heuristics which will recommend faster settings, using either padding or cropping. 
* autotune the heuristics proposals and thus obtainign the best setting possible

## Publication
_pending_

## Build
Simply run `make` in the root directory. Make sure that you have NVCC and Cuda libraries in the PATH

## Example
`./build/cuFFTAdvisor -help`

`./build/cuFFTAdvisor -benchmark -device 0 -x 85 -y 1`

`./build/cuFFTAdvisor -recommend 30 -x 2273 -y 2273 -n 6 --outOfPlaceOnly --realOnly --floatOnly --batchOnly  --forwardOnly --maxSignalInc 10`

`./build/cuFFTAdvisor -find 20 -x 445105 -n 5 --inPlaceOnly --complexOnly --doubleOnly --forwardOnly`

## Known limitations
Internal heuristics is based on Cuda v.8, so it might not be very accurate for other versions.
