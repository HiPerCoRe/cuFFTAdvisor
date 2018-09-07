CC=g++
SRC_DIR := .
OBJ_DIR := obj
BUILD_DIR := build/
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))
LDLIBFLAGS := --cudart static -shared -link -Wno-deprecated-gpu-targets
LDFLAGS := -lcufft -Wno-deprecated-gpu-targets
CPPFLAGS := ...
CXXFLAGS := -O3 --use_fast_math -Xcompiler -fPIC -std=c++11 -Wno-deprecated-gpu-targets

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: libcuFFTAdvisor.so cuFFTAdvisor

# Tool invocations
libcuFFTAdvisor.so: $(OBJ_FILES)
	mkdir -p $(BUILD_DIR)
	nvcc $(LDLIBFLAGS) -o $(BUILD_DIR)$@ $^

cuFFTAdvisor: $(OBJ_FILES)
	mkdir -p $(BUILD_DIR)
	nvcc -o $(BUILD_DIR)$@ $^ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	mkdir -p $(@D)
	nvcc $(CXXFLAGS) -c -o $@ $<

# Other Targets
clean:
	rm -rf $(OBJ_DIR) $(BUILD_DIR)
