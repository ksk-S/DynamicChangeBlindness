################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../decisional_states-1.0/examples/CellularAutomaton.cpp \
../decisional_states-1.0/examples/EvenProcess.cpp \
../decisional_states-1.0/examples/ImageFilter.cpp \
../decisional_states-1.0/examples/SymbolicSeries.cpp \
../decisional_states-1.0/examples/TimeSeries.cpp 

OBJS += \
./decisional_states-1.0/examples/CellularAutomaton.o \
./decisional_states-1.0/examples/EvenProcess.o \
./decisional_states-1.0/examples/ImageFilter.o \
./decisional_states-1.0/examples/SymbolicSeries.o \
./decisional_states-1.0/examples/TimeSeries.o 

CPP_DEPS += \
./decisional_states-1.0/examples/CellularAutomaton.d \
./decisional_states-1.0/examples/EvenProcess.d \
./decisional_states-1.0/examples/ImageFilter.d \
./decisional_states-1.0/examples/SymbolicSeries.d \
./decisional_states-1.0/examples/TimeSeries.d 


# Each subdirectory must supply rules for building sources it contributes
decisional_states-1.0/examples/%.o: ../decisional_states-1.0/examples/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C++ Compiler'
	g++ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


