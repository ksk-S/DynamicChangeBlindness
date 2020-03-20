################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../CSSR-v0.1.1/AllStates.cpp \
../CSSR-v0.1.1/G_Array.cpp \
../CSSR-v0.1.1/Hash.cpp \
../CSSR-v0.1.1/Hash2.cpp \
../CSSR-v0.1.1/Machine.cpp \
../CSSR-v0.1.1/Main.cpp \
../CSSR-v0.1.1/ParseTree.cpp \
../CSSR-v0.1.1/States.cpp \
../CSSR-v0.1.1/Test.cpp \
../CSSR-v0.1.1/TransTable.cpp 

OBJS += \
./CSSR-v0.1.1/AllStates.o \
./CSSR-v0.1.1/G_Array.o \
./CSSR-v0.1.1/Hash.o \
./CSSR-v0.1.1/Hash2.o \
./CSSR-v0.1.1/Machine.o \
./CSSR-v0.1.1/Main.o \
./CSSR-v0.1.1/ParseTree.o \
./CSSR-v0.1.1/States.o \
./CSSR-v0.1.1/Test.o \
./CSSR-v0.1.1/TransTable.o 

CPP_DEPS += \
./CSSR-v0.1.1/AllStates.d \
./CSSR-v0.1.1/G_Array.d \
./CSSR-v0.1.1/Hash.d \
./CSSR-v0.1.1/Hash2.d \
./CSSR-v0.1.1/Machine.d \
./CSSR-v0.1.1/Main.d \
./CSSR-v0.1.1/ParseTree.d \
./CSSR-v0.1.1/States.d \
./CSSR-v0.1.1/Test.d \
./CSSR-v0.1.1/TransTable.d 


# Each subdirectory must supply rules for building sources it contributes
CSSR-v0.1.1/%.o: ../CSSR-v0.1.1/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C++ Compiler'
	g++ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


