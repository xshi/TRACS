################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_UPPER_SRCS += \
../myApp/analysis_FIT/DoAll.C 

C_UPPER_DEPS += \
./myApp/analysis_FIT/DoAll.d 

OBJS += \
./myApp/analysis_FIT/DoAll.o 


# Each subdirectory must supply rules for building sources it contributes
myApp/analysis_FIT/%.o: ../myApp/analysis_FIT/%.C
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D __cplusplus=201103L -I/usr/include/eigen3/ -I/usr/include/boost/thread/ -I/usr/include/ -I/usr/include/boost -I/home/jcalvopi/TRACS_Concurrency/include/ -I/usr/local/root/include/ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


