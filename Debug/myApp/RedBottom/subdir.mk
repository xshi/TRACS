################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_UPPER_SRCS += \
../myApp/RedBottom/CompareTRACS2Meas.C \
../myApp/RedBottom/Polot.C 

C_UPPER_DEPS += \
./myApp/RedBottom/CompareTRACS2Meas.d \
./myApp/RedBottom/Polot.d 

OBJS += \
./myApp/RedBottom/CompareTRACS2Meas.o \
./myApp/RedBottom/Polot.o 


# Each subdirectory must supply rules for building sources it contributes
myApp/RedBottom/%.o: ../myApp/RedBottom/%.C
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D __cplusplus=201103L -I/usr/include/eigen3/ -I/usr/include/boost/thread/ -I/usr/include/ -I/usr/include/boost -I/home/jcalvopi/TRACS_Concurrency/include/ -I/usr/local/root/include/ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


