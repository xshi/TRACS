################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_UPPER_SRCS += \
../src/AutoPilot_eTCT.C \
../src/H1DConvolution.C \
../src/PlotCarrierFile.C \
../src/TMeasDict.C \
../src/TMeasHeaderDict.C \
../src/TWaveDict.C 

CPP_SRCS += \
../src/Carrier.cpp \
../src/CarrierCollection.cpp \
../src/CarrierMobility.cpp \
../src/CarrierTransport.cpp \
../src/DoTRACSFit.cpp \
../src/DoTracsOnly.cpp \
../src/Edge_tree.cpp \
../src/Global.cpp \
../src/MfgTRACSFit.cpp \
../src/SMSDSubDomains.cpp \
../src/SMSDetector.cpp \
../src/TMeas.cpp \
../src/TMeasHeader.cpp \
../src/TRACSFit.cpp \
../src/TRACSInterface.cpp \
../src/TScan.cpp \
../src/TWaveform.cpp \
../src/Threading.cpp \
../src/Utilities.cpp \
../src/programaPi.cpp 

C_UPPER_DEPS += \
./src/AutoPilot_eTCT.d \
./src/H1DConvolution.d \
./src/PlotCarrierFile.d \
./src/TMeasDict.d \
./src/TMeasHeaderDict.d \
./src/TWaveDict.d 

OBJS += \
./src/AutoPilot_eTCT.o \
./src/Carrier.o \
./src/CarrierCollection.o \
./src/CarrierMobility.o \
./src/CarrierTransport.o \
./src/DoTRACSFit.o \
./src/DoTracsOnly.o \
./src/Edge_tree.o \
./src/Global.o \
./src/H1DConvolution.o \
./src/MfgTRACSFit.o \
./src/PlotCarrierFile.o \
./src/SMSDSubDomains.o \
./src/SMSDetector.o \
./src/TMeas.o \
./src/TMeasDict.o \
./src/TMeasHeader.o \
./src/TMeasHeaderDict.o \
./src/TRACSFit.o \
./src/TRACSInterface.o \
./src/TScan.o \
./src/TWaveDict.o \
./src/TWaveform.o \
./src/Threading.o \
./src/Utilities.o \
./src/programaPi.o 

CPP_DEPS += \
./src/Carrier.d \
./src/CarrierCollection.d \
./src/CarrierMobility.d \
./src/CarrierTransport.d \
./src/DoTRACSFit.d \
./src/DoTracsOnly.d \
./src/Edge_tree.d \
./src/Global.d \
./src/MfgTRACSFit.d \
./src/SMSDSubDomains.d \
./src/SMSDetector.d \
./src/TMeas.d \
./src/TMeasHeader.d \
./src/TRACSFit.d \
./src/TRACSInterface.d \
./src/TScan.d \
./src/TWaveform.d \
./src/Threading.d \
./src/Utilities.d \
./src/programaPi.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.C
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D __cplusplus=201103L -I/usr/include/eigen3/ -I/usr/include/boost/thread/ -I/usr/include/ -I/usr/include/boost -I/home/jcalvopi/TRACS_Concurrency/include/ -I/usr/local/root/include/ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D __cplusplus=201103L -I/usr/include/eigen3/ -I/usr/include/boost/thread/ -I/usr/include/ -I/usr/include/boost -I/home/jcalvopi/TRACS_Concurrency/include/ -I/usr/local/root/include/ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


