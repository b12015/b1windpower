# b1windpower
# Created by Daniel Haverås and Adam Richert in April 2015

This repository contains files generated from Simulink. Together with the file Finalmodel_main.c and Getwind.c, they create a
wind power plant simulation model for the Raspberry Pi. An Enercon E44 is simulated by this program. The file named ModelCopy2.c contains the program functions. ModelCopy2.h is a header file which links the necessary files from the project into the program.

Getwind is a C function that uses the wiringPi library. It reads data from 5 GPIO pins A, B, C, D and E to retrieve a 5-bit number between 0 and 31. This is used as input for the program.

To compile and run the program, follow these steps:
1.  Install a Linux OS on a RPi, patched with a real-time kernel. For this project, we used Raspbian with the RT_PREEMPT patch
2.  Install the wiringPi library. For more info on this, see http://www.wiringpi.com
3.  Connect an A/D-converter to GPIO pin 17,18,7,22,23 on the RPi. Here, GPIO Pin 17 represent bit 4 (A), adding 16 to the input, while pin 23 represents bit 0 (E), adding 1 to the input. To make the A/D-converter, an Arduino was used for this project.
4.  Aquire superuser access. This is necessary for several aspects of the program, for example wiringPi.
5.  Place all files from this repository in a common folder on the RPi
6.  Compile the program with the command >> gcc Finalmodel_main.c -DRT -o Windpower -lrt -lm -lwiringPi
7.  Run the program by writing >> sudo ./Windpower
