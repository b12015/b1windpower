#include <stdio.h>
#include <wiringPi.h>
int GetWind(void)
{
  int A=digitalRead(17);
  int B=digitalRead(18);
  int C=digitalRead(7);
  int D=digitalRead(22);
  int E=digitalRead(23);

  int windspeed=(A<<4)|(B<<3)|(C<<2)|(D<<1)|E;
  return windspeed;
}
