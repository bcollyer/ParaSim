#ifndef PERSON_H
#define PERSON_H
#include<Burden.h>

struct Person
{
  double age;
  double risk;
  int home;
  Burden burden;
};

void increment_age(Person &person);


#endif
