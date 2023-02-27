#pragma once
#include<burden.h>

struct Person
{
	double age;
	double risk;
	int home;
	Burden burden;
};

void increment_age(Person &person);