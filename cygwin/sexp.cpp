#include <math.h>

double runif();
double sexp()
{
	return -log(runif());
}
